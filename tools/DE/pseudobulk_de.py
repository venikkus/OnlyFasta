#!/usr/bin/env python3
"""Pseudobulk JMML vs Control: HSC, GMP and CD14 Mono.

Based on data_look_and_psbulk_DE-Copy1.ipynb (PyDESeq2 0.5.4).
Python API: tables = run_pseudobulk_de(adata, outdir="DEGs")
CLI: python pseudobulk_de.py --input atlas.h5ad --outdir DEGs
Returns three DEG DataFrames and writes exactly three CSV tables.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

GROUP_MAP = {
    "PTPN11": "JMML", "NF1": "JMML", "NRAS": "JMML",
    "KRAS": "JMML", "CBL": "JMML", "RALD": "RALD",
    "Control": "Control",
}
CELLTYPES = {"HSC": "HSC", "GMP": "GMP", "CD14_Mono": "CD14 Mono"}


def run_pseudobulk_de(
    adata, outdir="DEGs", *, sample_col="batch",
    celltype_col="celltype_l2_merged", group_col="sample_groups",
    counts_layer="counts", min_cells=10, min_total_counts=1000,
    min_replicates=2, min_gene_count=10, fdr=0.05,
    lfc_cutoff=1.0, n_cpus=4,
):
    """Aggregate raw counts per sample and cell type; fit ~ condition.

    sample_col must identify independent biological replicates, not technical
    batches. group_col contains mutation labels in GROUP_MAP, or JMML/Control.
    adata is not modified. Returned dictionary keys: HSC, GMP, CD14_Mono.
    Only significant genes are exported; empty DEG tables retain column headers.
    Validation/model failures raise errors rather than silently skipping types.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    from pydeseq2.default_inference import DefaultInference

    if min_cells < 1 or min_total_counts < 1 or min_replicates < 2:
        raise ValueError("min_cells/counts >= 1; min_replicates >= 2 required")
    if min_gene_count < 1 or not 0 < fdr < 1 or lfc_cutoff < 0 or n_cpus < 1:
        raise ValueError("Invalid gene count, FDR, LFC cutoff or CPU count")
    missing = set([sample_col, celltype_col, group_col]) - set(adata.obs.columns)
    if missing:
        raise KeyError(f"Missing adata.obs columns: {sorted(missing)}")
    if counts_layer not in adata.layers:
        raise KeyError(f"Missing raw counts layer: {counts_layer!r}")
    if not adata.var_names.is_unique or adata.var_names.hasnans:
        raise ValueError("Gene identifiers must be unique and non-missing")

    obs = adata.obs[[sample_col, celltype_col, group_col]].copy()
    mapping = {**GROUP_MAP, "JMML": "JMML"}
    obs["condition"] = obs[group_col].astype(object).map(mapping)
    target = obs[celltype_col].isin(CELLTYPES.values())
    relevant = target & obs["condition"].isin(["JMML", "Control"])
    if obs.loc[relevant, sample_col].isna().any():
        raise ValueError("Missing sample identifiers in selected cells")
    consistency = (obs.loc[obs[sample_col].notna()]
                   .groupby(sample_col, observed=True)[group_col].nunique())
    if consistency.gt(1).any():
        raise ValueError(f"Samples assigned multiple groups: {consistency[consistency.gt(1)].index.tolist()}")
    excluded = obs.loc[target & ~relevant, group_col].value_counts(dropna=False)
    if len(excluded):
        print(f"Excluded groups (not JMML/Control): {excluded.to_dict()}")

    prepared = {}
    for key, label in CELLTYPES.items():
        positions = np.flatnonzero((relevant & obs[celltype_col].eq(label)).to_numpy())
        if not len(positions):
            raise ValueError(f"{label}: no JMML/Control cells")
        cells = obs.iloc[positions].copy()
        codes, samples = pd.factorize(cells[sample_col].astype(str), sort=True)
        index = pd.Index(samples, name=sample_col)
        X = adata.layers[counts_layer][positions, :]
        values = X.data if sparse.issparse(X) else np.asarray(X)
        if (not np.isfinite(values).all() or (values < 0).any()
                or not np.equal(values, np.rint(values)).all()):
            raise ValueError(f"{label}: counts must be finite nonnegative integers")
        X = X.astype(np.int64)
        aggregation = sparse.csr_matrix(
            (np.ones(len(codes), dtype=np.int64), (codes, np.arange(len(codes)))),
            shape=(len(samples), len(codes)),
        )
        summed = aggregation @ X
        summed = summed.toarray() if sparse.issparse(summed) else np.asarray(summed)
        counts = pd.DataFrame(summed, index=index, columns=adata.var_names)
        meta = pd.DataFrame({
            "condition": [cells.loc[codes == i, "condition"].iloc[0]
                          for i in range(len(samples))],
            "n_cells": np.bincount(codes, minlength=len(samples)),
            "total_counts": counts.sum(axis=1).to_numpy(),
        }, index=index)
        keep = meta.n_cells.ge(min_cells) & meta.total_counts.ge(min_total_counts)
        meta, counts = meta.loc[keep].copy(), counts.loc[keep].copy()
        n = meta.condition.value_counts().reindex(["Control", "JMML"], fill_value=0)
        if n.min() < min_replicates:
            raise ValueError(f"{label}: insufficient replicates after QC: {n.to_dict()}")
        counts = counts.loc[:, counts.ge(min_gene_count).sum(axis=0).ge(int(n.min()))]
        if counts.shape[1] == 0 or counts.sum(axis=1).eq(0).any():
            raise ValueError(f"{label}: empty gene set or sample after gene filtering")
        meta = meta[["condition"]]
        meta["condition"] = pd.Categorical(meta.condition, categories=["Control", "JMML"])
        prepared[key] = (counts, meta)
        print(f"{label}: replicates={n.to_dict()}, genes={counts.shape[1]}")

    results = {}
    for key, (counts, meta) in prepared.items():
        inference = DefaultInference(n_cpus=n_cpus)
        dds = DeseqDataSet(counts=counts, metadata=meta, design="~ condition",
                           refit_cooks=True, inference=inference)
        dds.deseq2()
        stats = DeseqStats(dds, contrast=["condition", "JMML", "Control"],
                          alpha=fdr, inference=inference)
        stats.summary()
        res = stats.results_df
        degs = res.loc[res.padj.lt(fdr) & res.log2FoldChange.abs().ge(lfc_cutoff)].copy()
        degs = degs.sort_values(["padj", "pvalue"])
        degs.index.name = "gene"
        results[key] = degs
        print(f"{key}: {len(degs)} DEGs")

    destination = Path(outdir)
    destination.mkdir(parents=True, exist_ok=True)
    for key, table in results.items():
        table.to_csv(destination / f"{key}_DEGs_JMML_vs_Control.csv")
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input AnnData .h5ad")
    parser.add_argument("--outdir", default="DEGs")
    parser.add_argument("--sample-col", default="batch")
    parser.add_argument("--celltype-col", default="celltype_l2_merged")
    parser.add_argument("--group-col", default="sample_groups")
    parser.add_argument("--counts-layer", default="counts")
    parser.add_argument("--min-cells", type=int, default=10)
    parser.add_argument("--min-total-counts", type=int, default=1000)
    parser.add_argument("--min-replicates", type=int, default=2)
    parser.add_argument("--min-gene-count", type=int, default=10)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--lfc-cutoff", type=float, default=1.0)
    parser.add_argument("--n-cpus", type=int, default=4)
    args = vars(parser.parse_args())
    import anndata as ad
    adata = ad.read_h5ad(args.pop("input"))
    run_pseudobulk_de(adata, **args)


if __name__ == "__main__":
    main()
