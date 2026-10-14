#!/usr/bin/env python3
"""Config-driven pseudobulk differential expression with PyDESeq2.

Run with the project YAML, for example:
    python tools/differential_expression/pseudobulk_de.py --config config/config.yaml

`pseudobulk_de` in the YAML supplies metadata columns, group and cell-type maps,
contrast labels, filters, and output location. The AnnData input is read-only.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from lazy_h5ad import open_counts  # noqa: E402


def _validate_mapping(name, mapping):
    if not isinstance(mapping, dict) or not mapping:
        raise ValueError(f"{name} must be a non-empty mapping")
    normalized = {str(k): str(v) for k, v in mapping.items()}
    if any(not k or not v for k, v in normalized.items()):
        raise ValueError(f"{name} contains an empty key or value")
    return normalized


def _validate_cell_types(cell_types):
    if not isinstance(cell_types, dict) or not cell_types:
        raise ValueError("cell_types must be a non-empty mapping")
    result = {}
    for key, values in cell_types.items():
        key = str(key)
        labels = values if isinstance(values, (list, tuple)) else [values]
        labels = tuple(str(value) for value in labels)
        if not key or not labels or any(not label for label in labels):
            raise ValueError("cell_types contains an empty output name or annotation value")
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", key):
            raise ValueError("cell_types keys may contain only letters, digits, dot, underscore and hyphen")
        result[key] = labels
    return result


def run_pseudobulk_de(
    adata,
    outdir,
    *,
    sample_col,
    celltype_col,
    group_col,
    block_col=None,
    counts_layer,
    group_map,
    cell_types,
    case_label,
    control_label,
    min_cells=10,
    min_total_counts=1000,
    min_replicates=2,
    min_gene_count=10,
    fdr=0.05,
    lfc_cutoff=1.0,
    n_cpus=4,
):
    """Aggregate raw counts per pseudobulk/population and test case vs control.

    `group_map` maps observed values in `group_col` to normalized condition
    labels; `cell_types` maps output names to values in `celltype_col`.
    Unmapped groups are excluded. Each `sample_col` ID must represent one
    pseudobulk and map consistently to one condition unless `block_col` is set.
    With a block column, counts are separately aggregated for each sample ×
    condition and the same patient can contribute both sides of the contrast.
    Set `block_col` to a
    patient/donor ID for paired designs (`~ subject + condition`). Results are
    written as `<population>_DEGs.csv`.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    from pydeseq2.default_inference import DefaultInference

    group_map = _validate_mapping("group_map", group_map)
    cell_types = _validate_cell_types(cell_types)
    if not case_label or not control_label or case_label == control_label:
        raise ValueError("case_label and control_label must be distinct non-empty strings")
    if not {case_label, control_label} <= set(group_map.values()):
        raise ValueError("group_map must map at least one observed group to each contrast label")
    if min_cells < 1 or min_total_counts < 1 or min_replicates < 2:
        raise ValueError("min_cells/counts >= 1; min_replicates >= 2 required")
    if min_gene_count < 1 or not 0 < fdr < 1 or lfc_cutoff < 0 or n_cpus < 1:
        raise ValueError("Invalid gene count, FDR, LFC cutoff or CPU count")
    required_obs = {sample_col, celltype_col, group_col}
    if block_col:
        required_obs.add(block_col)
    missing = required_obs - set(adata.obs.columns)
    if missing:
        raise KeyError(f"Missing adata.obs columns: {sorted(missing)}")
    if counts_layer != "X" and counts_layer not in adata.layers:
        raise KeyError(f"Missing raw counts layer: {counts_layer!r}")
    if not adata.var_names.is_unique or adata.var_names.hasnans:
        raise ValueError("Gene identifiers must be unique and non-missing")

    # A dataset may use the same obs field for cell selection and contrast
    # labels (e.g. sarcoma CellType_subset); avoid duplicate column names.
    obs_columns = list(dict.fromkeys([sample_col, celltype_col, group_col, *([block_col] if block_col else [])]))
    obs = adata.obs[obs_columns].copy()
    obs["condition"] = obs[group_col].astype(object).map(group_map)
    selected_labels = {label for labels in cell_types.values() for label in labels}
    target = obs[celltype_col].isin(selected_labels)
    relevant = target & obs["condition"].isin([case_label, control_label])
    if obs.loc[relevant, sample_col].isna().any():
        raise ValueError("Missing sample identifiers in selected cells")
    if block_col and obs.loc[relevant, block_col].isna().any():
        raise ValueError(f"Missing block/patient identifiers in adata.obs[{block_col!r}]")
    valid_obs = obs.loc[relevant & obs[sample_col].notna()]
    sample_conditions = valid_obs.groupby(sample_col, observed=True)["condition"].nunique()
    if not block_col and sample_conditions.gt(1).any():
        bad = sample_conditions[sample_conditions.gt(1)].index.tolist()
        raise ValueError(f"Samples assigned to multiple contrast groups; configure block_col for paired data: {bad}")
    if block_col:
        block_consistency = (valid_obs.groupby(sample_col, observed=True)[block_col].nunique())
        if block_consistency.gt(1).any():
            bad = block_consistency[block_consistency.gt(1)].index.tolist()
            raise ValueError(f"Pseudobulk sample IDs map to multiple patient/block IDs: {bad}")
    excluded = obs.loc[target & ~relevant, group_col].value_counts(dropna=False)
    if len(excluded):
        print(f"Excluded groups (not in contrast): {excluded.to_dict()}")

    prepared = {}
    for key, labels in cell_types.items():
        positions = np.flatnonzero((relevant & obs[celltype_col].isin(labels)).to_numpy())
        if not len(positions):
            raise ValueError(f"{key} ({', '.join(labels)}): no selected case/control cells")
        cells = obs.iloc[positions].copy()
        cell_groups = pd.MultiIndex.from_arrays(
            [cells[sample_col].astype(str), cells["condition"].astype(str)],
            names=[sample_col, "condition"])
        codes, sample_conditions_unique = cell_groups.factorize(sort=True)
        index = pd.Index([f"pseudobulk_{i}" for i in range(len(sample_conditions_unique))],
                         name="pseudobulk_id")
        matrix = adata.X if counts_layer == "X" else adata.layers[counts_layer]
        X = matrix[positions, :]
        values = X.data if sparse.issparse(X) else np.asarray(X)
        if (not np.isfinite(values).all() or (values < 0).any()
                or not np.equal(values, np.rint(values)).all()):
            raise ValueError(f"{key}: counts must be finite nonnegative integers")
        X = X.astype(np.int64)
        aggregation = sparse.csr_matrix(
            (np.ones(len(codes), dtype=np.int64), (codes, np.arange(len(codes)))),
            # One row per unique sample × condition pseudobulk. `samples` is
            # not defined in this branch; use the MultiIndex factorization
            # result that also defines the counts and metadata row order.
            shape=(len(sample_conditions_unique), len(codes)),
        )
        summed = aggregation @ X
        summed = summed.toarray() if sparse.issparse(summed) else np.asarray(summed)
        counts = pd.DataFrame(summed, index=index, columns=adata.var_names)
        meta = pd.DataFrame({
            # MultiIndex.factorize currently drops level names; condition is
            # the second component of the sample × condition key.
            "condition": sample_conditions_unique.get_level_values(1).to_list(),
            "n_cells": np.bincount(codes, minlength=len(sample_conditions_unique)),
            "total_counts": counts.sum(axis=1).to_numpy(),
        }, index=index)
        if block_col:
            meta["block"] = [cells.loc[codes == i, block_col].iloc[0]
                             for i in range(len(sample_conditions_unique))]
        keep = meta.n_cells.ge(min_cells) & meta.total_counts.ge(min_total_counts)
        meta, counts = meta.loc[keep].copy(), counts.loc[keep].copy()
        if block_col:
            complete_blocks = meta.groupby("block", observed=True)["condition"].agg(
                lambda values: set(values) >= {case_label, control_label})
            keep_blocks = complete_blocks.index[complete_blocks]
            meta = meta.loc[meta.block.isin(keep_blocks)].copy()
            counts = counts.loc[meta.index].copy()
            if meta.empty:
                raise ValueError(f"{key}: no complete paired blocks after sample QC")
        n = meta.condition.value_counts().reindex([control_label, case_label], fill_value=0)
        if n.min() < min_replicates:
            raise ValueError(f"{key}: insufficient replicates after QC: {n.to_dict()}")
        counts = counts.loc[:, counts.ge(min_gene_count).sum(axis=0).ge(int(n.min()))]
        if counts.shape[1] == 0 or counts.sum(axis=1).eq(0).any():
            raise ValueError(f"{key}: empty gene set or sample after gene filtering")
        if block_col:
            meta = meta[["condition", "block"]].rename(columns={"block": "subject"})
            meta["subject"] = meta["subject"].astype("category")
        else:
            meta = meta[["condition"]]
        meta["condition"] = pd.Categorical(meta.condition, categories=[control_label, case_label])
        prepared[key] = (counts, meta)
        print(f"{key}: replicates={n.to_dict()}, genes={counts.shape[1]}")

    results = {}
    for key, (counts, meta) in prepared.items():
        inference = DefaultInference(n_cpus=n_cpus)
        design = "~ subject + condition" if block_col else "~ condition"
        dds = DeseqDataSet(counts=counts, metadata=meta, design=design,
                           refit_cooks=True, inference=inference)
        dds.deseq2()
        stats = DeseqStats(dds, contrast=["condition", case_label, control_label],
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
        table.to_csv(destination / f"{key}_DEGs.csv")
    return results


def run_from_h5ad(path, outdir, *, counts_layer, **options):
    """run_pseudobulk_de on an .h5ad file, reading the counts matrix row-wise from disk."""
    handle, adata = open_counts(path, counts_layer)
    try:
        return run_pseudobulk_de(adata, outdir, counts_layer=counts_layer, **options)
    finally:
        handle.close()


def _resolve_project_path(project_root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (project_root / path).resolve()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Project YAML containing adata and pseudobulk_de sections")
    parser.add_argument("--input", help="Override adata.path from the config")
    parser.add_argument("--outdir", help="Override pseudobulk_de.output_dir from the config")
    args = parser.parse_args()

    import yaml
    config_path = Path(args.config).expanduser().resolve()
    config = yaml.safe_load(config_path.read_text())
    project_root = Path(config_path.parent.parent)
    paths = config.get("paths", {})
    configured_root = Path(paths.get("project_root", ".")).expanduser()
    project_root = (configured_root if configured_root.is_absolute()
                    else project_root / configured_root).resolve()
    adata_cfg = config["adata"]
    de_cfg = config["pseudobulk_de"]
    input_path = args.input or adata_cfg["path"]
    output_path = args.outdir or de_cfg.get("output_dir")
    if output_path is None:
        output_path = Path(paths.get("results_root", "results/snakemake")) / config["run_name"] / "pseudobulk_de"
    input_path = _resolve_project_path(project_root, str(input_path))
    output_path = _resolve_project_path(project_root, str(output_path))

    options = {k: de_cfg[k] for k in (
        "sample_col", "celltype_col", "group_col", "block_col", "counts_layer", "group_map",
        "cell_types", "case_label", "control_label", "min_cells",
        "min_total_counts", "min_replicates", "min_gene_count", "fdr",
        "lfc_cutoff", "n_cpus") if k in de_cfg}
    run_from_h5ad(input_path, output_path, **options)


if __name__ == "__main__":
    main()
