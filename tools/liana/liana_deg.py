#!/usr/bin/env python3
"""Compare LIANA CellChat scores in JMML vs Control for DEG receptors."""
import argparse
import logging
from pathlib import Path

import anndata as ad
import liana as li
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

LOG = logging.getLogger("liana_deg_diff")
KEYS = ["source", "target", "ligand", "receptor"]
SCORES = ["lr_probs", "cellchat_pvals"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True)
    p.add_argument("--deg", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--gene-col", default="gene")
    p.add_argument("--deg-celltype-col", default=None)
    p.add_argument("--condition-col", default="sample_groups_wide")
    p.add_argument("--celltype-col", default="celltype_l2_merged_rename")
    p.add_argument("--patient-col", default="batch")
    p.add_argument("--jmml-label", default="jmml")
    p.add_argument("--control-label", default="control")
    p.add_argument("--receivers", nargs="+", default=["GMP"])
    p.add_argument("--min-cells", type=int, default=10)
    p.add_argument("--expr-prop", type=float, default=0.1)
    p.add_argument("--n-perms", type=int, default=100)
    p.add_argument("--n-jobs", type=int, default=4)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--resource", default="consensus")
    p.add_argument("--max-patients-per-condition", type=int, default=None,
                   help="Testing only: select patients with the most cells")
    return p.parse_args()


def read_deg(path, gene_col, celltype_col):
    df = pd.read_csv(path, sep=None, engine="python")
    # Some DESeq2 exports put gene symbols in the unnamed first column.
    if gene_col not in df.columns:
        first = str(df.columns[0])
        if first.startswith("Unnamed:") or first in {"X", "index", "row.names"}:
            LOG.warning("Treating first DEG column %r as %r", first, gene_col)
            df = df.rename(columns={df.columns[0]: gene_col})
        else:
            raise ValueError(f"Missing DEG gene column {gene_col!r}; found {list(df.columns)}")
    if celltype_col and celltype_col not in df.columns:
        raise ValueError(f"Missing DEG cell type column {celltype_col!r}")
    cols = [gene_col] + ([celltype_col] if celltype_col else [])
    df = df[cols].dropna().copy()
    for col in cols:
        df[col] = df[col].astype(str).str.strip()
    df = df[df[gene_col].ne("")]
    return df.drop_duplicates()


def run_condition(backed, obs, label, args, outdir):
    conditions = obs[args.condition_col].astype("string").str.strip().str.lower()
    mask = (conditions.eq(label.strip().lower()).fillna(False)
            & obs[args.celltype_col].notna())

    if args.max_patients_per_condition is not None:
        if args.max_patients_per_condition < 1:
            raise ValueError("--max-patients-per-condition must be >= 1")
        if args.patient_col not in obs:
            raise ValueError(f"Missing patient column {args.patient_col!r}")
        selected = (obs.loc[mask, args.patient_col]
                    .value_counts().head(args.max_patients_per_condition).index)
        mask &= obs[args.patient_col].isin(selected)
        LOG.info("%s selected patients: %s", label, list(selected))

    idx = np.flatnonzero(mask.to_numpy(dtype=bool))
    if len(idx) == 0:
        raise ValueError(f"No cells found for condition {label!r}")
    LOG.info("%s: reading counts for %d cells", label, len(idx))
    counts = backed.layers["counts"][idx, :]
    counts = (counts.tocsr() if sp.issparse(counts) else sp.csr_matrix(counts))
    counts = counts.astype(np.float32)

    meta = obs.iloc[idx][[args.celltype_col]].copy()
    meta.index = pd.Index([f"cell_{i}" for i in idx])
    a = ad.AnnData(X=counts, obs=meta,
                   var=pd.DataFrame(index=backed.var_names.copy()))
    a.var_names_make_unique()
    a.obs[args.celltype_col] = a.obs[args.celltype_col].astype(str)
    sizes = a.obs[args.celltype_col].value_counts()
    valid = sizes[sizes >= args.min_cells].index
    a = a[a.obs[args.celltype_col].isin(valid)].copy()
    LOG.info("%s: %d cells, %d cell types", label, a.n_obs, len(valid))
    if not set(args.receivers).intersection(valid):
        raise ValueError(f"Requested receivers not found for {label}: {args.receivers}")

    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    li.mt.cellchat(
        a, groupby=args.celltype_col, resource_name=args.resource,
        expr_prop=args.expr_prop, min_cells=args.min_cells,
        n_perms=args.n_perms, n_jobs=args.n_jobs,
        use_raw=False, verbose=True, seed=args.seed,
    )
    lr = a.uns["liana_res"].copy()
    required = KEYS + SCORES
    missing = set(required) - set(lr.columns)
    if missing:
        raise ValueError(f"LIANA returned no columns: {missing}")
    lr = lr.loc[lr["target"].isin(args.receivers), required].copy()

    dup = lr.loc[lr.duplicated(subset=KEYS, keep=False)]
    if not dup.empty:
        dup.to_csv(outdir / f"{label}_duplicate_LR_rows.csv", index=False)
        conflicts = dup.groupby(KEYS, dropna=False)[SCORES].nunique(dropna=False)
        if conflicts.gt(1).any().any():
            raise ValueError(f"{label}: duplicate keys have different CellChat scores")
        LOG.warning("%s: dropping %d repeated LR rows with identical scores",
                    label, lr.duplicated(subset=KEYS).sum())
        lr = lr.drop_duplicates(subset=KEYS).copy()

    lr.to_csv(outdir / f"{label}_incoming_all.csv", index=False)
    return lr


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    deg = read_deg(args.deg, args.gene_col, args.deg_celltype_col)
    LOG.info("Loaded %d DEG records, %d unique genes",
             len(deg), deg[args.gene_col].nunique())

    backed = sc.read_h5ad(args.input, backed="r")
    try:
        obs = backed.obs.copy()
        for col in (args.condition_col, args.celltype_col):
            if col not in obs:
                raise ValueError(f"Missing obs column: {col}")
        if "counts" not in backed.layers:
            raise ValueError('Input AnnData has no layers["counts"]')
        jmml = run_condition(backed, obs, args.jmml_label, args, outdir)
        control = run_condition(backed, obs, args.control_label, args, outdir)
    finally:
        backed.file.close()

    jmml = jmml.rename(columns={c: f"{c}_jmml" for c in SCORES})
    control = control.rename(columns={c: f"{c}_control" for c in SCORES})
    result = jmml.merge(control, on=KEYS, how="outer", indicator=True,
                        validate="one_to_one")
    result["delta_lr_probs"] = (
        result["lr_probs_jmml"] - result["lr_probs_control"]
    )
    result["p_value"] = np.nan  # No JMML-vs-Control test has been performed.

    if args.deg_celltype_col:
        matched = deg.rename(columns={args.gene_col: "receptor",
                                      args.deg_celltype_col: "target"})
        result = result.merge(matched, on=["receptor", "target"],
                              how="inner", validate="many_to_one")
    else:
        deg_genes = set(deg[args.gene_col])
        result = result.loc[result["receptor"].isin(deg_genes)].copy()

    result = result.rename(columns={"receptor": "DEG_as_receptor"})
    cols = ["source", "target", "ligand", "DEG_as_receptor",
            "lr_probs_jmml", "lr_probs_control", "delta_lr_probs",
            "cellchat_pvals_jmml", "cellchat_pvals_control", "p_value", "_merge"]
    result = result[cols].sort_values("delta_lr_probs", ascending=False,
                                     na_position="last")
    result.to_csv(outdir / "JMML_vs_Control_DEG_receptors.csv", index=False)
    comparable = result.dropna(subset=["delta_lr_probs"])
    comparable.to_csv(outdir / "JMML_vs_Control_DEG_receptors_comparable.csv",
                      index=False)
    LOG.info("Saved %d DEG LR rows; %d comparable", len(result), len(comparable))
    LOG.info("Absent LR in one condition remains NaN; p_value is not a differential test")


if __name__ == "__main__":
    main()