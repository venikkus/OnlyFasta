#!/usr/bin/env python3
"""Merge an existing DEG table with precomputed DepMap features.

Example:
    python merge_deg_depmap.py --deg cd14_mono_deg.csv --depmap depmap_features.csv
"""

import argparse
from pathlib import Path
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Annotate DEG with precomputed DepMap features")
    parser.add_argument("--deg", required=True, help="DEG CSV/TSV file")
    parser.add_argument("--depmap", required=True, help="Precomputed DepMap features CSV/TSV file")
    parser.add_argument("--gene-col", default="gene", help="Gene column in DEG (default: gene)")
    parser.add_argument("--out", default="deg_depmap.csv", help="Output merged CSV")
    parser.add_argument("--strength", type=float, default=0.5, help="Minimum DepMap strength")
    parser.add_argument("--prevalence", type=float, default=0.5, help="Minimum dependent fraction")
    parser.add_argument("--selectivity", type=float, default=0.0, help="Minimum AML-vs-other selectivity")
    args = parser.parse_args()

    def read_table(path):
        return pd.read_csv(path, sep="\t" if Path(path).suffix.lower() in (".tsv", ".txt") else ",")

    deg = read_table(args.deg)
    depmap = read_table(args.depmap)
    deg.columns=['gene', 'baseMean', 'log2FC', 'lfcSE', 'stat', 'pvalue', 'padj']
    if args.gene_col not in deg.columns:
        raise ValueError(f"DEG has no column '{args.gene_col}'. Available: {list(deg.columns)}")
    required = ["gene", "depmap_strength", "depmap_fraction_dependent", "depmap_selectivity"]
    missing = [col for col in required if col not in depmap.columns]
    if missing:
        raise ValueError(f"DepMap table is missing: {missing}")

    deg = deg.copy()
    depmap = depmap.copy()
    deg["gene"] = deg[args.gene_col].astype("string").str.strip().str.replace(r"\s*\(\d+\)$", "", regex=True)
    depmap["gene"] = depmap["gene"].astype("string").str.strip().str.replace(r"\s*\(\d+\)$", "", regex=True)

    if depmap["gene"].duplicated().any():
        examples = depmap.loc[depmap["gene"].duplicated(), "gene"].head().tolist()
        raise ValueError(f"Duplicate genes in DepMap after cleaning: {examples}")

    merged = deg.merge(depmap, on="gene", how="left", validate="many_to_one", indicator=True)
    matched = int((merged["_merge"] == "both").sum())
    merged = merged.drop(columns="_merge")

    for col in required[1:]:
        merged[col] = pd.to_numeric(merged[col], errors="coerce")
    candidates = merged.loc[
        (merged["depmap_strength"] > args.strength)
        & (merged["depmap_fraction_dependent"] > args.prevalence)
        & (merged["depmap_selectivity"] > args.selectivity)
    ].copy()
    candidates = candidates.sort_values("depmap_strength", ascending=False)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    candidates_path = out.with_name(out.stem + "_candidates.csv")
    merged.to_csv(out, index=False)
    candidates.to_csv(candidates_path, index=False)

    print(f"DEG rows: {len(deg)}; matched with DepMap: {matched}; candidates: {len(candidates)}")
    print(f"Merged: {out}\nCandidates: {candidates_path}")
    cols = [c for c in ["gene", "depmap_strength", "depmap_fraction_dependent", "depmap_selectivity", "common_essential"] if c in candidates.columns]
    print("\nTop 10 candidates:")
    print(candidates[cols].head(10).to_string(index=False))


if __name__ == "__main__":
    main()