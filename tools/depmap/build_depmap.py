#!/usr/bin/env python3
"""Build DepMap gene-level dependency features for selected cancer models.

Examples:
  python build_depmap_features.py --data-dir tools/depmap --list-diseases
  python build_depmap_features.py --data-dir tools/depmap \
      --diseases 'Acute Myeloid Leukemia' 'Myeloproliferative Neoplasms' \
      --out depmap_features.csv
  python build_depmap_features.py --data-dir tools/depmap \
      --field OncotreeSubtype --diseases 'AML, NOS' 'Acute Monoblastic/Monocytic Leukemia' \
      --out aml_subtypes_depmap.csv

The background for selectivity comprises all screened models not selected as targets.
Positive selectivity means stronger dependency in selected models than background.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

MODEL_FIELDS = ("OncotreePrimaryDisease", "OncotreeSubtype", "OncotreeLineage")


def read_genes(path):
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str).str.strip()
    df.columns = df.columns.astype(str).str.strip().str.replace(r"\s+\(\d+\)$", "", regex=True)
    if df.columns.duplicated().any():
        duplicate = df.columns[df.columns.duplicated()].unique().tolist()
        raise ValueError(f"Duplicate gene symbols after removing Entrez IDs in {path}: {duplicate[:10]}")
    if df.index.duplicated().any():
        raise ValueError(f"Duplicate ModelID values in {path}")
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", type=Path, required=True,
                    help="Directory containing DepMap CSV inputs (provided from config by the workflow)")
    ap.add_argument("--field", choices=MODEL_FIELDS, default="OncotreePrimaryDisease",
                    help="Model.csv field used to select the target models")
    ap.add_argument("--diseases", nargs="+", help="Exact values in the selected Model.csv field")
    ap.add_argument("--subtypes", nargs="+", help="Exact values in OncotreeSubtype; combined with --diseases using OR")
    ap.add_argument("--list-diseases", action="store_true", help="Print available values and model counts, then exit")
    ap.add_argument("--out", type=Path, required=True,
                    help="Output feature table path (provided from config by the workflow)")
    ap.add_argument("--min-models", type=int, default=1,
                    help="Minimum models with measured Gene Effect for a gene to retain its metrics")
    args = ap.parse_args()

    models = pd.read_csv(args.data_dir / "Model.csv")
    if "ModelID" not in models.columns or args.field not in models.columns:
        raise ValueError(f"Model.csv must contain ModelID and {args.field}")
    models["ModelID"] = models["ModelID"].astype(str).str.strip()

    if args.list_diseases:
        for field in MODEL_FIELDS:
            if field in models:
                print(f"\n=== {field} ===")
                print(models[field].value_counts(dropna=False).to_string())
        return

    if not args.diseases and not args.subtypes:
        ap.error("Specify --diseases and/or --subtypes, or use --list-diseases")

    selected_mask = pd.Series(False, index=models.index)
    for field, values in ((args.field, args.diseases), ("OncotreeSubtype", args.subtypes)):
        if not values:
            continue
        if field not in models.columns:
            raise ValueError(f"Model.csv has no column {field}")
        available = set(models[field].dropna().astype(str))
        missing = sorted(set(values) - available)
        if missing:
            raise ValueError(f"Unknown values for {field}: {missing}. Run --list-diseases")
        selected_mask |= models[field].isin(values)

    selected = models.loc[selected_mask].drop_duplicates(subset="ModelID").copy()

    effect = read_genes(args.data_dir / "CRISPRGeneEffect.csv")
    dependency = read_genes(args.data_dir / "CRISPRGeneDependency.csv")
    common = pd.read_csv(args.data_dir / "CRISPRInferredCommonEssentials.csv")

    shared_genes = effect.columns.intersection(dependency.columns, sort=False)
    shared_ids = effect.index.intersection(dependency.index, sort=False)
    selected_ids = shared_ids.intersection(selected.ModelID, sort=False)
    other_ids = shared_ids.difference(selected_ids, sort=False)
    if len(selected_ids) == 0:
        raise ValueError("No selected models appear in both CRISPR matrices; check disease selection")
    if len(other_ids) == 0:
        raise ValueError("No background models available for selectivity")

    E = effect.loc[selected_ids, shared_genes].apply(pd.to_numeric, errors="coerce")
    P = dependency.loc[selected_ids, shared_genes].apply(pd.to_numeric, errors="coerce")
    E_other = effect.loc[other_ids, shared_genes].apply(pd.to_numeric, errors="coerce")

    n_effect = E.notna().sum(axis=0)
    n_dependency = P.notna().sum(axis=0)
    mean_effect = E.mean(axis=0).where(n_effect >= args.min_models)
    median_effect = E.median(axis=0).where(n_effect >= args.min_models)
    other_mean = E_other.mean(axis=0)
    prevalence = (P.gt(0.5).sum(axis=0) / n_dependency.replace(0, np.nan))
    prevalence = prevalence.where(n_dependency >= args.min_models)
    selectivity = (other_mean - mean_effect).where(E_other.notna().sum(axis=0) >= args.min_models)

    # Most DepMap releases ship a single-column list of common essential genes.
    common_genes = set(common.iloc[:, 0].dropna().astype(str).str.strip()
                       .str.replace(r"\s+\(\d+\)$", "", regex=True))

    features = pd.DataFrame({
        "gene": shared_genes,
        "depmap_mean_effect": mean_effect.to_numpy(),
        "depmap_median_effect": median_effect.to_numpy(),
        "depmap_strength": (-mean_effect).clip(lower=0).to_numpy(),
        "depmap_fraction_dependent": prevalence.to_numpy(),
        "depmap_n_models": n_effect.to_numpy(),
        "depmap_n_dependency_models": n_dependency.to_numpy(),
        "depmap_effect_sd": E.std(axis=0).where(n_effect >= args.min_models).to_numpy(),
        "depmap_other_mean_effect": other_mean.to_numpy(),
        "depmap_selectivity": selectivity.to_numpy(),
    })
    features["common_essential"] = features.gene.isin(common_genes)
    features["selective_dependency"] = (
        (features.depmap_selectivity > 0.3)
        & (features.depmap_fraction_dependent > 0.5)
    )
    features = features.sort_values("depmap_strength", ascending=False, na_position="last")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    features.to_csv(args.out, index=False)
    selected_file = args.out.with_name(args.out.stem + "_models.csv")
    selected.loc[selected.ModelID.isin(selected_ids)].to_csv(selected_file, index=False)

    print(f"Selected field: {args.field}")
    print(f"Selected diseases: {args.diseases or []}")
    print(f"Selected subtypes: {args.subtypes or []}")
    print("Selection rule: OR (union), duplicates removed")
    print(f"Models in Model.csv: {len(selected)}")
    print(f"Screened selected models: {len(selected_ids)}")
    print(f"Screened background models: {len(other_ids)}")
    print(f"Shared genes: {len(shared_genes)}")
    print(f"Common essentials matched: {features.common_essential.sum()}")
    print(f"Features: {args.out}")
    print(f"Selected model metadata: {selected_file}")
    print("\nSelected models by primary disease:")
    print(selected.loc[selected.ModelID.isin(selected_ids), "OncotreePrimaryDisease"].value_counts().to_string())


if __name__ == "__main__":
    try:
        main()
    except (ValueError, FileNotFoundError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)