"""Merge DEG with precomputed DepMap features and prioritize selective dependencies.

Example:
    python merge_deg_depmap.py --deg cd14_mono_deg.csv --depmap depmap_features.csv \
        --out cd14_mono_depmap.csv
"""

import argparse
from pathlib import Path
import pandas as pd


def read_table(path):
    return pd.read_csv(path, sep='\t' if Path(path).suffix.lower() in ('.tsv', '.txt') else ',')


def main():
    parser = argparse.ArgumentParser(description='Annotate DEG with precomputed DepMap features')
    parser.add_argument('--deg', required=True, help='DEG CSV/TSV file')
    parser.add_argument('--depmap', required=True, help='Precomputed DepMap features CSV/TSV file')
    parser.add_argument('--gene-col', default='gene', help='Gene column in DEG (default: gene)')
    parser.add_argument('--out', required=True, help='Output merged CSV path, supplied by the workflow/config')
    parser.add_argument('--strength', type=float, default=0.5, help='Minimum DepMap strength')
    parser.add_argument('--prevalence', type=float, default=0.5, help='Minimum dependent fraction')
    parser.add_argument('--selectivity', type=float, default=0.0, help='Minimum AML-vs-other selectivity')
    args = parser.parse_args()

    deg = read_table(args.deg)
    depmap = read_table(args.depmap)
    if args.gene_col not in deg.columns:
        deg.columns=['gene', 'baseMean', 'log2FC', 'lfcSE', 'stat', 'pvalue', 'padj']

    required = ['gene', 'depmap_strength', 'depmap_fraction_dependent', 'depmap_selectivity', 'common_essential']
    missing = [col for col in required if col not in depmap.columns]
    if missing:
        raise ValueError(f'DepMap table is missing: {missing}')

    deg = deg.copy()
    depmap = depmap.copy()
    deg['gene'] = deg[args.gene_col].astype('string').str.strip().str.replace(r'\s*\(\d+\)$', '', regex=True)
    depmap['gene'] = depmap['gene'].astype('string').str.strip().str.replace(r'\s*\(\d+\)$', '', regex=True)
    # Preserve one gene column even when the input gene column had another name.
    if depmap['gene'].duplicated().any():
        duplicates = depmap.loc[depmap['gene'].duplicated(keep=False), 'gene'].head(10).tolist()
        raise ValueError(f'Duplicate genes in DepMap after cleaning: {duplicates}')

    # Robustly parse boolean flags written by pandas/other CSV exporters.
    common = depmap['common_essential'].astype('string').str.strip().str.lower()
    valid = common.isin(['true', 'false', '1', '0']) | common.isna()
    if (~valid).any():
        raise ValueError(f"Unknown common_essential values: {common.loc[~valid].unique().tolist()}")
    depmap['common_essential'] = common.map({'true': True, '1': True, 'false': False, '0': False}).astype('boolean')

    merged = deg.merge(depmap, on='gene', how='left', validate='many_to_one', indicator=True)
    matched = int(merged['_merge'].eq('both').sum())
    merged = merged.drop(columns='_merge')

    for col in ['depmap_strength', 'depmap_fraction_dependent', 'depmap_selectivity']:
        merged[col] = pd.to_numeric(merged[col], errors='coerce')

    candidates = merged.loc[
        merged['common_essential'].eq(False).fillna(False)
        & merged['depmap_strength'].gt(args.strength)
        & merged['depmap_fraction_dependent'].gt(args.prevalence)
        & merged['depmap_selectivity'].gt(args.selectivity)
    ].copy()
    candidates = candidates.sort_values(
        ['depmap_selectivity', 'depmap_strength'], ascending=[False, False]
    )
    merged = merged.sort_values(['depmap_selectivity', 'depmap_strength'], ascending=[False, False])
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    candidates_path = out.with_name(out.stem + '_candidates.csv')
    print(merged[[
    'gene',
    'depmap_selectivity',
    'depmap_strength']].head(15).to_string(index=False))
    merged.to_csv(out, index=False)
    candidates.to_csv(candidates_path, index=False)

    print(f'DEG rows: {len(deg)}; matched with DepMap: {matched}; candidates: {len(candidates)}')
    print(f'Merged: {out}\nCandidates: {candidates_path}')
    cols = [c for c in ['gene', 'log2FC', 'padj', 'depmap_strength', 'depmap_fraction_dependent', 'depmap_selectivity', 'common_essential'] if c in candidates.columns]
    print('\nTop 10 candidates (ranked by selectivity, then strength):')
    print(candidates[cols].head(10).to_string(index=False))


if __name__ == '__main__':
    main()
