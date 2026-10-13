"""CellOracle graph metrics and candidates extracted from celloracle_candidates.ipynb.

Python API: run_celloracle(adata, output_dir="celloracle_results")
CLI: python celloracle_candidates.py --input data.h5ad --output-dir results
The input AnnData is not modified. Exactly four CSV tables are written.
Dependencies: anndata, numpy, pandas, scipy, celloracle (and its dependencies).
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

LOGGER = logging.getLogger(__name__)
CELL_TYPES = {"HSC": "HSC", "GMP": "GMP", "CD14_Mono": "CD14 Mono"}


def rank_candidates(
    metrics: pd.DataFrame,
    celltype: str,
    *,
    metric: str = "eigenvector_centrality",
    min_jmml_rank: float = 0.90,
    min_delta_rank: float | None = None,
) -> pd.DataFrame:
    """Match the final notebook blocks: union of genes, missing scores -> zero.

    min_delta_rank=None preserves the notebook's DISABLED delta-rank filter.
    These thresholds are exploratory, not statistical significance tests.
    """
    required = {"gene", "cluster", metric}
    missing = required.difference(metrics.columns)
    if missing:
        raise ValueError(f"Missing metric columns: {sorted(missing)}")
    if not 0 <= min_jmml_rank <= 1:
        raise ValueError("min_jmml_rank must be between 0 and 1.")
    control, jmml = f"{celltype}__Control", f"{celltype}__JMML"
    absent = {control, jmml}.difference(metrics["cluster"].dropna())
    if absent:
        raise ValueError(f"Missing network groups: {sorted(absent)}")
    selected = metrics.loc[metrics["cluster"].isin([control, jmml]),
                           ["gene", "cluster", metric]]
    if selected.duplicated(["gene", "cluster"]).any():
        raise ValueError("Duplicate gene/cluster pairs in network scores.")
    comparison = (selected.pivot(index="gene", columns="cluster", values=metric)
                  .reindex(columns=[control, jmml])
                  .rename(columns={control: "Control", jmml: "JMML"}))
    comparison.columns.name = None
    comparison["Control_was_missing"] = comparison["Control"].isna()
    comparison["JMML_was_missing"] = comparison["JMML"].isna()
    comparison[["Control", "JMML"]] = comparison[["Control", "JMML"]].fillna(0)
    comparison["rank_Control"] = comparison["Control"].rank(method="average", pct=True)
    comparison["rank_JMML"] = comparison["JMML"].rank(method="average", pct=True)
    comparison["delta_score"] = comparison["JMML"] - comparison["Control"]
    comparison["delta_rank"] = comparison["rank_JMML"] - comparison["rank_Control"]
    mask = ((comparison["JMML"] > 0)
            & (comparison["rank_JMML"] >= min_jmml_rank)
            & (comparison["delta_score"] > 0))
    if min_delta_rank is not None:
        mask &= comparison["delta_rank"] >= min_delta_rank
    comparison["selected_by_rank"] = mask
    return (comparison.loc[mask]
            .sort_values(["delta_rank", "JMML"], ascending=False, kind="mergesort")
            .reset_index())


def _prepare_adata(adata: Any, celltype_col: str, condition_col: str,
                   counts_layer: str, hvg_col: str, embedding_name: str) -> Any:
    import anndata as ad

    if not isinstance(adata, ad.AnnData):
        raise TypeError("adata must be an anndata.AnnData object.")
    for col in (celltype_col, condition_col):
        if col not in adata.obs:
            raise ValueError(f"Missing adata.obs[{col!r}].")
    if counts_layer not in adata.layers:
        raise ValueError(f"Missing adata.layers[{counts_layer!r}] containing raw counts.")
    if hvg_col not in adata.var:
        raise ValueError(f"Missing adata.var[{hvg_col!r}].")
    if embedding_name not in adata.obsm:
        raise ValueError(f"Missing adata.obsm[{embedding_name!r}].")
    hvg = adata.var[hvg_col]
    if not hvg.dropna().isin([True, False]).all():
        raise ValueError("HVG column must contain boolean values (or 0/1).")
    mask = (adata.obs[celltype_col].isin(CELL_TYPES.values())
            & adata.obs[condition_col].isin(["Control", "JMML"]))
    selected = adata[mask.to_numpy(), hvg.fillna(False).astype(bool).to_numpy()]
    if selected.n_obs < 3 or selected.n_vars < 2:
        raise ValueError("Too few selected cells or highly variable genes.")
    if not selected.var_names.is_unique:
        raise ValueError("Gene names must be unique before running CellOracle.")
    group_counts = selected.obs.groupby([celltype_col, condition_col], observed=True).size()
    for celltype in CELL_TYPES.values():
        for condition in ("Control", "JMML"):
            count = int(group_counts.get((celltype, condition), 0))
            if count < 2:
                raise ValueError(f"Group {celltype}__{condition} has {count} cells; need at least 2.")
    counts = selected.layers[counts_layer]
    values = counts.data if sparse.issparse(counts) else np.asarray(counts)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Counts must be finite and nonnegative.")
    if not np.allclose(values, np.rint(values)):
        raise ValueError("Counts layer must contain unnormalized integer-valued counts.")
    LOGGER.info("Selected %d cells and %d HVGs; dense counts require approximately %.2f GiB",
                selected.n_obs, selected.n_vars,
                selected.n_obs * selected.n_vars * counts.dtype.itemsize / 1024**3)
    dense = counts.toarray() if sparse.issparse(counts) else np.array(counts, copy=True)
    result = ad.AnnData(X=dense, obs=selected.obs.copy(), var=selected.var.copy(),
                       obsm={embedding_name: selected.obsm[embedding_name].copy()})
    result.obs_names_make_unique()
    result.obs[celltype_col] = result.obs[celltype_col].astype("category")
    result.obs[celltype_col] = result.obs[celltype_col].cat.remove_unused_categories()
    return result


def _select_pca_dims(ratios: np.ndarray, requested: int | None) -> int:
    available = len(ratios)
    if available < 1:
        raise ValueError("PCA returned no components.")
    if requested is not None:
        if not 1 <= requested <= available:
            raise ValueError(f"n_pca_dims must be between 1 and {available}.")
        return requested
    # Preserve the original heuristic, but handle missing/zero cutoffs safely.
    hits = np.where(np.diff(np.diff(np.cumsum(ratios)) > 0.002))[0]
    if len(hits) and hits[0] > 0:
        return min(int(hits[0]), 50, available)
    LOGGER.warning("Notebook PCA cutoff not usable; falling back to min(50, available PCs).")
    return min(50, available)


def run_celloracle(
    adata: Any,
    output_dir: str | Path = "celloracle_results",
    *,
    celltype_col: str = "celltype_l2_merged",
    condition_col: str = "sample_groups_wide",
    counts_layer: str = "counts",
    hvg_col: str = "highly_variable",
    embedding_name: str = "X_umap",
    base_grn: pd.DataFrame | None = None,
    n_pca_dims: int | None = None,
    k: int | None = None,
    n_jobs: int = 4,
    alpha: float = 10,
    link_p: float = 0.001,
    threshold_number: int = 2000,
    metric: str = "eigenvector_centrality",
    min_jmml_rank: float = 0.90,
    min_delta_rank: float | None = None,
    verbose_level: int = 10,
    overwrite: bool = False,
) -> dict[str, pd.DataFrame]:
    """Run the notebook's graph/candidate workflow, save and return four tables.

    Return keys are output CSV filenames; values are DataFrames (gene is a column).
    No plots, Oracle/Links checkpoints, base-GRN CSVs or knockout simulations.
    Default base GRN is CellOracle's human promoter GRN; input gene identifiers
    must match that GRN. Pass base_grn to supply a preloaded/custom TF matrix.
    Existing target files are protected unless overwrite=True.
    """
    if n_jobs < 1 or alpha <= 0 or not 0 < link_p <= 1 or threshold_number < 1:
        raise ValueError("Invalid n_jobs, alpha, link_p or threshold_number.")
    if not 0 <= min_jmml_rank <= 1:
        raise ValueError("min_jmml_rank must be between 0 and 1.")
    directory = Path(output_dir)
    names = ["celloracle_graph_metrics.csv"] + [f"{key}_rank_candidates.csv" for key in CELL_TYPES]
    existing = [str(directory / name) for name in names if (directory / name).exists()]
    if existing and not overwrite:
        raise FileExistsError(f"Output files already exist: {existing}; use overwrite=True.")
    prepared = _prepare_adata(adata, celltype_col, condition_col, counts_layer,
                              hvg_col, embedding_name)
    import celloracle as co

    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(adata=prepared,
                                      cluster_column_name=celltype_col,
                                      embedding_name=embedding_name)
    oracle.import_TF_data(TF_info_matrix=(co.data.load_human_promoter_base_GRN()
                                         if base_grn is None else base_grn))
    oracle.perform_PCA()
    dims = _select_pca_dims(oracle.pca.explained_variance_ratio_, n_pca_dims)
    n_cells = oracle.adata.n_obs
    neighbors = max(1, int(0.025 * n_cells)) if k is None else k
    if not isinstance(neighbors, (int, np.integer)) or not 1 <= neighbors < n_cells:
        raise ValueError(f"k must be an integer between 1 and {n_cells - 1}.")
    LOGGER.info("KNN imputation: %d PCs, k=%d", dims, neighbors)
    oracle.knn_imputation(n_pca_dims=dims, k=neighbors, balanced=True,
                          b_sight=neighbors * 8, b_maxl=neighbors * 4, n_jobs=n_jobs)
    group_col = "celltype_condition"
    if group_col in (celltype_col, condition_col):
        raise ValueError("Metadata column names must differ from 'celltype_condition'.")
    obs = oracle.adata.obs
    obs[group_col] = (obs[celltype_col].astype(str) + "__"
                      + obs[condition_col].astype(str)).astype("category")
    links = oracle.get_links(cluster_name_for_GRN_unit=group_col,
                             alpha=alpha, verbose_level=verbose_level)
    links.filter_links(p=link_p, weight="coef_abs", threshold_number=threshold_number)
    links.get_network_score()
    metrics = links.merged_score.copy()
    if "gene" not in metrics.columns:
        metrics = metrics.rename_axis("gene").reset_index()
    tables = {"celloracle_graph_metrics.csv": metrics}
    for prefix, celltype in CELL_TYPES.items():
        tables[f"{prefix}_rank_candidates.csv"] = rank_candidates(
            metrics, celltype, metric=metric, min_jmml_rank=min_jmml_rank,
            min_delta_rank=min_delta_rank)
    directory.mkdir(parents=True, exist_ok=True)
    for name, table in tables.items():
        table.to_csv(directory / name, index=False)
        LOGGER.info("Saved %s (%d rows)", directory / name, len(table))
    return tables


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input AnnData .h5ad file")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory path, supplied by the workflow/config")
    parser.add_argument("--celltype-col", default="celltype_l2_merged")
    parser.add_argument("--condition-col", default="sample_groups_wide")
    parser.add_argument("--counts-layer", default="counts")
    parser.add_argument("--hvg-col", default="highly_variable")
    parser.add_argument("--embedding-name", default="X_umap")
    parser.add_argument("--n-pca-dims", type=int, default=None)
    parser.add_argument("--k", type=int, default=None)
    parser.add_argument("--n-jobs", type=int, default=4)
    parser.add_argument("--min-jmml-rank", type=float, default=0.90)
    parser.add_argument("--min-delta-rank", type=float, default=None)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    import anndata as ad
    adata = ad.read_h5ad(args.input)
    options = vars(args).copy()
    options.pop("input")
    run_celloracle(adata, **options)


if __name__ == "__main__":
    main()
