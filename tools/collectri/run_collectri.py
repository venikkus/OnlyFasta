#!/usr/bin/env python3
"""CollecTRI/ULM analysis of JMML vs Control from an AnnData counts layer.

Outputs sample-level TF activities, differential activities, DEG-TF intersection,
and optional TF->target network plots.

Rank DEG TFs by condition-level ULM score in JMML, requiring JMML > Control.

Stage 1 (expensive, only once):
  python run_collectri.py --input data/jmml.h5ad --outdir output/collectri_results

Stage 2 (fast, whenever your DEG candidates change):
  python run_collectri.py --reuse --reuse-dir output/collectri_results \
      --deg data/CD14_Mono_candidates_GTEx_v11.csv --outdir output/collectri_cd14

Important: --sample-col must identify INDEPENDENT BIOLOGICAL samples for
statistical testing. A purely technical sequencing batch is not sufficient.
"""

from __future__ import annotations

import argparse
import json
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from lazy_h5ad import open_counts  # noqa: E402
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

LOG = logging.getLogger("collectri")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", type=Path, help="AnnData .h5ad input for activity stage")
    p.add_argument("--deg", type=Path, help="Candidate/DEG CSV or TSV; optional in activity stage")
    p.add_argument("--reuse", action="store_true", help="Reuse saved activity results; input h5ad not needed")
    p.add_argument("--reuse-dir", type=Path, help="Directory containing previous CollecTRI results (read-only); required with --reuse")
    p.add_argument("--outdir", type=Path, required=True,
                   help="Output directory path, supplied by the workflow/config")
    p.add_argument("--counts-layer", default="counts", help="Raw counts: a name in adata.layers, or X for adata.X")
    p.add_argument("--sample-col", default="batch", help="Independent biological sample ID")
    p.add_argument("--group-col", default="sample_groups_wide")
    p.add_argument("--group-map", default=None,
                   help="JSON mapping of observed metadata labels to analysis groups")
    p.add_argument("--celltype-col", default="celltype_l2_merged_rename")
    p.add_argument("--groups", nargs=2, default=["jmml", "control"], metavar=("CASE", "CONTROL"))
    p.add_argument("--celltypes", nargs="+", default=["hsc", "gmp", "cd14 mono"])
    p.add_argument("--gene-col", default="gene", help="Gene-symbol column in DEG file")
    p.add_argument("--deg-celltype-col", default="celltype", help="Column used if present in DEG CSV")
    p.add_argument("--deg-celltype", default="cd14 mono", help="Assigned cell type when DEG file has no celltype column")
    p.add_argument("--min-cells", type=int, default=20, help="Minimum cells per sample x celltype profile")
    p.add_argument("--min-samples", type=int, default=2, help="Minimum samples per condition for testing")
    p.add_argument("--min-targets", type=int, default=5, help="Minimum measured targets per TF for ULM")
    p.add_argument("--batch-size", type=int, default=256)
    p.add_argument("--plot-celltype", default="cd14 mono")
    p.add_argument("--top", type=int, default=None, help="Deprecated: all DEG TFs with JMML > Control are saved; use --plot-top to limit plots")
    p.add_argument("--plot-top", type=int, default=5, help="Number of top-ranked eligible DEG TFs to visualize")
    p.add_argument("--network-targets", type=int, default=15)
    p.add_argument("--plot-tfs", nargs="+", default=None, help="Explicit TFs to plot (e.g. ZEB1 GATA2); overrides --plot-top")
    p.add_argument("--no-plots", action="store_true", help="Disable PDF network plots")
    p.add_argument("--backed", action="store_true",
                   help="Keep the counts matrix on disk and read it row-wise (low RAM, potentially slower)")
    return p.parse_args()


def norm(s: pd.Series) -> pd.Series:
    return s.astype("string").str.strip().str.lower()


def select_groups(adata, args: argparse.Namespace):
    for col in (args.sample_col, args.group_col, args.celltype_col):
        if col not in adata.obs:
            raise KeyError(f"Missing adata.obs column: {col!r}. Available: {list(adata.obs.columns)}")

    sample = adata.obs[args.sample_col].astype("string").str.strip()
    group = norm(adata.obs[args.group_col])
    if args.group_map:
        mapping = {str(k).strip().lower(): str(v).strip().lower()
                   for k, v in json.loads(args.group_map).items()}
        group = group.map(mapping).fillna(group)
    ct = norm(adata.obs[args.celltype_col])
    keep = (
        group.isin(args.groups)
        & ct.isin(args.celltypes)
        & sample.notna()
        & sample.ne("")
    ).fillna(False).to_numpy(dtype=bool)
    cell_idx = np.flatnonzero(keep)
    if len(cell_idx) == 0:
        raise ValueError("No matching cells. Check group/celltype labels and metadata columns.")

    meta = pd.DataFrame({
        args.sample_col: sample.iloc[cell_idx].to_numpy(),
        args.group_col: group.iloc[cell_idx].to_numpy(),
        args.celltype_col: ct.iloc[cell_idx].to_numpy(),
    }).reset_index(drop=True)

    # A single sample cannot be simultaneously a JMML and control replicate.
    collision = meta.groupby(args.sample_col, observed=True)[args.group_col].nunique()
    if (collision > 1).any():
        bad = collision[collision > 1].index.tolist()[:10]
        raise ValueError(
            f"Some {args.sample_col!r} identifiers occur in both conditions: {bad}. "
            "Provide a unique biological-sample column (--sample-col)."
        )

    group_cols = [args.sample_col, args.group_col, args.celltype_col]
    meta["gid"] = meta.groupby(group_cols, sort=True, observed=True).ngroup()
    sizes = meta.groupby("gid").size()
    good = meta["gid"].isin(sizes.index[sizes >= args.min_cells]).to_numpy()
    removed = int((~good).sum())
    cell_idx, meta = cell_idx[good], meta.iloc[np.flatnonzero(good)].reset_index(drop=True)
    if meta.empty:
        raise ValueError(f"All sample x celltype groups have fewer than {args.min_cells} cells")

    meta["gid"] = meta.groupby(group_cols, sort=True, observed=True).ngroup()
    group_codes = meta["gid"].to_numpy(dtype=np.int32)
    group_info = (meta.drop_duplicates("gid")[["gid"] + group_cols]
                      .sort_values("gid").reset_index(drop=True))
    group_info["group_id"] = "G" + group_info["gid"].astype(str)
    group_info["n_cells"] = np.bincount(group_codes, minlength=len(group_info))
    LOG.info("Selected %s cells in %s profiles; filtered %s cells in small profiles",
             len(cell_idx), len(group_info), removed)
    return cell_idx, group_codes, group_info


def load_network(adata, args, dc):
    # Note: .var_names must contain HGNC gene symbols rather than Ensembl IDs.
    if not adata.var_names.is_unique:
        raise ValueError("adata.var_names must contain unique gene names")
    network = dc.op.collectri(organism="human", remove_complexes=True)
    network = network.dropna(subset=["source", "target"]).copy()
    network["source"] = network["source"].astype(str).str.strip()
    network["target"] = network["target"].astype(str).str.strip()

    gene_names = pd.Index(adata.var_names.astype(str))
    target_mask = gene_names.isin(network["target"])
    target_idx = np.flatnonzero(target_mask)
    if len(target_idx) == 0:
        raise ValueError("No CollecTRI target genes matched adata.var_names. Are these HGNC symbols or Ensembl IDs?")
    covered = network[network["target"].isin(gene_names)].copy()

    total = network.groupby("source")["target"].nunique().rename("total_targets")
    observed = covered.groupby("source")["target"].nunique().rename("available_targets")
    coverage = total.to_frame().join(observed).fillna({"available_targets": 0})
    coverage["available_targets"] = coverage["available_targets"].astype(int)
    coverage["coverage"] = coverage["available_targets"] / coverage["total_targets"]
    coverage["passes_min_targets"] = coverage["available_targets"] >= args.min_targets
    coverage.index.name = "gene"
    coverage = coverage.reset_index()
    LOG.info("CollecTRI: %s edges, %s TFs, %s/%s targets in adata",
             len(network), network["source"].nunique(), len(target_idx), network["target"].nunique())
    return network, covered, coverage, target_idx


def make_profiles(adata, cell_idx, group_codes, group_info, target_idx, args):
    if args.counts_layer != "X" and args.counts_layer not in adata.layers:
        raise KeyError(f"Missing adata.layers[{args.counts_layer!r}]; available: {list(adata.layers.keys())}")
    counts = adata.X if args.counts_layer == "X" else adata.layers[args.counts_layer]
    n_groups = len(group_info)
    n_targets = len(target_idx)
    sums = np.zeros((n_groups, n_targets), dtype=np.float64)
    n_cells = group_info["n_cells"].to_numpy()
    for start in range(0, len(cell_idx), args.batch_size):
        end = min(start + args.batch_size, len(cell_idx))
        idx = cell_idx[start:end]
        x = counts[idx, :]
        x = x.tocsr().astype(np.float32) if sp.issparse(x) else sp.csr_matrix(x, dtype=np.float32)
        if x.nnz and (np.any(x.data < 0) or not np.all(np.isfinite(x.data))):
            raise ValueError("Counts layer contains non-finite or negative entries")
        total_counts = np.asarray(x.sum(axis=1)).ravel()
        if np.any(total_counts == 0):
            LOG.warning("Found %d cells with zero library size; they contribute zero normalized counts", int(np.sum(total_counts == 0)))
        scale = (1e4 / np.maximum(total_counts, 1)).astype(np.float32)
        # Denominator includes ALL genes, even though only network targets are stored.
        xt = x[:, target_idx].multiply(scale[:, None]).tocsr()
        xt.data = np.log1p(xt.data)
        local_codes = group_codes[start:end]
        membership = sp.csr_matrix(
            (np.ones(len(local_codes), dtype=np.float32),
             (local_codes, np.arange(len(local_codes)))),
            shape=(n_groups, len(local_codes)),
        )
        sums += (membership @ xt).toarray()
        if start == 0 or end == len(cell_idx) or end % 12800 == 0:
            LOG.info("Processed %s / %s cells", end, len(cell_idx))

    profiles = pd.DataFrame(
        sums / n_cells[:, None],
        index=group_info["group_id"],
        columns=pd.Index(adata.var_names[target_idx].astype(str)),
    )
    profiles.index.name = "group_id"
    return profiles


def compute_activity(profiles, covered, args, dc):
    net_usable = covered.copy()
    counts = net_usable.groupby("source")["target"].nunique()
    net_usable = net_usable[net_usable["source"].isin(counts.index[counts >= args.min_targets])]
    if net_usable.empty:
        raise ValueError("No TF has enough measured targets for --min-targets")
    tf_scores, tf_pvals = dc.mt.ulm(data=profiles, net=net_usable, tmin=args.min_targets)
    if not tf_scores.index.equals(profiles.index):
        raise RuntimeError("Unexpected row index/order from decoupler.mt.ulm")
    if tf_scores.empty:
        raise RuntimeError("ULM returned an empty activity matrix")
    LOG.info("ULM activities: %s profiles x %s TFs", *tf_scores.shape)
    return tf_scores, tf_pvals, net_usable


def summarize_activity(tf_scores, group_info, args):
    ids = group_info.set_index("group_id")
    rows = []
    for ct in args.celltypes:
        ct_info = ids[ids[args.celltype_col] == ct]
        case_ids = ct_info.index[ct_info[args.group_col] == args.groups[0]]
        ctrl_ids = ct_info.index[ct_info[args.group_col] == args.groups[1]]
        for tf in tf_scores.columns:
            a = tf_scores.loc[case_ids, tf].dropna().to_numpy(dtype=float)
            b = tf_scores.loc[ctrl_ids, tf].dropna().to_numpy(dtype=float)
            pval = np.nan
            if len(a) >= args.min_samples and len(b) >= args.min_samples:
                pval = float(mannwhitneyu(a, b, alternative="two-sided").pvalue)
            mean_a = float(a.mean()) if len(a) else np.nan
            mean_b = float(b.mean()) if len(b) else np.nan
            rows.append({
                "celltype": ct,
                "gene": tf,
                "n_jmml": len(a),
                "n_control": len(b),
                "mean_jmml": mean_a,
                "mean_control": mean_b,
                "delta_activity": mean_a - mean_b,
                "pvalue": pval,
            })
    summary = pd.DataFrame(rows)
    summary["padj"] = np.nan
    for ct, pos in summary.groupby("celltype", sort=False).groups.items():
        idx = pd.Index(pos)
        ok = idx[summary.loc[idx, "pvalue"].notna()]
        if len(ok):
            summary.loc[ok, "padj"] = multipletests(summary.loc[ok, "pvalue"], method="fdr_bh")[1]
    return summary.sort_values(["celltype", "padj", "gene"], na_position="last").reset_index(drop=True)



def compute_condition_ulm(profiles, group_info, net, args, dc):
    """ULM of each condition's mean expression profile, NOT a group-comparison p-value.

    One row per cell type and condition; each biological sample has equal weight.
    The ULM p-value is conditional on the regulator network and gene expression
    across targets. It DOES NOT account for patient-to-patient variability.
    """
    info = group_info.set_index("group_id")
    rows = []
    labels = []
    for ct in args.celltypes:
        for condition in args.groups:
            ids = info.index[(info[args.celltype_col] == ct) & (info[args.group_col] == condition)]
            if len(ids) == 0:
                continue
            labels.append((ct, condition, len(ids)))
            rows.append(profiles.loc[ids].mean(axis=0))
    if not rows:
        raise ValueError("No profiles available for condition ULM")
    mat = pd.DataFrame(rows, index=[f"condition_{i}" for i in range(len(rows))])
    scores, pvals = dc.mt.ulm(data=mat, net=net, tmin=args.min_targets)
    long = pd.concat([
        scores.stack().rename("ulm_group_score"),
        pvals.stack().rename("ulm_group_pvalue"),
    ], axis=1).reset_index()
    long.columns = ["profile_id", "gene", "ulm_group_score", "ulm_group_pvalue"]
    map_df = pd.DataFrame({
        "profile_id": mat.index,
        "celltype": [x[0] for x in labels],
        "condition": [x[1] for x in labels],
        "n_samples": [x[2] for x in labels],
    })
    long = long.merge(map_df, on="profile_id", validate="many_to_one")
    return long, mat


def attach_condition_ulm(summary, condition_ulm, args):
    # Safe when reuse-dir already contains a summary enriched by a previous run.
    stale = [c for c in summary.columns if c.startswith("ulm_group_") or c.startswith("ulm_score_") or c.startswith("ulm_pvalue_")]
    result = summary.drop(columns=stale).copy()
    for condition in args.groups:
        short = condition_ulm[condition_ulm["condition"] == condition]
        short = short[["celltype", "gene", "ulm_group_score", "ulm_group_pvalue"]].rename(columns={
            "ulm_group_score": f"ulm_score_{condition}",
            "ulm_group_pvalue": f"ulm_pvalue_{condition}",
        })
        result = result.merge(short, on=["celltype", "gene"], how="left", validate="one_to_one")
    return result


def match_candidates(summary, network, args):
    """Rank all CollecTRI TFs, annotate supplied DEGs, and select the top N.

    Ranks are assigned across all TFs with available condition-level ULM
    scores *before* intersecting with DEG genes, so they are comparable even
    when the DEG list changes. Larger ULM score receives smaller rank (1=best).
    """
    deg = pd.read_csv(args.deg, sep=None, engine="python")
    if args.gene_col not in deg.columns:
        raise KeyError(f"Missing DEG gene column {args.gene_col!r}; found: {list(deg.columns)}")
    deg = deg.rename(columns={args.gene_col: "gene"})
    if args.deg_celltype_col in deg.columns:
        deg["celltype"] = norm(deg[args.deg_celltype_col])
    else:
        deg["celltype"] = args.deg_celltype.strip().lower()
    deg["gene"] = deg["gene"].astype("string").str.strip().str.upper()
    deg = deg.dropna(subset=["gene", "celltype"])

    score_j = f"ulm_score_{args.groups[0]}"
    score_c = f"ulm_score_{args.groups[1]}"
    if score_j not in summary or score_c not in summary:
        raise KeyError(f"Condition ULM scores missing: {score_j}, {score_c}")

    ranked = summary.copy()
    ranked["gene"] = ranked["gene"].astype(str).str.upper()
    for col in (score_j, score_c):
        ranked[col] = pd.to_numeric(ranked[col], errors="coerce")
    ranked["rank_jmml"] = ranked.groupby("celltype")[score_j].rank(method="min", ascending=False)
    ranked["rank_control"] = ranked.groupby("celltype")[score_c].rank(method="min", ascending=False)
    ranked["delta_ulm_score"] = ranked[score_j] - ranked[score_c]
    ranked["jmml_higher"] = ranked[score_j].gt(ranked[score_c]) & ranked[[score_j, score_c]].notna().all(axis=1)

    # Keep only TFs known to the original CollecTRI network (including those
    # with insufficient observed targets to obtain scores).
    tf_sources = set(network["source"].astype(str).str.upper())
    deg_tfs = deg[deg["gene"].isin(tf_sources)].copy()
    collisions = set(deg_tfs.columns) & (set(ranked.columns) - {"gene", "celltype"})
    if collisions:
        # Protect original DEG columns by labeling them explicitly.
        deg_tfs = deg_tfs.rename(columns={c: f"deg_{c}" for c in collisions})
    results = deg_tfs.merge(ranked, on=["celltype", "gene"], how="left", validate="many_to_one")
    results["activity_available"] = results[[score_j, score_c]].notna().all(axis=1)
    results["jmml_higher"] = results["jmml_higher"].fillna(False).astype(bool)
    results["selected_top"] = False
    # Only candidates whose JMML score is strictly higher than Control enter
    # the shortlist; do not filter on differential activity p-value.
    eligible = results.loc[results["jmml_higher"]].copy()
    eligible = eligible.sort_values(["celltype", score_j, "gene"], ascending=[True, False, True])
    eligible = eligible.drop_duplicates(["celltype", "gene"])
    # Save ALL eligible TFs, not just the first 20. --plot-top only controls visualization.
    shortlisted = eligible.copy()
    results.loc[shortlisted.index, "selected_top"] = True
    results = results.sort_values(
        ["celltype", "selected_top", score_j, "gene"],
        ascending=[True, False, False, True], na_position="last",
    ).reset_index(drop=True)
    shortlisted = results.loc[results["selected_top"]].copy()
    LOG.info("DEG rows: %d; CollecTRI TF rows: %d; JMML > Control: %d; shortlisted: %d",
             len(deg), len(results), int(results["jmml_higher"].sum()), len(shortlisted))
    return results, shortlisted, ranked


def safe_filename(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_")


def plot_networks(profiles, group_info, summary, candidates, net, args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as pe
    import decoupler as dc

    ct = args.plot_celltype.strip().lower()
    info = group_info.set_index("group_id")
    info = info.loc[info[args.celltype_col] == ct]
    j_ids = info.index[info[args.group_col] == args.groups[0]]
    c_ids = info.index[info[args.group_col] == args.groups[1]]
    if len(j_ids) == 0 or len(c_ids) == 0:
        LOG.warning("Skipping plots: need both JMML and Control profiles for %s", ct)
        return

    jmml_expr = profiles.loc[j_ids].mean(axis=0)
    ctrl_expr = profiles.loc[c_ids].mean(axis=0)
    expr_change = pd.DataFrame([jmml_expr - ctrl_expr], index=["JMML_vs_Control"])

    # Compute scores from condition-average profiles for condition-specific networks.
    condition_expr = pd.DataFrame([jmml_expr, ctrl_expr], index=["JMML", "Control"])
    condition_scores, _ = dc.mt.ulm(data=condition_expr, net=net, tmin=args.min_targets)
    subsummary = summary[(summary["celltype"] == ct) & summary["delta_activity"].notna()]
    score_change = subsummary.set_index("gene")["delta_activity"].to_frame().T
    score_change.index = ["JMML_vs_Control"]

    if args.plot_tfs:
        top = [x.strip().upper() for x in args.plot_tfs]
    else:
        sub = candidates[(candidates["celltype"] == ct) & candidates["selected_top"]]
        top = (sub.sort_values(f"ulm_score_{args.groups[0]}", ascending=False)
               ["gene"].drop_duplicates().head(args.plot_top).tolist())
    permitted = set(net["source"]) & set(condition_scores.columns)
    top = list(dict.fromkeys([tf for tf in top if tf in permitted]))
    if not top:
        LOG.warning("No candidate TF with estimated activity to plot for %s", ct)
        return

    out = args.outdir / "networks" / safe_filename(ct)
    out.mkdir(parents=True, exist_ok=True)
    plot_specs = [(tf, [tf]) for tf in top]
    if len(top) > 1:
        plot_specs.append(("combined", top))
    for label, selected in plot_specs:
        for condition, expr_data, score_data in [
            ("JMML", condition_expr.loc[["JMML"]], condition_scores.loc[["JMML"]]),
            ("Control", condition_expr.loc[["Control"]], condition_scores.loc[["Control"]]),
            ("JMML_vs_Control", expr_change, score_change),
        ]:
            try:
                is_delta = condition == "JMML_vs_Control"
                fig = dc.pl.network(
                    net=net, data=expr_data, score=score_data,
                    sources=selected, targets=args.network_targets,
                    size_node=10, size_label=4.5, figsize=(12, 9),
                    s_cmap="RdBu_r", t_cmap="RdBu_r" if is_delta else "viridis",
                    vcenter=is_delta,
                    s_label="Delta TF activity" if is_delta else f"{condition} TF activity (ULM)",
                    t_label="Delta target expression" if is_delta else f"{condition} mean log1p(CP10K)",
                    return_fig=True,
                )
                # Make gene labels readable over colored network nodes and edges.
                # A white halo increases contrast without obscuring the node color.
                gene_labels = set(net["source"]) | set(net["target"])
                for ax in fig.axes:
                    for artist in ax.texts:
                        if artist.get_text().strip() in gene_labels:
                            artist.set_fontweight("bold")
                            artist.set_path_effects([
                                pe.Stroke(linewidth=3.0, foreground="white"),
                                pe.Normal(),
                            ])
                            artist.set_zorder(20)
                path = out / f"{safe_filename(label)}_{condition}.pdf"
                fig.savefig(path, bbox_inches="tight")
                plt.close(fig)
                LOG.info("Saved network: %s", path)
            except (ValueError, ImportError, RuntimeError, AttributeError, TypeError) as exc:
                LOG.warning("Network plot failed for %s (%s): %s", label, condition, exc)


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args.groups = [v.strip().lower() for v in args.groups]
    args.celltypes = [v.strip().lower() for v in args.celltypes]
    if args.min_cells < 1 or args.min_samples < 1 or args.min_targets < 1 or args.batch_size < 1:
        raise ValueError("min-cells, min-samples, min-targets and batch-size must be positive")
    if args.groups[0] == args.groups[1]:
        raise ValueError("--groups must be two different conditions")
    if args.deg is not None and not args.deg.is_file():
        raise FileNotFoundError(f"DEG file does not exist: {args.deg}")
    if args.reuse:
        if args.reuse_dir is None:
            raise ValueError("--reuse requires --reuse-dir to specify where existing results are read from")
        if args.reuse_dir.resolve() == args.outdir.resolve():
            raise ValueError("--reuse-dir and --outdir must be different directories")
        required = [
            "expression_profiles.csv", "sample_celltype_info.csv",
            "tf_activity_summary.csv", "collectri_network_full.csv",
            "collectri_network.csv"
        ]
        missing = [name for name in required if not (args.reuse_dir / name).is_file()]
        if missing:
            raise FileNotFoundError(f"Missing files in --reuse-dir {args.reuse_dir}: {', '.join(missing)}")
    elif args.reuse_dir is not None:
        raise ValueError("--reuse-dir can only be used together with --reuse")

    args.outdir.mkdir(parents=True, exist_ok=True)

    if args.reuse:
        if args.deg is None:
            raise ValueError("--reuse requires --deg to select candidate TFs")
        LOG.info("Reusing previously calculated CollecTRI profiles from %s", args.reuse_dir)
        profiles = pd.read_csv(args.reuse_dir / "expression_profiles.csv", index_col=0)
        group_info = pd.read_csv(args.reuse_dir / "sample_celltype_info.csv")
        summary = pd.read_csv(args.reuse_dir / "tf_activity_summary.csv")
        net_full = pd.read_csv(args.reuse_dir / "collectri_network_full.csv")
        net_usable = pd.read_csv(args.reuse_dir / "collectri_network.csv")
    else:
        if args.input is None or not args.input.is_file():
            raise FileNotFoundError("Supply an existing --input .h5ad, or run with --reuse")
        import anndata as ad
        import decoupler as dc
        LOG.info("Reading %s", args.input)
        # anndata's backed mode would still load every layer into RAM; open_counts does not.
        handle, adata = (open_counts(args.input, args.counts_layer) if args.backed
                         else (None, ad.read_h5ad(args.input)))
        try:
            cell_idx, group_codes, group_info = select_groups(adata, args)
            net_full, net_covered, coverage, target_idx = load_network(adata, args, dc)
            profiles = make_profiles(adata, cell_idx, group_codes, group_info, target_idx, args)
            tf_scores, tf_pvals, net_usable = compute_activity(profiles, net_covered, args, dc)
            summary = summarize_activity(tf_scores, group_info, args)

            group_info.to_csv(args.outdir / "sample_celltype_info.csv", index=False)
            coverage.to_csv(args.outdir / "tf_network_coverage.csv", index=False)
            net_full.to_csv(args.outdir / "collectri_network_full.csv", index=False)
            net_usable.to_csv(args.outdir / "collectri_network.csv", index=False)
            profiles.to_csv(args.outdir / "expression_profiles.csv")
            # Wide format: one row per sample x celltype.
            activity = group_info.set_index("group_id").join(tf_scores).reset_index()
            activity.to_csv(args.outdir / "tf_activity.csv", index=False)
            group_info.set_index("group_id").join(tf_pvals).reset_index().to_csv(
                args.outdir / "tf_ulm_pvalues_by_sample.csv", index=False
            )
            summary.to_csv(args.outdir / "tf_activity_summary.csv", index=False)
        finally:
            if handle is not None:
                handle.close()

    # Calculate group-profile ULM p-values for BOTH fresh and --reuse runs.
    # These are within-condition enrichment tests, not p-values of JMML vs Control.
    import decoupler as dc
    condition_ulm, condition_profiles = compute_condition_ulm(profiles, group_info, net_usable, args, dc)
    summary = attach_condition_ulm(summary, condition_ulm, args)
    condition_ulm.to_csv(args.outdir / "tf_ulm_by_condition.csv", index=False)
    condition_profiles.to_csv(args.outdir / "expression_by_condition.csv")
    summary.to_csv(args.outdir / "tf_activity_summary.csv", index=False)

    if args.deg is not None:
        candidates, shortlisted, ranked = match_candidates(summary, net_full, args)
        ranked.to_csv(args.outdir / "all_tf_ranked.csv", index=False)
        candidates.to_csv(args.outdir / "candidate_tf_activity.csv", index=False)
        shortlisted.to_csv(args.outdir / "jmml_higher_deg_tf.csv", index=False)
        # Backward-compatible filename; now also contains ALL eligible TFs.
        shortlisted.to_csv(args.outdir / "top_deg_tf_jmml.csv", index=False)
        if not args.no_plots and args.plot_top > 0:
            plot_networks(profiles, group_info, summary, candidates, net_usable, args)
    LOG.info("Done. Outputs: %s", args.outdir.resolve())


if __name__ == "__main__":
    main()
