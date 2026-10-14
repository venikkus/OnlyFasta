#!/usr/bin/env python3
r"""Final candidate table from the outputs of one OnlyFasta workflow run.

Candidates per population are the union of two branches:
  DEG branch - surface Candidate = Yes (upregulated, membrane, low normal-tissue GTEx expression);
  TF branch  - DEG transcription factors with higher CollecTRI activity in case than control,
               plus CellOracle rank candidates.
Modules that were switched off leave their columns empty.

Outputs in --outdir:
  final_candidates.csv              one row per population x candidate gene
  final_candidates.xlsx             sheet 'all' + one sheet per population
  opentargets_report_examples.md    Open Targets report for the top candidates (if available)
  README.md                         column descriptions

CLI: python build_summary.py --results-dir results/<run_name> --populations HSC GMP \
         --case JMML --control Control --outdir results/<run_name>/summary
"""
import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DEFAULT_OT_REPORT = HERE.parent / "opentarget" / "ot_report_module.py"
LIANA_FILE = "JMML_vs_Control_DEG_receptors.csv"


def read_optional(path):
    path = Path(path)
    if not path.is_file():
        return None
    table = pd.read_csv(path)
    if "gene" in table:
        # CollecTRI upper-cases symbols; match all tables on an upper-case key.
        table["key"] = table["gene"].astype(str).str.strip().str.upper()
    return table


def by_key(table, column):
    if table is None or column not in table:
        return {}
    return table.drop_duplicates("key").set_index("key")[column].to_dict()


def yes_no(value):
    if pd.isna(value):
        return np.nan
    return "Yes" if str(value).strip().lower() in {"true", "yes", "1"} else "No"


def liana_top(path, n=5):
    """Top-n interactions per DEG receptor by gain in case vs control (delta_lr_probs)."""
    table = read_optional(path)
    if table is None or table.empty:
        return {}
    table = table.dropna(subset=["DEG_as_receptor", "delta_lr_probs"])
    table = table.sort_values("delta_lr_probs", ascending=False, kind="mergesort")
    table["pair"] = table["source"].astype(str) + ": " + table["ligand"].astype(str)
    top = table.groupby(table["DEG_as_receptor"].astype(str).str.upper(), sort=False)["pair"]
    return top.apply(lambda pairs: ", ".join(pairs.drop_duplicates().head(n))).to_dict()


def population_table(root, population, case_label):
    root = Path(root)
    de = read_optional(root / "de" / population / "genes.csv")
    if de is None:
        raise FileNotFoundError(f"{population}: no DE table in {root / 'de' / population}")
    surface = read_optional(root / "surface" / population / "DEGs_surface_GTEx.csv")
    collectri = read_optional(root / "collectri" / population / "candidate_tf_activity.csv")
    celloracle = read_optional(root / "celloracle" / f"{population}_rank_candidates.csv")
    depmap = read_optional(root / "depmap" / population / "features.csv")
    opentargets = read_optional(root / "opentargets" / population / "opentargets_features.csv")
    liana = liana_top(root / "liana" / population / LIANA_FILE)

    deg_branch = [] if surface is None else surface.loc[surface["Candidate"].eq("Yes"), "key"].tolist()
    tf_branch = []
    if collectri is not None:
        selected = "selected_top" if "selected_top" in collectri else "jmml_higher"
        tf_branch += collectri.loc[collectri[selected].astype(str).str.lower().eq("true"), "key"].tolist()
    if celloracle is not None:
        tf_branch += celloracle["key"].tolist()
    keys = list(dict.fromkeys(deg_branch + tf_branch))

    # Display the symbol as written in the DE table when available.
    symbol = {**by_key(celloracle, "gene"), **by_key(collectri, "gene"), **by_key(de, "gene")}
    status = by_key(opentargets, "ot_request_status")
    found = by_key(opentargets, "ot_target_found")
    atlases = by_key(surface, "Atlases")
    tumor_specific = by_key(surface, "Tumor_specific")
    surface_flag = by_key(surface, "Surface")
    column = {name: by_key(table, col) for name, table, col in [
        ("name", opentargets, "ot_approved_name"), ("pathways", opentargets, "ot_pathway_top_level"),
        ("lfc", de, "log2FoldChange"), ("padj", de, "padj"), ("tpm", surface, "GTEx_max_TPM"),
        ("depmap", depmap, "depmap_selectivity"), ("co_rank", celloracle, "rank_case"),
        ("co_delta", celloracle, "delta_score"), ("tf_delta", collectri, "delta_activity"),
        ("tf_p", collectri, "pvalue")]}

    def values(name):
        return [column[name].get(k, np.nan) for k in keys]

    def surface_text(key):
        if key not in surface_flag:
            return np.nan
        return f"Yes: {atlases[key]}" if surface_flag[key] == "Membrane" else "No"

    table = pd.DataFrame({
        "Gene": [symbol.get(k, k) for k in keys],
        "Gene_name": values("name"),
        "Reactome_pathways": values("pathways"),
        "DE_log2FC": values("lfc"),
        "DE_padj": values("padj"),
        "Surface": [surface_text(k) for k in keys],
        "Tumor_specific": [{"Specific": "Yes", "Not_specific": "No", "No_GTEx_data": "No GTEx data"}
                           .get(tumor_specific.get(k), np.nan) for k in keys],
        "GTEx_max_TPM": values("tpm"),
        "LIANA_top5": [liana.get(k, np.nan) for k in keys],
        "DepMap_tumor_vs_other": values("depmap"),
        "OpenTargets": [(yes_no(found.get(k)) if status.get(k) in {"ok", "not_found"} else np.nan)
                        if k in status else np.nan for k in keys],
        "TF": ["Yes" if k in tf_branch else "No" for k in keys],
        f"CellOracle_rank_{case_label}": values("co_rank"),
        "CellOracle_delta_score": values("co_delta"),
        "TF_activity_delta": values("tf_delta"),
        "TF_activity_pvalue": values("tf_p"),
    })
    table = table.sort_values(["DE_padj", "Gene"], na_position="last", kind="mergesort")
    return table.reset_index(drop=True), opentargets


def load_report_module(path):
    spec = importlib.util.spec_from_file_location("ot_report_module", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def column_notes(case_label, control_label):
    return f"""# Final candidate table

`final_candidates.csv` / `final_candidates.xlsx`: one row per population and candidate gene.
Empty cells mean the module was switched off or has no data for the gene, not a negative result.

| Column | Meaning |
|---|---|
| Gene | Candidate gene symbol |
| Population | Population from the config (combined table only) |
| Gene_name | Approved gene name (Open Targets) |
| Reactome_pathways | Reactome top-level pathways (Open Targets) |
| DE_log2FC, DE_padj | Pseudobulk DE, {case_label} vs {control_label} |
| Surface | Yes + atlases supporting membrane localisation (HPA, CSPA, TCSA), or No |
| Tumor_specific | Yes if max normal-tissue GTEx TPM is below the threshold (surface.gtex_max_tpm) |
| GTEx_max_TPM | Max median TPM across included normal GTEx tissues |
| LIANA_top5 | Up to 5 interactions (source cell type: ligand) with this gene as receptor, ranked by gain in {case_label} vs {control_label} |
| DepMap_tumor_vs_other | CRISPR dependency selectivity: mean gene effect in other cancer models minus in the selected disease models; > 0 = stronger dependency in the selected models. DepMap has no normal controls |
| OpenTargets | Yes/No: target found in Open Targets; details for top genes in `opentargets_report_examples.md` |
| TF | Yes if the gene comes from the regulatory branch: a DEG transcription factor with higher CollecTRI activity in {case_label}, or a CellOracle GRN hub (hubs are not necessarily TFs). No if only from the DEG/surface branch |
| CellOracle_rank_{case_label}, CellOracle_delta_score | GRN centrality percentile rank in {case_label} and {case_label} - {control_label} score difference |
| TF_activity_delta, TF_activity_pvalue | CollecTRI/decoupler TF activity difference {case_label} - {control_label} between samples and its p-value |
"""


def build_summary(results_dir, populations, case_label, control_label, outdir,
                  report_genes=2, ot_report_path=DEFAULT_OT_REPORT):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    tables, features = {}, []
    for population in populations:
        table, opentargets = population_table(results_dir, population, case_label)
        tables[population] = table
        if opentargets is not None:
            features.append(opentargets.assign(_population=population))
    combined = pd.concat([t.assign(Population=p) for p, t in tables.items()], ignore_index=True)
    combined = combined[["Gene", "Population", *[c for c in combined if c not in {"Gene", "Population"}]]]
    combined.to_csv(outdir / "final_candidates.csv", index=False)
    with pd.ExcelWriter(outdir / "final_candidates.xlsx") as writer:
        combined.to_excel(writer, sheet_name="all", index=False)
        for population, table in tables.items():
            table.to_excel(writer, sheet_name=population[:31], index=False)
    (outdir / "README.md").write_text(column_notes(case_label, control_label), encoding="utf-8")

    # Open Targets report for the best-ranked candidates found in Open Targets.
    if report_genes > 0 and features:
        ranked = combined.loc[combined["OpenTargets"].eq("Yes")].sort_values("DE_padj", na_position="last")
        picked = ranked.drop_duplicates("Gene").head(report_genes)
        if len(picked):
            pool = pd.concat(features, ignore_index=True)
            rows = [pool.loc[(pool["key"] == g.upper()) & (pool["_population"] == p)].head(1)
                    for g, p in zip(picked["Gene"], picked["Population"])]
            example = pd.concat(rows, ignore_index=True).drop(columns=["key", "_population"])
            example_csv = outdir / "opentargets_report_examples.csv"
            example.to_csv(example_csv, index=False)
            load_report_module(ot_report_path).make_report(
                example_csv, outdir / "opentargets_report_examples.md", max_genes=report_genes)
    return combined


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", type=Path, required=True, help="results/<run_name> of a workflow run")
    ap.add_argument("--populations", nargs="+", required=True)
    ap.add_argument("--case", required=True, help="case label, e.g. JMML")
    ap.add_argument("--control", required=True, help="control label, e.g. Control")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--report-genes", type=int, default=2, help="genes in the Open Targets example report; 0 = none")
    args = ap.parse_args()
    table = build_summary(args.results_dir, args.populations, args.case, args.control, args.outdir,
                          report_genes=args.report_genes)
    print(f"{len(table)} candidate rows -> {args.outdir / 'final_candidates.csv'}")


if __name__ == "__main__":
    main()
