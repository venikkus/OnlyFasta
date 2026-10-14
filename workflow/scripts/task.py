"""Snakemake adapters: original tools stay unchanged; inputs remain read-only."""
import os
import sys
from pathlib import Path

# Set limits BEFORE importing numerical libraries, including child processes.
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[name] = str(snakemake.threads)
os.environ["MPLBACKEND"] = "Agg"

import contextlib
import hashlib
import importlib.util
import importlib.metadata
import json
import shutil
import subprocess
import tempfile
import traceback
from datetime import datetime, timezone
import numpy as np
import pandas as pd

C = snakemake.params.cfg
STAGE = snakemake.params.stage
POP = getattr(snakemake.wildcards, "population", None)
DS = C["datasets"].get(POP, {})
I = snakemake.input
QC = {}
WARNINGS = []
DECOLS = ["gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]


def table(path, required=(), unique=None, allow_empty=True):
    d = pd.read_csv(path)
    missing = set(required) - set(d.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    if not allow_empty and d.empty:
        raise ValueError(f"{path}: unexpectedly empty table")
    if unique:
        if d[unique].isna().any().any() or d.duplicated(unique).any():
            raise ValueError(f"{path}: null or duplicate key {unique}")
    return d


def tool_module(path):
    spec = importlib.util.spec_from_file_location("onlyfasta_tool", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run(*args):
    args = [str(a) for a in args]
    print("COMMAND:", json.dumps(args, ensure_ascii=False), flush=True)
    subprocess.run(args, check=True, stdout=sys.stdout, stderr=sys.stderr)


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, default=str) + "\n")


def fingerprint(path):
    p = Path(path)
    result = {"path": str(p.resolve()), "size": p.stat().st_size, "mtime_ns": p.stat().st_mtime_ns}
    if p.is_file() and p.stat().st_size <= 64 * 1024**2:
        result["sha256"] = hashlib.file_digest(p.open("rb"), "sha256").hexdigest()
    elif p.is_dir():
        result["files"] = [fingerprint(x) for x in sorted(p.iterdir()) if x.is_file() and not x.name.startswith(".")]
    else:
        result["sha256"] = None
        result["hash_note"] = "Large input: path/size/mtime recorded; no expensive full read."
    return result


def norm_deg(w):
    if C["pseudobulk_de"]["mode"] == "compute":
        source = Path(I.source) / f"{POP}_DEGs.csv"
        gene_column = "gene"
    else:
        source = I.source
        gene_column = DS["gene_column"]
    d = table(source, [gene_column, *DECOLS[1:]], allow_empty=False)
    d = d.rename(columns={gene_column: "gene"})[DECOLS].copy()
    if d.gene.isna().any():
        raise ValueError("Missing gene symbols in DEG input")
    d["gene"] = d.gene.astype(str).str.strip()
    if d.gene.eq("").any() or d.gene.duplicated().any():
        raise ValueError("Empty or duplicate gene symbols in DEG input; resolve explicitly")
    for col in DECOLS[1:]:
        d[col] = pd.to_numeric(d[col], errors="raise")
    for col in ["pvalue", "padj"]:
        if not d[col].dropna().between(0, 1).all():
            raise ValueError(f"Invalid {col} outside [0,1]")
    if not np.isfinite(d.log2FoldChange).all():
        raise ValueError("Non-finite log2FoldChange")
    d.to_csv(w / "genes.csv", index=False)
    QC["n_genes"] = len(d)
    QC["n_padj_lt_005"] = int(d.padj.lt(.05).sum())


def pseudobulk_de(w):
    cfg = C["pseudobulk_de"]
    module = tool_module(I.tool)
    import anndata as ad
    # The integrated JMML atlas is ~33 GB. Backed mode keeps the full matrix
    # on disk; pseudobulk_de only indexes selected cells from the raw counts layer.
    adata = ad.read_h5ad(I.adata, backed="r")
    options = {key: cfg[key] for key in (
        "sample_col", "celltype_col", "group_col", "block_col", "counts_layer", "group_map",
        "cell_types", "case_label", "control_label", "min_cells", "min_total_counts",
        "min_replicates", "min_gene_count", "fdr", "lfc_cutoff")}
    options["n_cpus"] = snakemake.threads
    try:
        results = module.run_pseudobulk_de(adata, w, **options)
    finally:
        adata.file.close()
    for population, result in results.items():
        if population not in C["populations"]:
            raise ValueError(f"DE configured output {population!r} is not in populations")
        table(w / f"{population}_DEGs.csv", ["gene", *DECOLS[1:]], allow_empty=True)
        QC[population] = {"n_significant_genes": len(result)}
    missing = set(C["populations"]) - set(results)
    if missing:
        raise ValueError(f"No pseudobulk DE result for configured populations: {sorted(missing)}")


def depmap_features(w):
    cfg = C["depmap"]
    if cfg["mode"] == "reuse":
        shutil.copy2(I.source[0], w / "features.csv")
    else:
        cmd = [sys.executable, I.tool, "--data-dir", cfg["raw_dir"], "--field", cfg["field"],
               "--out", w / "features.csv", "--min-models", cfg["min_models"]]
        if cfg["diseases"]:
            cmd += ["--diseases", *cfg["diseases"]]
        if cfg["subtypes"]:
            cmd += ["--subtypes", *cfg["subtypes"]]
        run(*cmd)
    d = table(w / "features.csv", ["gene", "depmap_strength", "depmap_fraction_dependent", "depmap_selectivity", "common_essential"], "gene", False)
    QC["n_features"] = len(d)
    QC["release"] = cfg["release"]
    if cfg["release"].startswith("unknown"):
        WARNINGS.append("DepMap release not documented; input checksum identifies this snapshot.")


def depmap(w):
    cfg = C["depmap"]
    run(sys.executable, I.tool, "--deg", Path(I.deg) / "genes.csv", "--depmap", Path(I.features) / "features.csv",
        "--gene-col", "gene", "--out", w / "features.csv", "--strength", cfg["strength"],
        "--prevalence", cfg["prevalence"], "--selectivity", cfg["selectivity"])
    d = table(w / "features.csv", ["gene", "depmap_strength"], "gene", False)
    QC["n_genes"] = len(d)
    QC["n_matched"] = int(d.depmap_strength.notna().sum())
    QC["n_tool_candidates"] = len(table(w / "features_candidates.csv", ["gene"]))
    if not QC["n_matched"]:
        raise ValueError("No DEG symbols matched DepMap")


def surface(w):
    cfg = C["surface"]
    genes = table(Path(I.deg) / "genes.csv", ["gene", "log2FoldChange", "padj"], "gene", False)
    if cfg["mode"] == "compute":
        module = tool_module(I.tool)
        # Import entry point avoids the original script's automatic pip/venv setup;
        # atlases are read in place from surface.atlas_dir.
        exclude = cfg.get("exclude_gtex_tissues")
        module.run_surface(genes.set_index("gene"), POP, I.source, w,
                           gtex_max_tpm=cfg["gtex_max_tpm"], exclude_tissue_keys=exclude,
                           paper_check=cfg.get("paper_check", False))
        excluded = module.HEMATO_KEYS if exclude is None else exclude
    else:
        snapshot = table(I.source, ["gene", "Surface", "GTEx_max_TPM", "Candidate"], "gene", False)
        snapshot[snapshot.gene.isin(genes.gene)].to_csv(w / "DEGs_surface_GTEx.csv", index=False)
        excluded = "as in the reused snapshot"
    d = table(w / "DEGs_surface_GTEx.csv", ["gene", "Surface", "GTEx_max_TPM", "Candidate"], "gene", False)
    if set(d.gene) != set(genes.gene):
        raise ValueError(f"{POP}: surface table does not cover exactly the input DEG genes")
    QC["rows"] = len(d)
    QC["surface_candidates"] = int(d.Candidate.eq("Yes").sum())
    WARNINGS.append(f"GTEx normal-tissue filter excludes tissues: {excluded}; Candidate is not proven tumor specificity or safety.")


def collectri_args():
    c, a = C["collectri"], C["adata"]
    return ["--sample-col", a["sample_column"], "--group-col", a["condition_column"],
        "--group-map", json.dumps(C["pseudobulk_de"]["group_map"]),
        "--celltype-col", a["celltype_column"], "--groups", a["case"], a["control"],
        "--celltypes", *[C["datasets"][p]["celltype"] for p in C["populations"]],
        "--min-cells", c["min_cells"], "--min-samples", c["min_samples"],
        "--min-targets", c["min_targets"], "--batch-size", c["batch_size"]]


def collectri(w):
    c = C["collectri"]
    source = Path(c["reuse_dir"]) if c["mode"] == "reuse" else Path(I.profiles[0])
    info = table(source / "sample_celltype_info.csv", ["group_id", C["adata"]["sample_column"], C["adata"]["celltype_column"], C["adata"]["condition_column"]], "group_id", False)
    ct = DS["celltype"].strip().lower()
    selected = info.loc[info[C["adata"]["celltype_column"]].str.lower().eq(ct)]
    counts = selected.groupby(C["adata"]["condition_column"])[C["adata"]["sample_column"]].nunique().to_dict()
    if not {C["adata"]["case"].lower(), C["adata"]["control"].lower()} <= {str(k).lower() for k in counts}:
        raise ValueError(f"{POP}: missing case/control profiles")
    QC["donors_per_condition"] = counts
    if min(counts.values()) < c["min_samples"]:
        WARNINGS.append(f"{POP}: insufficient biological replicates for condition comparison: {counts}; activity contrasts are exploratory.")
    run(sys.executable, I.tool, "--reuse", "--reuse-dir", source, "--outdir", w,
        "--deg", Path(I.deg) / "genes.csv", "--gene-col", "gene", "--deg-celltype", DS["celltype"],
        "--plot-celltype", DS["celltype"], "--plot-top", c["plot_top"], "--network-targets", c["network_targets"], *collectri_args())
    d = table(w / "candidate_tf_activity.csv", ["gene", "celltype", "delta_activity", "padj"])
    if not d.empty and set(d.celltype.str.lower()) != {ct}:
        raise ValueError("CollecTRI output has the wrong cell type")
    QC["n_tf_candidates"] = len(d)
    WARNINGS.append("Reuse preserves the existing donor-level summary; it does not recalculate its differential p-values. Within-condition ULM p-values are not JMML-vs-Control tests.")
    shutil.copy2(source / "sample_celltype_info.csv", w / "input_sample_celltype_info.csv")


def liana(w):
    c, a = C["liana"], C["adata"]
    filename = "JMML_vs_Control_DEG_receptors.csv"
    if c["mode"] == "reuse":
        shutil.copy2(I.source, w / filename)
    else:
        run(sys.executable, I.tool, "--input", I.source, "--deg", Path(I.deg) / "genes.csv",
            "--outdir", w, "--gene-col", "gene", "--receivers", DS["celltype"],
            "--condition-col", a["condition_column"], "--celltype-col", a["celltype_column"],
            "--counts-layer", a["counts_layer"],
            "--group-map", json.dumps(C["pseudobulk_de"]["group_map"]),
            "--patient-col", a["sample_column"], "--jmml-label", a["case"].lower(), "--control-label", a["control"].lower(),
            "--min-cells", c["min_cells"], "--expr-prop", c["expr_prop"], "--n-perms", c["n_perms"],
            "--n-jobs", snakemake.threads, "--seed", c["seed"], "--resource", c["resource"])
    d = table(w / filename, ["target", "DEG_as_receptor", "delta_lr_probs"])
    targets = set(d.target.dropna().astype(str))
    if targets - {DS["celltype"]}:
        raise ValueError(f"{POP}: LIANA receiver mismatch: {targets}, expected {DS['celltype']!r}")
    QC["n_lr_rows"] = len(d)
    QC["receivers"] = sorted(targets)
    WARNINGS.append("LIANA pools cells by condition; delta_lr_probs and CellChat p-values are not a between-patient differential test.")


def network_api(w):
    c = C[STAGE]
    if STAGE == "opentargets" and c["mode"] == "reuse":
        shutil.copy2(I.source[0], w / "opentargets_features.csv")
    else:
        cache = Path(C["paths"]["cache_root"]) / C["run_name"] / STAGE / POP
        cache.mkdir(parents=True, exist_ok=True)
        (w / "cache").symlink_to(cache.resolve(), target_is_directory=True)
        if STAGE == "opentargets":
            mapping = cache / "gene_mapping.csv"
            if mapping.exists():
                shutil.copy2(mapping, w / "gene_mapping.csv")
            try:
                run(sys.executable, I.tool, "--input", Path(I.deg) / "genes.csv", "--outdir", w,
                    "--workers", snakemake.threads, "--retries", c["retries"], "--timeout", c["timeout"])
            finally:
                if (w / "gene_mapping.csv").exists():
                    shutil.copy2(w / "gene_mapping.csv", mapping)
        else:
            # Original ChEMBL reader expects exactly seven DESeq-style columns.
            cmd = [sys.executable, I.tool, "--deg", Path(I.deg) / "genes.csv", "--out", w,
                "--workers", snakemake.threads, "--min-pchembl", c["min_pchembl"], "--max-genes", c["max_genes"]]
            if c.get("mapping"):
                cmd += ["--mapping", c["mapping"]]
            run(*cmd)
            status = json.loads((w / "query_status.json").read_text())
            # The legacy parser logs query failures but can return exit status 0.
            if any(status.values()):
                dump(cache / "last_failed_query_status.json", status)
                raise RuntimeError(f"Incomplete ChEMBL queries (cached for retry): {status}")
    filename = "opentargets_features.csv" if STAGE == "opentargets" else "chembl_gene_features.csv"
    d = table(w / filename, ["gene"], "gene", False)
    QC["n_genes"] = len(d)
    if STAGE == "opentargets":
        if "ot_request_status" not in d:
            raise ValueError("Open Targets snapshot lacks request status")
        if d.ot_request_status.eq("error").any():
            raise ValueError("Open Targets snapshot contains API errors; retry them explicitly")
        QC["status_counts"] = d.ot_request_status.value_counts().to_dict()


def integrate(w):
    d = table(Path(I.deg) / "genes.csv", ["gene"], "gene", False)
    original = set(d.gene)
    def join(path, prefix, cols=None):
        nonlocal d
        x = table(path, ["gene"], "gene")
        if cols is not None:
            x = x[["gene", *[c for c in cols if c in x]]]
        else:
            x = x[[c for c in x if c == "gene" or c not in DECOLS]]
        x = x.rename(columns={c: prefix + c for c in x if c != "gene"})
        d = d.merge(x, on="gene", how="left", validate="one_to_one")
    join(Path(I.depmap) / "features.csv", "")
    join(Path(I.surface) / "DEGs_surface_GTEx.csv", "surface_")
    if "collectri" in I.keys():
        # CollecTRI names score columns after the lower-cased case/control labels.
        case, control = C["adata"]["case"].strip().lower(), C["adata"]["control"].strip().lower()
        join(Path(I.collectri) / "candidate_tf_activity.csv", "collectri_",
             ["delta_activity", "pvalue", "padj", "n_jmml", "n_control", f"ulm_score_{case}", f"ulm_score_{control}",
              "activity_available", "jmml_higher"])
    for name in ["opentargets", "chembl"]:
        if name in I.keys():
            filename = "opentargets_features.csv" if name == "opentargets" else "chembl_gene_features.csv"
            join(Path(I[name]) / filename, "")
    if "liana" in I.keys():
        lr = table(Path(I.liana) / "JMML_vs_Control_DEG_receptors.csv")
        x = lr.groupby("DEG_as_receptor", as_index=False).agg(
            liana_n_interactions=("target", "size"), liana_n_comparable=("delta_lr_probs", "count"),
            liana_max_delta=("delta_lr_probs", "max"), liana_mean_delta=("delta_lr_probs", "mean"))
        d = d.merge(x.rename(columns={"DEG_as_receptor": "gene"}), on="gene", how="left", validate="one_to_one")
    # Display independent fields; do not fabricate a combined therapeutic score.
    d.insert(0, "population", POP)
    d["depmap_pass_tool_thresholds"] = d.gene.isin(table(Path(I.depmap) / "features_candidates.csv").gene)
    d["depmap_sufficient_models"] = (d.depmap_n_models.ge(C["depmap"]["min_models"]) & d.depmap_n_dependency_models.ge(C["depmap"]["min_models"])) if "depmap_n_models" in d else pd.NA
    if set(d.gene) != original or len(d) != len(original):
        raise ValueError("Integration lost or duplicated input genes")
    d.sort_values(["padj", "gene"], na_position="last").to_csv(w / "target_evidence.csv", index=False)
    sources = {k: json.loads((Path(v) / "manifest.json").read_text()) for k, v in I.items() if (Path(v) / "manifest.json").is_file()}
    for value in sources.values():
        WARNINGS.extend(value["warnings"])
    QC.update({"n_genes": len(d), "depmap_matched": int(d.depmap_strength.notna().sum()),
        "depmap_threshold_candidates": int(d.depmap_pass_tool_thresholds.sum()),
        "collectri_measured": int(d.collectri_delta_activity.notna().sum()) if "collectri_delta_activity" in d else None})
    dump(w / "source_manifests.json", sources)


def report(w):
    tables, rows = [], []
    for p, path in zip(C["populations"], I):
        tables.append(table(Path(path) / "target_evidence.csv"))
        m = json.loads((Path(path) / "manifest.json").read_text())
        rows.append({"population": p, **m["qc"]})
        WARNINGS.extend(m["warnings"])
    pd.concat(tables, ignore_index=True).to_csv(w / "all_target_evidence.csv", index=False)
    pd.DataFrame(rows).to_csv(w / "run_summary.csv", index=False)
    dump(w / "config.resolved.json", C)
    lines = ["# OnlyFasta — результат запуска", "", f"Run: `{C['run_name']}`", "",
             "## Выполненные блоки", "", "| Блок | Режим |", "|---|---|"]
    lines += [f"| {k} | {C[k]['mode']} |" for k in ["depmap", "surface", "collectri", "liana", "opentargets", "chembl"]]
    lines += ["", "`reuse` означает использование существующего снимка. DepMap пересекается с текущими DEG; CollecTRI заново вычисляет ULM групповых профилей и выбор TF, сохраняя ранее рассчитанный донорский контраст.", "",
              "## Таблицы", "", "- `all_target_evidence.csv` — одна строка на популяцию и ген; пропуски не заменяются нулями.",
              "- `run_summary.csv` — покрытие входных генов.",
              "- `config.resolved.json` — параметры запуска; manifest.json каждого этапа содержит версии и контрольные суммы.",
              "", "Порядок строк определяется исходным DE padj, это не рейтинг пригодности для терапии. Сумма баллов не рассчитывается: зависимость, активность TF, поверхностная локализация и tractability являются разными критериями.", "", "## Ограничения", ""]
    lines += ["- " + x for x in sorted(set(WARNINGS))]
    lines += ["", "Технически успешный workflow не подтверждает корректность исходного DE, состава пациентов или выбора контроля. Поле Candidate старого surface-инструмента не доказывает селективность в отношении нормальных HSPC.", ""]
    (w / "REPORT.md").write_text("\n".join(lines))
    QC["total_gene_population_rows"] = sum(len(t) for t in tables)


def summary(w):
    module = tool_module(I.tool)
    table = module.build_summary(Path(I.report).parent, C["populations"], C["adata"]["case"],
                                 C["adata"]["control"], w,
                                 report_genes=C.get("summary", {}).get("opentargets_report_genes", 2),
                                 ot_report_path=I.ot_report)
    QC["n_candidate_rows"] = len(table)
    QC["per_population"] = table.Population.value_counts().to_dict()
    QC["tf_branch"] = int(table.TF.eq("Yes").sum())


def execute(w):
    if STAGE == "fetch_surface":
        module = tool_module(I.tool)
        module.fetch_data(w, update_hpa=False)
        QC["source_urls"] = module.URLS
    elif STAGE == "collectri_profiles":
        run(sys.executable, I.tool, "--input", I.adata, "--outdir", w, "--backed", "--no-plots",
            "--counts-layer", C["adata"]["counts_layer"], *collectri_args())
        table(w / "tf_activity_summary.csv", ["gene", "celltype", "padj"], allow_empty=False)
    else:
        funcs = {"normalize": norm_deg, "depmap_features": depmap_features, "depmap": depmap,
            "surface": surface, "collectri": collectri, "liana": liana,
            "opentargets": network_api, "chembl": network_api, "integrate": integrate,
            "report": report, "summary": summary, "pseudobulk_de": pseudobulk_de}
        funcs[STAGE](w)


out = Path(str(snakemake.output[0]))
out.parent.mkdir(parents=True, exist_ok=True)
log = Path(str(snakemake.log[0]))
log.parent.mkdir(parents=True, exist_ok=True)
with log.open("w", buffering=1) as stream, contextlib.redirect_stdout(stream), contextlib.redirect_stderr(stream):
    print(f"START {STAGE} {POP or ''}", flush=True)
    try:
        with tempfile.TemporaryDirectory(prefix=f".{out.name}.partial-", dir=out.parent) as tmp:
            work = Path(tmp)
            execute(work)
            manifest = {"stage": STAGE, "population": POP, "completed_utc": datetime.now(timezone.utc).isoformat(),
                "config": C, "inputs": [fingerprint(p) for p in I], "qc": QC, "warnings": sorted(set(WARNINGS)),
                "python": sys.version, "packages": {d.metadata["Name"]: d.version for d in importlib.metadata.distributions()}}
            dump(work / "manifest.json", manifest)
            # Snakemake owns this declared output directory, never an input path.
            if out.exists():
                shutil.rmtree(out)
            work.replace(out)
        print("COMPLETE", QC, flush=True)
    except BaseException:
        traceback.print_exc()
        raise
