#!/usr/bin/env python3
r"""DEG × атласы поверхностных белков (HPA + CSPA + TCSA) × GTEx v11.

Одна популяция за запуск:
  python3 surface_gtex_pipeline.py --deg DEGs.csv --name HSC --atlas-dir resources/surface_gtex \
      --outdir results/surface/HSC [--gene-col gene] [--gtex-max-tpm 20] \
      [--exclude-tissues "whole blood" spleen] [--paper-check] [--no-venv]
Только скачать атласы:  python3 surface_gtex_pipeline.py --atlas-dir DIR --fetch-only
Python API:  run_surface(deg_df, "HSC", atlas_dir, outdir, gtex_max_tpm=20, exclude_tissue_keys=[...])
"""
import argparse
import importlib.util
import os
import shutil
import subprocess
import sys
import venv
import zipfile
from pathlib import Path

REQS = ["pandas", "numpy", "openpyxl", "certifi"]
HERE = Path(__file__).resolve().parent
VENV_DIR = HERE / ".venv"

URLS = {
    "hpa_zip": "https://www.proteinatlas.org/download/proteinatlas.tsv.zip",
    "cspa": "https://wollscheidlab.org/CSPA/data/S2_File.xlsx",
    "tcsa": "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs43018-021-00282-w/MediaObjects/43018_2021_282_MOESM2_ESM.xlsx",
    "gtex_v11": "https://storage.googleapis.com/adult-gtex/bulk-gex/v11/rna-seq/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz",
    "gtex_v8": "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
}

HPA_FILE = "proteinatlas.tsv"
CSPA_FILE = "S2_File.xlsx"
TCSA_FILE = "43018_2021_282_MOESM2_ESM.xlsx"
GTEX_FILE = "GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz"
GTEX_V8_FILE = "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"  # только для сопоставления названий
GTEX_VERSION = "v11"

HPA_LOCATIONS = {"Plasma membrane", "Cell Junctions"}
HPA_CLASS_CD = "CD markers"
HPA_SECRETOME = "Intracellular and membrane"
CSPA_CATEGORIES = ["1 - high confidence", "2 - putative", "3 - unspecific"]
TCSA_EXCLUDE_MONOTOPIC = False
GTEX_MAX_TPM = 20
# Ткани GTEx, которые не учитываются при оценке экспрессии в норме. Значение по умолчанию —
# кроветворные ткани (для гемобластозов вроде JMML); для солидных опухолей обычно [].
HEMATO_KEYS = ["whole blood", "spleen", "ebv-transformed lymphocytes"]

# Сверка с JMML-статьёй Werner et al. 2025; только при --paper-check (surface.paper_check).
PAPER = {
    "bulk, 15 финальных": ["CLEC12A", "CD300LB", "CEACAM8", "CLEC5A", "ENTPD7", "FCAR", "FCGR1A", "ITGB7",
                           "LILRB1", "LRRC25", "PTAFR", "SIGLEC9", "SIRPB1", "SLC22A4", "SLC8A1"],
    "scRNA HSC (рис. 2D)": ["TNF", "CD7", "CD96", "ADGRE2", "SAMSN1", "SLC8A1", "CLEC12A", "ITGAL", "RELT",
                            "GPR132", "TNFRSF18", "ICOSLG", "CD70", "IL2RA", "CRLF2", "TNFRSF4", "IL1RAP",
                            "GP6", "GPR174", "UBASH3B", "CLECL1", "RASGRP1", "HAVCR2", "CD86", "CSF2RA",
                            "HLA-DQA2", "CD48", "CD79A", "LTB"],
    "РНК + белок": ["CLEC12A", "CD48", "ADGRE2"],
    "упомянуты отдельно": ["CD38", "CD33", "IL3RA"],
}
CHECK_GENES = ["CLEC12A", "CD48", "ADGRE2", "CD33", "IL3RA", "CD7", "CD96", "CD38"]
BOOL_COLS = ["HPA_PM_or_CJ", "HPA_CD_marker", "HPA_secretome_IM", "HPA_surface",
             "CSPA_surface", "TCSA_surface", "TCSA_monotopic"]


# ---------- окружение ----------
def venv_python():
    return VENV_DIR / ("Scripts/python.exe" if os.name == "nt" else "bin/python")


def ensure_env(use_venv):
    if not use_venv:
        missing = [m for m in REQS if importlib.util.find_spec(m) is None]
        if missing:
            print("Устанавливаю пакеты:", ", ".join(missing), flush=True)
            subprocess.check_call([sys.executable, "-m", "pip", "install", "-q",
                                   "--disable-pip-version-check", *missing])
        return
    if Path(sys.prefix).resolve() == VENV_DIR.resolve():
        return
    py = venv_python()
    if not py.exists():
        print(f"Создаю окружение {VENV_DIR}", flush=True)
        try:
            venv.create(VENV_DIR, with_pip=True)
        except Exception as e:
            sys.exit(f"Не удалось создать venv ({e}). Запустите с --no-venv.")
    check = subprocess.run([str(py), "-c", "import " + ", ".join(REQS)], capture_output=True)
    if check.returncode:
        print("Устанавливаю пакеты:", ", ".join(REQS), flush=True)
        subprocess.check_call([str(py), "-m", "pip", "install", "-q", "--disable-pip-version-check", *REQS])
    os.execv(str(py), [str(py), str(Path(__file__).resolve()), *sys.argv[1:]])


def download(url, path, force=False):
    if path.exists() and not force:
        return
    import ssl
    import urllib.request
    import certifi
    print(f"Скачиваю {path.name}", flush=True)
    ctx = ssl.create_default_context(cafile=certifi.where())
    tmp = path.with_suffix(path.suffix + ".part")
    with urllib.request.urlopen(url, context=ctx, timeout=300) as r, open(tmp, "wb") as f:
        shutil.copyfileobj(r, f)
    tmp.replace(path)


def fetch_data(d, update_hpa):
    hpa_zip = d / "proteinatlas.tsv.zip"
    download(URLS["hpa_zip"], hpa_zip, force=update_hpa)
    if update_hpa or not (d / HPA_FILE).exists():
        with zipfile.ZipFile(hpa_zip) as z:
            z.extract(HPA_FILE, d)
    download(URLS["cspa"], d / CSPA_FILE)
    download(URLS["tcsa"], d / TCSA_FILE)
    download(URLS["gtex_v11"], d / GTEX_FILE)
    download(URLS["gtex_v8"], d / GTEX_V8_FILE)


# ---------- анализ ----------
def read_deg(path, gene_col=None):
    """DEG-таблица с генами в колонке gene_col (по умолчанию — первая колонка)."""
    import pandas as pd
    deg = pd.read_csv(path, index_col=0) if gene_col is None else pd.read_csv(path).set_index(gene_col)
    deg.index = deg.index.astype(str)
    deg.index.name = "gene"
    return deg


def run_surface(deg, name, atlas_dir, outdir, *, gtex_max_tpm=GTEX_MAX_TPM,
                exclude_tissue_keys=None, paper_check=False):
    """Аннотирует DEG одной популяции и пишет таблицы в outdir.

    deg: DataFrame с генами в индексе и колонками DESeq2 (нужны log2FoldChange и padj).
    Выходы: DEGs_surface.csv, DEGs_surface_GTEx.csv, candidates.csv, DEGs_surface_GTEx.xlsx
    (+ paper_comparison.csv при paper_check).
    """
    import warnings
    import numpy as np
    import pandas as pd

    d, out = Path(atlas_dir), Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    exclude = [k.strip().lower() for k in (HEMATO_KEYS if exclude_tissue_keys is None else exclude_tissue_keys)]
    missing_cols = {"log2FoldChange", "padj"} - set(deg.columns)
    if missing_cols:
        sys.exit(f"{name}: в DEG-таблице нет колонок {sorted(missing_cols)}")
    missing = [f for f in (HPA_FILE, CSPA_FILE, TCSA_FILE, GTEX_FILE) if not (d / f).exists()]
    if missing:
        sys.exit(f"Нет атласов в {d}: {', '.join(missing)} (скачайте: --fetch-only)")

    warnings.filterwarnings("ignore", message="Unknown extension is not supported")
    pd.set_option("display.width", 250, "display.max_columns", 50, "display.max_colwidth", 120)

    def show(title, obj=None):
        print(f"\n=== {title} ===")
        if obj is not None:
            print(obj.to_string() if hasattr(obj, "to_string") else obj)

    degs = {name: deg}

    overview = pd.DataFrame({
        ct: {"DEG всего": len(x),
             "↑ повышены": int((x.log2FoldChange > 0).sum()),
             "↓ понижены": int((x.log2FoldChange < 0).sum())}
        for ct, x in degs.items()}).T
    show("Входные DEG (padj < 0.05, |log2FC| ≥ 1)", overview)

    # HPA
    hpa = pd.read_csv(d / HPA_FILE, sep="\t", low_memory=False)
    hpa["PM_or_CJ"] = hpa["Subcellular location"].fillna("").str.split(",").map(
        lambda ts: any(t.strip() in HPA_LOCATIONS for t in ts))
    hpa["CD_marker"] = hpa["Protein class"].fillna("").str.split(", ").map(lambda ts: HPA_CLASS_CD in ts)
    hpa["Secretome_IM"] = hpa["Secretome location"].fillna("").eq(HPA_SECRETOME)
    hpa_flags = hpa.groupby("Gene")[["PM_or_CJ", "CD_marker", "Secretome_IM"]].any()
    hpa_flags["any"] = hpa_flags.any(axis=1)
    hpa_surface = set(hpa_flags.index[hpa_flags["any"]])

    hpa_current = set(hpa["Gene"])
    syn2gene = {}
    for gene, syns in hpa[["Gene", "Gene synonym"]].dropna().itertuples(index=False):
        for s in syns.split(", "):
            syn2gene.setdefault(s.strip(), set()).add(gene)

    def to_hpa(g):
        if g in hpa_current:
            return g
        hits = syn2gene.get(g, set())
        return next(iter(hits)) if len(hits) == 1 else None

    # CSPA
    cspa = pd.read_excel(d / CSPA_FILE, sheet_name="Table A")
    cspa = cspa[(cspa["organism"] == "Human") & cspa["CSPA category"].isin(CSPA_CATEGORIES)].copy()
    cspa["symbol"] = cspa["ENTREZ gene symbol"].astype(str).str.strip()
    cspa = cspa[~cspa["symbol"].isin(["0", "nan", ""])]
    cspa_cat = cspa.groupby("symbol")["CSPA category"].min()
    cspa_surface = set(cspa_cat.index)

    # TCSA
    tcsa = pd.read_excel(d / TCSA_FILE, sheet_name="Table S2", header=2)
    tcsa.columns = [str(c).strip() for c in tcsa.columns]
    tcsa = tcsa[tcsa["ENSEMBL Gene ID"].astype(str).str.startswith("ENSG")].copy()
    tcsa["symbol"] = tcsa["HGNC Symbol"].astype(str).str.strip()
    tcsa_mono = tcsa.groupby("symbol")["Integral monotopic GESP"].agg(lambda s: (s == "Y").any())
    if TCSA_EXCLUDE_MONOTOPIC:
        tcsa_mono = tcsa_mono[~tcsa_mono]
    tcsa_surface = set(tcsa_mono.index)

    def annotate_surface(genes):
        a = pd.DataFrame(index=pd.Index(list(genes), name="gene"))
        a["HPA_symbol"] = [to_hpa(g) for g in a.index]
        f = hpa_flags.reindex(a["HPA_symbol"].tolist())
        a["HPA_PM_or_CJ"] = f["PM_or_CJ"].fillna(False).astype(bool).values
        a["HPA_CD_marker"] = f["CD_marker"].fillna(False).astype(bool).values
        a["HPA_secretome_IM"] = f["Secretome_IM"].fillna(False).astype(bool).values
        a["HPA_surface"] = a[["HPA_PM_or_CJ", "HPA_CD_marker", "HPA_secretome_IM"]].any(axis=1)

        def pick(series, g, h):
            if g in series.index:
                return series[g]
            if h is not None and h in series.index:
                return series[h]
            return np.nan

        a["CSPA_category"] = [pick(cspa_cat, g, h) for g, h in zip(a.index, a["HPA_symbol"])]
        a["CSPA_surface"] = a["CSPA_category"].notna()
        mono = [pick(tcsa_mono, g, h) for g, h in zip(a.index, a["HPA_symbol"])]
        a["TCSA_surface"] = pd.notna(mono)
        a["TCSA_monotopic"] = [bool(m) if pd.notna(m) else False for m in mono]
        a["n_atlases"] = a[["HPA_surface", "CSPA_surface", "TCSA_surface"]].sum(axis=1)
        a["Surface"] = np.where(a["n_atlases"] > 0, "Membrane", "Not_surface")
        return a

    def which_atlases(row):
        hpa_why = [n for n, c in [("PM/CJ", "HPA_PM_or_CJ"), ("CD", "HPA_CD_marker"),
                                  ("secretome", "HPA_secretome_IM")] if row[c]]
        parts = []
        if hpa_why:
            parts.append("HPA (" + ", ".join(hpa_why) + ")")
        if row["CSPA_surface"]:
            parts.append("CSPA " + str(row["CSPA_category"])[0])
        if row["TCSA_surface"]:
            parts.append("TCSA" + (" (monotopic)" if row["TCSA_monotopic"] else ""))
        return "; ".join(parts) or "—"

    def to_yes_no(df):
        out = df.copy()
        for col in BOOL_COLS:
            out[col] = out[col].map({True: "Yes", False: "No"})
        return out

    show("Атласы: сколько генов считаются мембранными", pd.Series({
        "HPA: Plasma membrane / Cell Junctions": int(hpa_flags["PM_or_CJ"].sum()),
        "HPA: CD markers": int(hpa_flags["CD_marker"].sum()),
        "HPA: Secretome — Intracellular and membrane": int(hpa_flags["Secretome_IM"].sum()),
        "HPA итого": len(hpa_surface),
        "CSPA": len(cspa_surface),
        "TCSA": len(tcsa_surface),
        "хотя бы в одном атласе": len(hpa_surface | cspa_surface | tcsa_surface),
    }, name="генов"))

    if paper_check:
        ctrl = annotate_surface(CHECK_GENES)
        mark = {True: "✓", False: "✗"}
        show("Контроль на мишенях из статьи", pd.DataFrame({
            "HPA": ctrl["HPA_surface"].map(mark), "CSPA": ctrl["CSPA_surface"].map(mark),
            "TCSA": ctrl["TCSA_surface"].map(mark), "итог": ctrl["Surface"],
            "основание": ctrl.apply(which_atlases, axis=1)}))

    # воронка по атласам
    funnel = []
    for ct, deg in degs.items():
        t = deg[deg.log2FoldChange > 0].join(annotate_surface(deg.index[deg.log2FoldChange > 0]))
        n_memb = int((t["Surface"] == "Membrane").sum())
        funnel.append({
            "тип клеток": ct, "повышены": len(t),
            "в HPA": int(t["HPA_surface"].sum()), "в CSPA": int(t["CSPA_surface"].sum()),
            "в TCSA": int(t["TCSA_surface"].sum()), "≥1 атлас": n_memb,
            "%": round(100 * n_memb / len(t), 1) if len(t) else 0.0,
            "в 1 / 2 / 3": f"{(t.n_atlases == 1).sum()} / {(t.n_atlases == 2).sum()} / {(t.n_atlases == 3).sum()}"})
    show("Повышенные гены → атласы", pd.DataFrame(funnel).set_index("тип клеток"))

    # GTEx v11
    gtex = pd.read_csv(d / GTEX_FILE, sep="\t", skiprows=2)
    norm = lambda s: s.lower().replace("_", " ")
    ens = lambda s: s.astype(str).str.split(".").str[0]
    hemato = [c for c in gtex.columns if any(k in norm(c) for k in exclude)]
    tissues = [c for c in gtex.columns if c not in ["Name", "Description", *hemato]]
    by_sym = gtex.groupby("Description")[tissues].max()
    by_ens = gtex.assign(ens=ens(gtex["Name"])).groupby("ens")[tissues].max()
    sym_max, sym_top = by_sym.max(axis=1), by_sym.idxmax(axis=1)
    ens_max, ens_top = by_ens.max(axis=1), by_ens.idxmax(axis=1)

    if (d / GTEX_V8_FILE).exists():
        v8 = pd.read_csv(d / GTEX_V8_FILE, sep="\t", skiprows=2, usecols=["Name", "Description"])
        old_sym2ens = v8.assign(ens=ens(v8["Name"])).drop_duplicates("Description").set_index("Description")["ens"]
    else:
        old_sym2ens = pd.Series(dtype=str)
        print("GTEx v8 не найден — сопоставление переименованных генов через Ensembl ID пропущено")

    def gtex_lookup(gene, hpa_symbol):
        for how, key in [("symbol", gene), ("HPA_symbol", hpa_symbol)]:
            if isinstance(key, str) and key in sym_max.index:
                return sym_max[key], sym_top[key], how
        e = old_sym2ens.get(gene)
        if e is not None and e in ens_max.index:
            return ens_max[e], ens_top[e], "Ensembl"
        return np.nan, np.nan, "not_found"

    show(f"GTEx {GTEX_VERSION}: {len(tissues)} тканей; исключены: {', '.join(hemato) or '—'}")

    # аннотация всех DEG
    surface_only, final, cands, funnel = {}, {}, {}, []
    for ct, deg in degs.items():
        t = deg.join(annotate_surface(deg.index))
        t.insert(len(deg.columns), "DE_status", np.where(t.log2FoldChange > 0, "Up", "Down"))
        t.index.name = "gene"
        surface_only[ct] = t.copy()

        res = [gtex_lookup(g, h) for g, h in zip(t.index, t["HPA_symbol"])]
        tpm = pd.Series([r[0] for r in res], index=t.index, dtype=float)
        t["GTEx_max_TPM"] = tpm.round(2)
        t["GTEx_max_tissue"] = [r[1] for r in res]
        t["GTEx_match"] = [r[2] for r in res]
        t["Tumor_specific"] = np.select([tpm.isna(), tpm < gtex_max_tpm],
                                        ["No_GTEx_data", "Specific"], default="Not_specific")
        t["GTEx_version"] = GTEX_VERSION
        t["Candidate"] = np.where((t.DE_status == "Up") & (t.Surface == "Membrane")
                                  & (t.Tumor_specific == "Specific"), "Yes", "No")
        t["Atlases"] = t.apply(which_atlases, axis=1)

        um = t[(t.DE_status == "Up") & (t.Surface == "Membrane")]
        c = t[t.Candidate == "Yes"].sort_values("log2FoldChange", ascending=False)
        final[ct], cands[ct] = t, c
        funnel.append({
            "тип клеток": ct,
            "повышены": int((t.DE_status == "Up").sum()),
            "+ мембранные": len(um),
            f"+ GTEx < {gtex_max_tpm} TPM → кандидаты": len(c),
            "мембранные: Not_specific": int((um.Tumor_specific == "Not_specific").sum()),
            "мембранные: No_GTEx_data": int((um.Tumor_specific == "No_GTEx_data").sum())})

    # сохранение: версия GTEx записана в колонке GTEx_version, не в имени файла
    to_yes_no(surface_only[name]).to_csv(out / "DEGs_surface.csv")
    to_yes_no(final[name]).to_csv(out / "DEGs_surface_GTEx.csv")
    to_yes_no(cands[name]).to_csv(out / "candidates.csv")
    with pd.ExcelWriter(out / "DEGs_surface_GTEx.xlsx") as xw:
        to_yes_no(surface_only[name]).to_excel(xw, sheet_name="surface")
        to_yes_no(final[name]).to_excel(xw, sheet_name="surface_GTEx")
        to_yes_no(cands[name]).to_excel(xw, sheet_name="candidates")

    show(f"Воронка: атласы → GTEx {GTEX_VERSION}", pd.DataFrame(funnel).set_index("тип клеток"))

    for ct, c in cands.items():
        show(f"{ct}: прошли фильтр GTEx — {len(c)} генов", ", ".join(c.index) or "—")
        print(pd.DataFrame({
            "log2FC": c["log2FoldChange"].round(2),
            "padj": c["padj"].map(lambda x: f"{x:.1e}"),
            "атласы": c["Atlases"],
            "TPM": c["GTEx_max_TPM"].round(1),
            "ткань": c["GTEx_max_tissue"]}).to_string())
        nodata = final[ct][(final[ct].DE_status == "Up") & (final[ct].Surface == "Membrane")
                           & (final[ct].Tumor_specific == "No_GTEx_data")]
        if len(nodata):
            print(f"Без данных GTEx (проверить вручную): {', '.join(nodata.index)}")

    all_genes = sorted(set().union(*[set(c.index) for c in cands.values()]))
    summary = pd.DataFrame(
        {ct: [f"↑ {cands[ct].loc[g, 'log2FoldChange']:.1f}" if g in cands[ct].index else "—" for g in all_genes]
         for ct in cands}, index=pd.Index(all_genes, name="ген"))
    first = {g: next(final[ct].loc[g] for ct in cands if g in cands[ct].index) for g in all_genes}
    summary["TPM"] = [round(first[g]["GTEx_max_TPM"], 1) for g in all_genes]
    summary["ткань"] = [first[g]["GTEx_max_tissue"] for g in all_genes]
    summary["атласы"] = [first[g]["Atlases"] for g in all_genes]
    show(f"Все типы клеток: {len(all_genes)} уникальных генов прошли фильтр GTEx {GTEX_VERSION}",
         ", ".join(all_genes))
    print(summary.to_string())

    if not paper_check:
        print(f"\nГотово. Файлы сохранены в {out}")
        return final[name]

    # сверка со статьёй
    rows = []
    for lst, genes in PAPER.items():
        for g in genes:
            rows.append({"список": lst, "ген": g,
                         **{ct: deg.loc[g, "log2FoldChange"] if g in deg.index else np.nan
                            for ct, deg in degs.items()}})
    cmp = pd.DataFrame(rows)
    cmp["повышен у нас"] = cmp[list(degs)].gt(0).any(axis=1)
    surf = annotate_surface(cmp["ген"].unique())
    cmp["атласы"] = cmp["ген"].map(surf.apply(which_atlases, axis=1))
    cmp.to_csv(out / "paper_comparison.csv", index=False)

    summ = []
    for lst, s in cmp.groupby("список", sort=False):
        found = s.loc[s["повышен у нас"], "ген"].tolist()
        summ.append({"список": lst, "генов": len(s), "повышены у нас": len(found),
                     "какие": ", ".join(found) or "—"})
    show("Сверка с Werner et al. 2025", pd.DataFrame(summ).set_index("список"))

    fmt = lambda x: "—" if pd.isna(x) else (f"↑ {x:.1f}" if x > 0 else f"↓ {abs(x):.1f}")
    hits = cmp[cmp[list(degs)].notna().any(axis=1)].drop_duplicates("ген").set_index("ген")
    show("Гены из статьи, значимые у нас", hits[list(degs)].map(fmt).join(hits["атласы"]))

    print(f"\nГотово. Файлы сохранены в {out}")
    return final[name]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--deg", type=Path, help="DEG-таблица одной популяции (CSV)")
    ap.add_argument("--name", help="имя популяции (подписи в выводе)")
    ap.add_argument("--gene-col", default=None, help="колонка с генами; по умолчанию первая колонка")
    ap.add_argument("--atlas-dir", type=Path, required=True, help="папка с атласами (скачиваются, если их нет)")
    ap.add_argument("--outdir", type=Path, help="папка для результатов")
    ap.add_argument("--gtex-max-tpm", type=float, default=GTEX_MAX_TPM)
    ap.add_argument("--exclude-tissues", nargs="*", default=HEMATO_KEYS,
                    help="подстроки названий тканей GTEx, которые не учитываются; без значений — не исключать ничего")
    ap.add_argument("--paper-check", action="store_true", help="сверка с JMML-статьёй Werner et al. 2025")
    ap.add_argument("--fetch-only", action="store_true", help="только скачать атласы в --atlas-dir")
    ap.add_argument("--no-venv", action="store_true", help="не создавать .venv, ставить пакеты в текущий Python")
    ap.add_argument("--update-hpa", action="store_true", help="перекачать Human Protein Atlas")
    args = ap.parse_args()
    if not args.fetch_only and not (args.deg and args.name and args.outdir):
        ap.error("нужны --deg, --name и --outdir (или --fetch-only)")

    ensure_env(use_venv=not args.no_venv)
    atlas_dir = args.atlas_dir.resolve()
    atlas_dir.mkdir(parents=True, exist_ok=True)
    fetch_data(atlas_dir, args.update_hpa)
    if not args.fetch_only:
        run_surface(read_deg(args.deg, args.gene_col), args.name, atlas_dir, args.outdir.resolve(),
                    gtex_max_tpm=args.gtex_max_tpm, exclude_tissue_keys=args.exclude_tissues,
                    paper_check=args.paper_check)
