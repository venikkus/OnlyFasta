#!/usr/bin/env python3
"""JMML DEG -> ChEMBL evidence, with limited I/O concurrency and resumable JSON caches.

Requires: requests, pandas, mygene. Uses official ChEMBL REST JSON endpoints.
"""
import argparse
import concurrent.futures as cf
import gzip
import hashlib
import json
import logging
import os
import re
import threading
import time
from pathlib import Path
from urllib.parse import urljoin

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

BASE = "https://www.ebi.ac.uk/chembl/api/data/"
_tls = threading.local()


def session():
    if not hasattr(_tls, "session"):
        s = requests.Session()
        s.headers.update({"Accept": "application/json", "User-Agent": "jmml-target-ranking/1.0"})
        retry = Retry(total=4, connect=4, read=3, status=4,
                      backoff_factor=1.2, status_forcelist=[429, 500, 502, 503, 504],
                      allowed_methods=frozenset(["GET"]))
        adapter = HTTPAdapter(max_retries=retry, pool_connections=4, pool_maxsize=4)
        s.mount("https://", adapter)
        _tls.session = s
    return _tls.session


def json_get(url, params=None, timeout=55):
    r = session().get(url, params=params, timeout=(15, timeout))
    r.raise_for_status()
    return r.json()


def pages(endpoint, filters, key, page_size=500, max_pages=1000):
    url = BASE + endpoint + ".json"
    params = dict(filters, limit=page_size, offset=0)
    rows = []
    for _ in range(max_pages):
        obj = json_get(url, params=params)
        rows.extend(obj.get(key, []))
        nxt = obj.get("page_meta", {}).get("next")
        if not nxt:
            return rows
        url = urljoin(BASE, nxt)
        params = None
    raise RuntimeError(f"Exceeded {max_pages} pages for {endpoint} query")


def read_cache(path):
    if not path.exists():
        return None
    with gzip.open(path, "rt", encoding="utf8") as f:
        return json.load(f)


def write_cache(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".{os.getpid()}.{threading.get_ident()}.tmp")
    with gzip.open(tmp, "wt", encoding="utf8") as f:
        json.dump(payload, f, ensure_ascii=False)
    os.replace(tmp, path)


def cached(path, fn):
    val = read_cache(path)
    if val is not None:
        return val
    val = fn()
    write_cache(path, val)  # [] is a successful empty response
    return val


def parallel(items, fn, workers, label):
    results, failed = {}, {}
    with cf.ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(fn, x): x for x in items}
        for i, f in enumerate(cf.as_completed(futures), 1):
            x = futures[f]
            try:
                results[x] = f.result()
            except Exception as e:
                failed[x] = str(e)
                logging.warning("%s failed for %s: %s", label, x, e)
            if i % 10 == 0 or i == len(items):
                logging.info("%s: %s/%s completed (%s errors)", label, i, len(items), len(failed))
    return results, failed


def batched(seq, n):
    seq = list(seq)
    return [seq[i:i+n] for i in range(0, len(seq), n)]


def target_lookup(accession, root):
    path = root / "target" / (accession + ".json.gz")
    def run():
        rows = pages("target", {
            "target_components__accession": accession,
            "target_type": "SINGLE PROTEIN", "organism": "Homo sapiens"
        }, "targets")
        # Exact UniProt component match, not just fuzzy search result
        ids = []
        for t in rows:
            components = t.get("target_components") or []
            if accession in [c.get("accession") for c in components]:
                ids.append(t["target_chembl_id"])
        return sorted(set(ids))
    return cached(path, run)


def get_activities(target, root):
    path = root / "activities" / (target + ".json.gz")
    def run():
        rows = pages("activity", {
            "target_chembl_id": target,
            "standard_type__in": "IC50,Ki,Kd",
            "standard_relation": "=",
            "pchembl_value__isnull": "false"
        }, "activities")
        # Strictly recheck server-side filters (API version differences)
        return [a for a in rows
                if a.get("target_chembl_id") == target
                and a.get("standard_type") in ("IC50", "Ki", "Kd")
                and a.get("standard_relation") == "="
                and a.get("pchembl_value") is not None]
    return cached(path, run)


def bulk_metadata(kind, ids, root, batch_size=30):
    key = {"assay": "assays", "molecule": "molecules", "mechanism": "mechanisms"}[kind]
    field = {"assay": "assay_chembl_id", "molecule": "molecule_chembl_id",
             "mechanism": "molecule_chembl_id"}[kind]
    batches = batched(sorted(set(ids)), batch_size)
    jobs = [(i, tuple(batch)) for i, batch in enumerate(batches)]
    def work(job):
        _, batch = job
        digest = hashlib.sha256("|".join(batch).encode()).hexdigest()[:24]
        path = root / kind / (digest + ".json.gz")
        return cached(path, lambda: pages(kind, {field + "__in": ",".join(batch)}, key))
    return jobs, work


def numeric(s):
    return pd.to_numeric(s, errors="coerce")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--deg", required=True, help="CSV/TSV with column gene")
    ap.add_argument("--out", default="chembl_results")
    ap.add_argument("--workers", type=int, default=8, help="concurrent API requests; suggested 2-4")
    ap.add_argument("--min-pchembl", type=float, default=6)
    ap.add_argument("--mapping", help="Optional CSV with columns gene,target_chembl_id (skip MyGene lookup)")
    ap.add_argument("--max-genes", type=int, default=0, help="testing only; 0 = all genes")
    args = ap.parse_args()
    if not 1 <= args.workers <= 8:
        ap.error("--workers must be 1..8 to avoid overloading the API")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    root = out / "cache"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(), logging.FileHandler(out / "pipeline.log")])

    deg = pd.read_csv(args.deg, sep=None, engine="python", names=['gene', 'baseMean', 'log2FC', 'lfcSE', 'stat', 'pvalue', 'padj'], header=0)
    if "gene" not in deg:
        raise ValueError("Input must have a column named 'gene'")
    genes = sorted(set(deg.gene.dropna().astype(str).str.strip()) - {""})
    if args.max_genes:
        genes = genes[:args.max_genes]
    logging.info("Input genes: %s", len(genes))

    # Mapping stage: gene -> reviewed UniProt accession -> exact ChEMBL SINGLE PROTEIN target
    mapping_rows = []
    mapping_failed = set()
    if args.mapping:
        tab = pd.read_csv(args.mapping)
        if not {"gene", "target_chembl_id"} <= set(tab):
            raise ValueError("--mapping requires columns gene,target_chembl_id")
        for row in tab.itertuples(index=False):
            if str(row.gene) in genes and pd.notna(row.target_chembl_id):
                mapping_rows.append({"gene": str(row.gene), "uniprot": None,
                                     "target_chembl_id": str(row.target_chembl_id), "mapping_status": "mapped"})
    else:
        import mygene
        mg_file = root / "mygene.json.gz"
        def get_mygene():
            mg = mygene.MyGeneInfo()
            return mg.querymany(genes, scopes="symbol", fields="symbol,taxid,uniprot.Swiss-Prot",
                                species="human", as_dataframe=False, verbose=False)
        # Cache depends on gene list to prevent stale mappings
        digest = hashlib.sha256("|".join(genes).encode()).hexdigest()[:20]
        entries = cached(root / f"mygene_{digest}.json.gz", get_mygene)
        accessions = {}
        for h in entries:
            if h.get("notfound") or h.get("taxid") != 9606:
                continue
            q = h.get("query")
            if q not in genes or h.get("symbol", "").upper() != q.upper():
                continue
            acc = (h.get("uniprot") or {}).get("Swiss-Prot", [])
            if isinstance(acc, str): acc = [acc]
            accessions.setdefault(q, set()).update(acc)
        all_accessions = sorted(set().union(*accessions.values())) if accessions else []
        acc_to_targets, acc_failed = parallel(all_accessions,
            lambda x: target_lookup(x, root), args.workers, "target mapping")
        for gene in genes:
            accs = sorted(accessions.get(gene, []))
            if any(a in acc_failed for a in accs):
                mapping_failed.add(gene)
            tids = sorted({tid for a in accs for tid in acc_to_targets.get(a, [])})
            for tid in tids:
                mapping_rows.append({"gene": gene, "uniprot": ";".join(accs),
                                     "target_chembl_id": tid, "mapping_status": "mapped"})
    target_map = pd.DataFrame(mapping_rows, columns=["gene","uniprot","target_chembl_id","mapping_status"])
    for gene in genes:
        if gene not in set(target_map.gene):
            target_map.loc[len(target_map)] = [gene, None, None,
               "query_failed" if gene in mapping_failed else "not_mapped"]
    target_map.to_csv(out / "chembl_target_mapping.csv", index=False)
    mapped = target_map.dropna(subset=["target_chembl_id"])
    tids = sorted(mapped.target_chembl_id.unique())
    logging.info("Mapped genes %s, targets %s", mapped.gene.nunique(), len(tids))

    activity_map, failed_tids = parallel(tids, lambda t: get_activities(t, root), args.workers, "activities")
    acts = [a for rows in activity_map.values() for a in rows]
    activity = pd.DataFrame(acts)
    if activity.empty:
        activity = pd.DataFrame(columns=["target_chembl_id","molecule_chembl_id","assay_chembl_id",
               "document_chembl_id","pchembl_value","standard_type","standard_value",
               "standard_units","standard_relation"])
    activity = activity.merge(mapped[["gene","target_chembl_id"]].drop_duplicates(),
                              on="target_chembl_id", how="inner")
    activity["pchembl_value"] = numeric(activity["pchembl_value"])
    activity.to_csv(out / "chembl_activities.csv", index=False)
    logging.info("Bioactivity records %d", len(activity))

    all_failed = {}
    for kind, ids, index_col in [
        ("assay", activity.assay_chembl_id.dropna().unique(), "assay_chembl_id"),
        ("molecule", activity.molecule_chembl_id.dropna().unique(), "molecule_chembl_id"),
        ("mechanism", activity.molecule_chembl_id.dropna().unique(), "molecule_chembl_id")
    ]:
        jobs, work = bulk_metadata(kind, ids, root)
        response, failures = parallel(jobs, work, args.workers, kind)
        failed_ids = set(x for job in failures for x in job[1])
        all_failed[kind] = failed_ids
        records = [r for batch in response.values() for r in batch]
        info = pd.DataFrame(records)
        if kind == "assay":
            fields = ["assay_chembl_id","confidence_score","assay_type","relationship_type"]
            info = info.reindex(columns=fields).drop_duplicates("assay_chembl_id")
            activity = activity.merge(info, on="assay_chembl_id", how="left")
            activity["confidence_score"] = numeric(activity.confidence_score)
        elif kind == "molecule":
            fields = ["molecule_chembl_id","pref_name","max_phase"]
            info = info.reindex(columns=fields).drop_duplicates("molecule_chembl_id")
            activity = activity.merge(info, on="molecule_chembl_id", how="left")
            activity["max_phase"] = numeric(activity.max_phase)
        else:
            fields = ["molecule_chembl_id","target_chembl_id","action_type","mechanism_of_action"]
            info = info.reindex(columns=fields).drop_duplicates()
            # Match BOTH compound and target. Avoid many-to-many multiplication
            info = info.groupby(["molecule_chembl_id","target_chembl_id"], as_index=False).agg(
                action_type=("action_type", lambda x: ";".join(sorted(set(x.dropna().astype(str))))),
                mechanism_of_action=("mechanism_of_action", lambda x: ";".join(sorted(set(x.dropna().astype(str)))))) if not info.empty else info
            # Убираем возможные дублирующиеся колонки
            activity = activity.drop(
                columns=["action_type", "mechanism_of_action"],
                errors="ignore"
            )

            # Добавляем аннотированные механизмы из ChEMBL
            activity = activity.merge(
                info,
                on=["molecule_chembl_id", "target_chembl_id"],
                how="left",
                validate="many_to_one"
            )
        info.to_csv(out / f"chembl_{kind}_metadata.csv", index=False)
        logging.info("%s records %d; failed IDs %d", kind, len(info), len(failed_ids))

    # Flags and summary. Missing ≠ 0 on API errors.
    activity["is_potent"] = activity.pchembl_value.ge(args.min_pchembl)
    activity["confidence_9"] = activity.confidence_score.eq(9)
    activity["confidence_7plus"] = activity.confidence_score.ge(7)
    activity["is_strict"] = activity.is_potent & activity.confidence_9
    activity["is_inhibitor"] = activity.action_type.fillna("").str.upper().str.contains(
        r"INHIBITOR|ANTAGONIST|DEGRADER", regex=True)
    activity.to_csv(out / "chembl_compound_target_pairs.csv", index=False)
    rows = []
    for gene in genes:
        gene_tids = set(mapped.loc[mapped.gene.eq(gene), "target_chembl_id"])
        sub = activity.loc[activity.gene.eq(gene)]
        mapped_ok = bool(gene_tids)
        act_ok = mapped_ok and not bool(gene_tids & set(failed_tids))
        a_ok = act_ok and not bool(set(sub.assay_chembl_id.dropna()) & all_failed["assay"])
        c_ok = act_ok and not bool(set(sub.molecule_chembl_id.dropna()) & all_failed["molecule"])
        m_ok = act_ok and not bool(set(sub.molecule_chembl_id.dropna()) & all_failed["mechanism"])
        d = dict(gene=gene, chembl_mapping_success=int(mapped_ok), chembl_activity_complete=int(act_ok),
                 chembl_assay_complete=int(a_ok), chembl_compound_complete=int(c_ok),
                 chembl_mechanism_complete=int(m_ok), chembl_n_targets=len(gene_tids) if mapped_ok else None)
        if act_ok:
            d.update(chembl_n_compounds=sub.molecule_chembl_id.nunique(),
                     chembl_n_potent_compounds=sub.loc[sub.is_potent, "molecule_chembl_id"].nunique(),
                     chembl_n_assays=sub.assay_chembl_id.nunique(),
                     chembl_n_documents=sub.document_chembl_id.nunique(),
                     chembl_best_pchembl=sub.pchembl_value.max(),
                     chembl_median_pchembl=sub.pchembl_value.median())
        if a_ok:
            d.update(chembl_n_compounds_conf9=sub.loc[sub.confidence_9,"molecule_chembl_id"].nunique(),
                     chembl_n_compounds_conf7plus=sub.loc[sub.confidence_7plus,"molecule_chembl_id"].nunique(),
                     chembl_n_strict_compounds=sub.loc[sub.is_strict,"molecule_chembl_id"].nunique(),
                     chembl_max_confidence=sub.confidence_score.max())
        if c_ok:
            d.update(chembl_max_phase=sub.max_phase.max(),
                     chembl_n_clinical_compounds=sub.loc[sub.max_phase.ge(1),"molecule_chembl_id"].nunique())
        if m_ok:
            inhib = sub.loc[sub.is_inhibitor,"molecule_chembl_id"].nunique()
            d.update(chembl_n_inhibitors=inhib, chembl_has_inhibitor=int(inhib > 0))
        rows.append(d)
    gene_features = pd.DataFrame(rows)
    gene_features.to_csv(out / "chembl_gene_features.csv", index=False)
    with open(out / "query_status.json", "w") as f:
        json.dump(dict(failed_targets=list(failed_tids),
                       failed_assay_ids=sorted(all_failed["assay"]),
                       failed_molecule_ids=sorted(all_failed["molecule"]),
                       failed_mechanism_ids=sorted(all_failed["mechanism"]),
                       mapping_failed=sorted(mapping_failed)), f, indent=2)
    logging.info("FINISHED: %d genes; results in %s", len(genes), out.resolve())

if __name__ == "__main__":
    main()