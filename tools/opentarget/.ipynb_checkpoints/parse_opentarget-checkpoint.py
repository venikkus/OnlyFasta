#!/usr/bin/env python3
"""Resumable, rate-limited Open Targets target annotation for DEG tables."""
import argparse
import json
import logging
import os
import random
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from tqdm.auto import tqdm

API = 'https://api.platform.opentargets.org/api/v4/graphql'
QUERY = '''
query TargetInfo($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id approvedSymbol biotype
    tractability { label modality value }
    geneticConstraint { constraintType score oe oeUpper }
  }
}
'''
THREAD = threading.local()
RETRY_CODES = {408, 425, 429, 500, 502, 503, 504}


def get_session():
    if not hasattr(THREAD, 'session'):
        s = requests.Session()
        s.headers.update({'User-Agent': 'OpenTargets-DEG-annotation/1.0', 'Accept': 'application/json'})
        adapter = HTTPAdapter(pool_connections=2, pool_maxsize=2, max_retries=0)
        s.mount('https://', adapter)
        THREAD.session = s
    return THREAD.session


def atomic_json(path, obj):
    temp = path.with_suffix(f'.{threading.get_ident()}.tmp')
    temp.write_text(json.dumps(obj, ensure_ascii=False), encoding='utf-8')
    os.replace(temp, path)


def mapping_for_genes(genes, mapping_file):
    existing = pd.DataFrame(columns=['gene', 'ensembl_id', 'mapping_status', 'candidate_ids'])
    if mapping_file.exists():
        existing = pd.read_csv(mapping_file, dtype=str).drop_duplicates('gene', keep='last')
        if not {'gene', 'ensembl_id', 'mapping_status', 'candidate_ids'}.issubset(existing.columns):
            raise ValueError('Existing gene_mapping.csv has incompatible columns')
    missing = sorted(set(genes) - set(existing['gene'].dropna()))
    if missing:
        new_records = []
        batch_size = 50
        for i in tqdm(range(0, len(missing), batch_size), desc='MyGene mapping batches'):
            batch = missing[i:i + batch_size]
            candidates = {g: set() for g in batch}
            for attempt in range(5):
                try:
                    response = requests.post(
                        'https://mygene.info/v3/query',
                        data={
                            'q': ','.join(batch),
                            'scopes': 'symbol',
                            'fields': 'symbol,taxid,ensembl.gene',
                            'species': 'human',
                        },
                        timeout=(10, 25),
                    )
                    response.raise_for_status()
                    hits = response.json()
                    if not isinstance(hits, list):
                        raise ValueError('Unexpected MyGene response: expected list')
                    break
                except (requests.RequestException, ValueError) as exc:
                    logging.warning('MyGene batch %d (%d genes), attempt %d/5: %s',
                                    i // batch_size + 1, len(batch), attempt + 1, str(exc)[:200])
                    if attempt == 4:
                        raise RuntimeError(f'MyGene mapping failed for batch starting at {i}; earlier batches saved') from exc
                    time.sleep(min(30, 2 ** attempt + random.random()))
            for h in hits:
                requested = str(h.get('query', ''))
                if requested not in candidates or str(h.get('symbol', '')).upper() != requested.upper():
                    continue
                if h.get('taxid') != 9606:
                    continue
                ens = h.get('ensembl') or []
                if isinstance(ens, dict):
                    ens = [ens]
                for entry in ens:
                    eid = entry.get('gene') if isinstance(entry, dict) else None
                    if eid and re.fullmatch(r'ENSG\d+', str(eid)):
                        candidates[requested].add(str(eid))
            for gene in batch:
                ids = sorted(candidates[gene])
                new_records.append({
                    'gene': gene,
                    'ensembl_id': ids[0] if len(ids) == 1 else None,
                    'mapping_status': 'mapped' if len(ids) == 1 else ('ambiguous' if ids else 'not_mapped'),
                    'candidate_ids': ';'.join(ids),
                })
            # Incremental mapping checkpoint per batch
            existing = pd.concat([existing, pd.DataFrame(new_records[-len(batch):])], ignore_index=True)
            existing.to_csv(mapping_file, index=False)
    result = pd.DataFrame({'gene': genes}).merge(existing, how='left', on='gene', validate='one_to_one')
    return result


def fetch_one(ensembl_id, cache_dir, retries, timeout):
    cache_path = cache_dir / f'{ensembl_id}.json'
    if cache_path.exists():
        try:
            payload = json.loads(cache_path.read_text(encoding='utf-8'))
            if isinstance(payload, dict) and payload.get('status') in ('ok', 'not_found'):
                return ensembl_id, payload['status'], payload.get('target'), None
            # Backward compatibility with the old script's cached target JSON.
            if isinstance(payload, dict) and 'id' in payload:
                return ensembl_id, 'ok', payload, None
            if payload is None:
                return ensembl_id, 'not_found', None, None
        except (ValueError, OSError):
            pass
    last_error = ''
    for attempt in range(retries):
        try:
            response = get_session().post(API, json={
                'query': QUERY, 'variables': {'ensemblId': ensembl_id}
            }, timeout=(10, timeout))
            if response.status_code in RETRY_CODES:
                retry_after = response.headers.get('Retry-After')
                if retry_after and retry_after.isdigit():
                    delay = min(120, int(retry_after))
                else:
                    delay = min(60, 1.5 * (2 ** attempt) + random.uniform(0, 1))
                last_error = f'HTTP {response.status_code}'
                time.sleep(delay)
                continue
            response.raise_for_status()
            payload = response.json()
            if payload.get('errors'):
                raise RuntimeError(str(payload['errors'])[:400])
            if 'data' not in payload or 'target' not in payload['data']:
                raise RuntimeError('Missing data.target in API response')
            target = payload['data']['target']
            status = 'not_found' if target is None else 'ok'
            atomic_json(cache_path, {'status': status, 'target': target})
            return ensembl_id, status, target, None
        except (requests.RequestException, ValueError, RuntimeError) as exc:
            last_error = f'{type(exc).__name__}: {str(exc)[:250]}'
            if attempt < retries - 1:
                time.sleep(min(60, 1.5 * 2 ** attempt + random.uniform(0, 1)))
    return ensembl_id, 'error', None, last_error


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', default='gmp_filt.csv')
    ap.add_argument('--outdir', default='opentargets_results')
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--retries', type=int, default=6)
    ap.add_argument('--timeout', type=float, default=35)
    args = ap.parse_args()
    if args.workers < 1 or args.retries < 1:
        ap.error('--workers and --retries must be >= 1')
    out = Path(args.outdir)
    cache = out / 'cache'
    cache.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')

    deg = pd.read_csv(args.input)
    deg.columns=['gene', 'baseMean', 'log2FC', 'lfcSE', 'stat', 'pvalue', 'padj']

    if 'gene' not in deg.columns:
        raise ValueError('Input CSV must contain a gene column')
    genes = sorted(set(deg['gene'].dropna().astype(str).str.strip()) - {''})
    logging.info('Input unique genes: %d', len(genes))
    mapping = mapping_for_genes(genes, out / 'gene_mapping.csv')
    logging.info('Mapping status:\n%s', mapping['mapping_status'].value_counts(dropna=False))

    ids = sorted(mapping['ensembl_id'].dropna().astype(str).unique())
    results = {}
    # No shared Session across threads. Each request retries independently.
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(fetch_one, eid, cache, args.retries, args.timeout): eid for eid in ids}
        for f in tqdm(as_completed(futures), total=len(futures), desc='Open Targets'):
            eid = futures[f]
            try:
                k, status, target, err = f.result()
            except Exception as exc:
                k, status, target, err = eid, 'error', None, str(exc)
            results[k] = (status, target, err)

    records, raw_tract = [], []
    for _, m in mapping.iterrows():
        gene, eid = m['gene'], m['ensembl_id']
        record = {
            'gene': gene, 'ensembl_id': eid, 'mapping_status': m['mapping_status'],
            'candidate_ids': m['candidate_ids'], 'ot_request_status': 'not_mapped',
            'ot_target_found': pd.NA, 'ot_error': None,
            'ot_approved_symbol': None, 'ot_biotype': None,
            'ot_lof_oe': None, 'ot_lof_oe_upper': None, 'ot_lof_constraint_score': None,
        }
        if pd.notna(eid):
            status, target, err = results[eid]
            record['ot_request_status'] = status
            record['ot_error'] = err
            if status == 'not_found':
                record['ot_target_found'] = False
            elif status == 'ok':
                record['ot_target_found'] = True
                record['ot_approved_symbol'] = target.get('approvedSymbol')
                record['ot_biotype'] = target.get('biotype')
                lof = next((c for c in (target.get('geneticConstraint') or [])
                            if str(c.get('constraintType', '')).lower() in ('lof', 'lof_hc')), None)
                if lof:
                    record['ot_lof_oe'] = lof.get('oe')
                    record['ot_lof_oe_upper'] = lof.get('oeUpper')
                    record['ot_lof_constraint_score'] = lof.get('score')
                for t in (target.get('tractability') or []):
                    raw_tract.append({
                        'gene': gene, 'ensembl_id': eid,
                        'modality': t.get('modality'), 'label': t.get('label'), 'value': t.get('value')
                    })
        records.append(record)

    features = pd.DataFrame(records)
    tractability = pd.DataFrame(raw_tract, columns=['gene', 'ensembl_id', 'modality', 'label', 'value'])
    tractability.to_csv(out / 'opentargets_tractability_long.csv', index=False)
    if not tractability.empty:
        t = tractability.dropna(subset=['modality', 'label']).copy()
        t['bucket'] = t['modality'].astype(str) + '__' + t['label'].astype(str)
        wide = t.pivot_table(index='gene', columns='bucket', values='value', aggfunc='max')
        wide.columns = ['ot_' + re.sub(r'[^a-z0-9]+', '_', str(c).lower()).strip('_') for c in wide.columns]
        wide = wide.reset_index()
        features = features.merge(wide, how='left', on='gene', validate='one_to_one')
    features.to_csv(out / 'opentargets_features.csv', index=False)
    errors = features.loc[features.ot_request_status.eq('error'), ['gene', 'ensembl_id', 'ot_error']]
    errors.to_csv(out / 'opentargets_errors.csv', index=False)
    logging.info('Complete. found=%d, not_found=%d, errors=%d, unmapped=%d',
                 features.ot_request_status.eq('ok').sum(),
                 features.ot_request_status.eq('not_found').sum(), len(errors),
                 features.ot_request_status.eq('not_mapped').sum())
    if len(errors):
        logging.warning('Rerun the same command to retry failed IDs (successful IDs are cached).')
        raise SystemExit(2)

if __name__ == '__main__':
    main()