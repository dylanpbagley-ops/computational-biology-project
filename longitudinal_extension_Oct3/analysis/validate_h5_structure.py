from __future__ import annotations
import csv, json, sys
from collections import Counter
from pathlib import Path
import numpy as np

STAGE = Path(__file__).resolve().parents[1]
ROOT = STAGE.parent
sys.path.insert(0, str(STAGE / 'dependencies'))
import h5py

manifest = json.loads((STAGE / 'data' / 'download_manifest.json').read_text(encoding='utf-8'))
protocol = STAGE / 'analysis' / 'pilot_protocol.md'
lock = json.loads((STAGE / 'analysis' / 'pilot_protocol.lock.json').read_text(encoding='utf-8'))
import hashlib
if hashlib.sha256(protocol.read_bytes()).hexdigest() != lock['protocol_sha256']:
    raise RuntimeError('Locked pilot protocol changed after lock')
targets = ['TNFRSF17', 'GPRC5D', 'SLAMF7']
rows = []
for item in manifest['files']:
    path = STAGE / item['relative_path']
    if not path.is_file() or path.stat().st_size != item['bytes']:
        raise RuntimeError(f'File missing or size mismatch: {item["gsm"]}')
    with h5py.File(path, 'r') as h5:
        if 'matrix' not in h5:
            raise RuntimeError(f'{item["gsm"]}: no /matrix group')
        m = h5['matrix']
        required = {'data','indices','indptr','shape','barcodes','features'}
        missing = sorted(required - set(m.keys()))
        if missing:
            raise RuntimeError(f'{item["gsm"]}: missing expected 10x matrix paths: {missing}')
        f = m['features']
        missing_features = sorted({'id','name'} - set(f.keys()))
        if missing_features:
            raise RuntimeError(f'{item["gsm"]}: missing feature metadata {missing_features}')
        def strings(ds):
            vals = ds[:]
            return [x.decode('utf-8') if isinstance(x, bytes) else str(x) for x in vals]
        ids, names, barcodes = strings(f['id']), strings(f['name']), strings(m['barcodes'])
        shape = tuple(int(x) for x in m['shape'][:])
        data_ds, indices_ds, indptr_ds = m['data'], m['indices'], m['indptr']
        if len(shape) != 2 or shape != (len(ids), len(barcodes)):
            raise RuntimeError(f'{item["gsm"]}: shape {shape} incompatible with feature/barcode metadata')
        if len(names) != len(ids):
            raise RuntimeError(f'{item["gsm"]}: feature id/name length mismatch')
        if len(set(barcodes)) != len(barcodes):
            raise RuntimeError(f'{item["gsm"]}: duplicate barcodes')
        if len(set(ids)) != len(ids) or any(not x for x in ids):
            raise RuntimeError(f'{item["gsm"]}: duplicate or empty feature identifiers')
        if data_ds.dtype.kind not in 'iu':
            raise RuntimeError(f'{item["gsm"]}: noninteger stored data dtype {data_ds.dtype}')
        if indices_ds.dtype.kind not in 'iu' or indptr_ds.dtype.kind not in 'iu':
            raise RuntimeError(f'{item["gsm"]}: indices/indptr are not integer arrays')
        if data_ds.shape != indices_ds.shape or len(indptr_ds) != len(barcodes) + 1:
            raise RuntimeError(f'{item["gsm"]}: CSC array lengths are inconsistent')
        ptr = indptr_ds[:]
        idx = indices_ds[:]
        values = data_ds[:]
        if len(ptr) == 0 or int(ptr[0]) != 0 or int(ptr[-1]) != len(values) or np.any(ptr[1:] < ptr[:-1]):
            raise RuntimeError(f'{item["gsm"]}: invalid CSC indptr')
        if len(idx) and (int(idx.min()) < 0 or int(idx.max()) >= shape[0]):
            raise RuntimeError(f'{item["gsm"]}: feature row index outside matrix shape')
        if any(len(np.unique(idx[int(ptr[j]):int(ptr[j+1])])) != int(ptr[j+1] - ptr[j]) for j in range(len(barcodes))):
            raise RuntimeError(f'{item["gsm"]}: duplicate feature-row coordinates within a barcode column')
        if len(values) and int(values.min()) < 0:
            raise RuntimeError(f'{item["gsm"]}: negative count values')
        if 'feature_type' not in f:
            raise RuntimeError(f'{item["gsm"]}: feature_type missing; cannot confirm Gene Expression-only matrix from source metadata')
        types = strings(f['feature_type'])
        if len(types) != len(ids):
            raise RuntimeError(f'{item["gsm"]}: feature_type length mismatch')
        type_counts = Counter(types)
        gene_mask = np.asarray([t == 'Gene Expression' for t in types])
        target_matches = {g: [i for i, (name, ft) in enumerate(zip(names, types)) if name == g and ft == 'Gene Expression'] for g in targets}
        bad = {g: len(ix) for g, ix in target_matches.items() if len(ix) != 1}
        if bad:
            raise RuntimeError(f'{item["gsm"]}: targets not uniquely represented as Gene Expression features: {bad}')
        rows.append({
            'gsm': item['gsm'], 'path': item['relative_path'], 'sha256': item['sha256'],
            'shape_features_cells': f'{shape[0]}x{shape[1]}', 'n_features': len(ids), 'n_barcodes': len(barcodes),
            'n_stored_values': len(values), 'data_dtype': str(data_ds.dtype), 'indices_dtype': str(indices_ds.dtype),
            'unique_barcode_count': len(set(barcodes)), 'unique_feature_id_count': len(set(ids)),
            'duplicate_feature_row_coordinates': False,
            'gene_expression_features': int(gene_mask.sum()), 'feature_type_counts_json': json.dumps(type_counts, sort_keys=True),
            'target_symbol_features_each': 1, 'raw_integer_nonnegative_values': True,
            'csc_structure_valid': True, 'status': 'PASS; no target expression outcomes computed'
        })
print('\t'.join(['gsm','n_features','n_barcodes','n_stored_values','data_dtype','gene_expression_features','feature_type_counts_json','status']))
for r in rows:
    print('\t'.join(str(r[k]) for k in ['gsm','n_features','n_barcodes','n_stored_values','data_dtype','gene_expression_features','feature_type_counts_json','status']))
out = STAGE / 'data' / 'input_structure_qc.tsv'
with out.open('w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
    w.writeheader(); w.writerows(rows)
print(f'Wrote {out.relative_to(STAGE)}; validated {len(rows)} files without calculating target outcomes')
