from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

STAGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STAGE / 'dependencies'))
import h5py

TARGETS = ('TNFRSF17', 'GPRC5D', 'SLAMF7')
PAIR_ORDER = ((0, 1), (0, 2), (1, 2))
PAIR_NAMES = (('TNFRSF17', 'GPRC5D'), ('TNFRSF17', 'SLAMF7'), ('GPRC5D', 'SLAMF7'))


def decode_strings(dataset) -> list[str]:
    vals = dataset[:]
    return [x.decode('utf-8') if isinstance(x, (bytes, np.bytes_)) else str(x) for x in vals]


def resolve_target_features(names, feature_types, targets=TARGETS) -> list[int]:
    if len(names) != len(feature_types):
        raise ValueError('feature name/type length mismatch')
    found = []
    for target in targets:
        matches = [i for i, (name, ftype) in enumerate(zip(names, feature_types))
                   if name == target and ftype == 'Gene Expression']
        if len(matches) != 1:
            raise ValueError(f'{target}: expected one Gene Expression feature, found {len(matches)}')
        found.append(matches[0])
    return found


def validate_csc_structure(indptr, indices, n_features, n_cells) -> None:
    ptr = np.asarray(indptr)
    idx = np.asarray(indices)
    if ptr.dtype.kind not in 'iu' or idx.dtype.kind not in 'iu':
        raise ValueError('CSC pointers and indices must be integers')
    if ptr.ndim != 1 or idx.ndim != 1 or len(ptr) != n_cells + 1:
        raise ValueError('CSC pointer/index dimensions are inconsistent')
    if len(ptr) == 0 or int(ptr[0]) != 0 or int(ptr[-1]) != len(idx) or np.any(ptr[1:] < ptr[:-1]):
        raise ValueError('invalid CSC indptr')
    if len(idx) and (int(idx.min()) < 0 or int(idx.max()) >= n_features):
        raise ValueError('CSC row index outside feature range')
    for col in range(n_cells):
        lo, hi = int(ptr[col]), int(ptr[col + 1])
        segment = idx[lo:hi]
        if len(np.unique(segment)) != len(segment):
            raise ValueError(f'duplicate feature-row coordinate in barcode column {col}')


def load_sample(path: Path, gsm: str) -> dict:
    with h5py.File(path, 'r') as h5:
        if 'matrix' not in h5:
            raise ValueError(f'{gsm}: missing /matrix group')
        matrix = h5['matrix']
        required = {'data', 'indices', 'indptr', 'shape', 'barcodes', 'features'}
        if not required.issubset(matrix.keys()):
            raise ValueError(f'{gsm}: missing 10x matrix datasets {sorted(required - set(matrix.keys()))}')
        feature_group = matrix['features']
        if not {'id', 'name', 'feature_type'}.issubset(feature_group.keys()):
            raise ValueError(f'{gsm}: missing feature id, name, or type metadata')
        ids = decode_strings(feature_group['id'])
        names = decode_strings(feature_group['name'])
        types = decode_strings(feature_group['feature_type'])
        barcodes = decode_strings(matrix['barcodes'])
        if len(set(barcodes)) != len(barcodes) or any(not b for b in barcodes):
            raise ValueError(f'{gsm}: duplicate or empty barcodes')
        if len(set(ids)) != len(ids) or any(not x for x in ids):
            raise ValueError(f'{gsm}: duplicate or empty feature IDs')
        shape = tuple(int(x) for x in matrix['shape'][:])
        if shape != (len(ids), len(barcodes)) or len(names) != len(ids) or len(types) != len(ids):
            raise ValueError(f'{gsm}: metadata/shape dimensions are inconsistent')
        if 'Gene Expression' not in types:
            raise ValueError(f'{gsm}: no Gene Expression feature type')
        feature_types = np.asarray(types, dtype=object)
        gene_mask = feature_types == 'Gene Expression'
        target_rows = resolve_target_features(names, types)
        values = matrix['data'][:]
        indices = matrix['indices'][:]
        indptr = matrix['indptr'][:]
        if values.dtype.kind not in 'iu' or (len(values) and int(values.min()) < 0):
            raise ValueError(f'{gsm}: counts are not nonnegative integer values')
        if values.shape != indices.shape:
            raise ValueError(f'{gsm}: CSC data/index lengths differ')
        validate_csc_structure(indptr, indices, shape[0], shape[1])
        totals = np.zeros(shape[1], dtype=np.int64)
        target_counts = np.zeros((shape[1], len(TARGETS)), dtype=np.int64)
        for col in range(shape[1]):
            lo, hi = int(indptr[col]), int(indptr[col + 1])
            rows = indices[lo:hi]
            vals = values[lo:hi].astype(np.int64, copy=False)
            in_gene_expression = gene_mask[rows]
            totals[col] = vals[in_gene_expression].sum(dtype=np.int64)
            for k, target_row in enumerate(target_rows):
                at = rows == target_row
                if int(at.sum()) > 1:
                    raise ValueError(f'{gsm}: duplicate {TARGETS[k]} coordinate in a cell')
                if at.any():
                    target_counts[col, k] = vals[at].sum(dtype=np.int64)
        if np.any(target_counts.sum(axis=1) > totals):
            raise ValueError(f'{gsm}: target counts exceed Gene Expression totals')
        return {
            'gsm': gsm,
            'barcodes': barcodes,
            'cell_keys': [f'{gsm}|{barcode}' for barcode in barcodes],
            'target_counts': target_counts,
            'total_umi': totals,
            'n_gene_expression_features': int(gene_mask.sum()),
            'target_feature_rows': target_rows,
            'all_feature_types': sorted(set(types)),
        }


def summarize_pair(counts: np.ndarray, i: int, j: int, n_cells: int | None = None) -> dict:
    counts = np.asarray(counts)
    if counts.ndim != 2 or max(i, j) >= counts.shape[1]:
        raise ValueError('target count array/pair indices are invalid')
    if counts.dtype.kind not in 'iu' or (counts.size and int(counts.min()) < 0):
        raise ValueError('target counts must be nonnegative integers')
    if n_cells is None:
        n_cells = counts.shape[0]
    if n_cells != counts.shape[0] or n_cells <= 0:
        raise ValueError('invalid cell denominator')
    a = counts[:, i] > 0
    b = counts[:, j] > 0
    n11 = int(np.count_nonzero(a & b))
    n10 = int(np.count_nonzero(a & ~b))
    n01 = int(np.count_nonzero(~a & b))
    n00 = int(np.count_nonzero(~a & ~b))
    if n11 + n10 + n01 + n00 != n_cells:
        raise AssertionError('joint category counts do not sum to cell count')
    return {
        'n_cells': n_cells, 'n11_both': n11, 'n10_a_only': n10,
        'n01_b_only': n01, 'n00_neither': n00,
        'n11_fraction': n11 / n_cells, 'n10_fraction': n10 / n_cells,
        'n01_fraction': n01 / n_cells, 'n00_fraction': n00 / n_cells,
        'a_detected_n': n11 + n10, 'a_detected_fraction': (n11 + n10) / n_cells,
        'b_detected_n': n11 + n01, 'b_detected_fraction': (n11 + n01) / n_cells,
        'or_coverage_fraction': (n11 + n10 + n01) / n_cells,
        'double_undetected_fraction': n00 / n_cells,
        'b_given_a_undetected': (n01 / (n01 + n00)) if n01 + n00 else None,
        'a_given_b_undetected': (n10 / (n10 + n00)) if n10 + n00 else None,
    }


def joint_thin_counts(target_counts: np.ndarray, total_umi: np.ndarray, depth: int, seed: int) -> np.ndarray:
    counts = np.asarray(target_counts)
    totals = np.asarray(total_umi)
    if depth < 0 or int(depth) != depth:
        raise ValueError('depth must be a nonnegative integer')
    if counts.ndim != 2 or totals.ndim != 1 or len(totals) != len(counts):
        raise ValueError('target and total count dimensions differ')
    if counts.dtype.kind not in 'iu' or totals.dtype.kind not in 'iu':
        raise ValueError('counts and totals must be integer arrays')
    if (counts.size and int(counts.min()) < 0) or (totals.size and int(totals.min()) < 0):
        raise ValueError('counts and totals must be nonnegative')
    other = totals.astype(np.int64) - counts.astype(np.int64).sum(axis=1)
    if np.any(other < 0):
        raise ValueError('target counts exceed total RNA counts')
    rng = np.random.default_rng(int(seed))
    thinned = np.zeros_like(counts, dtype=np.int64)
    for row in range(len(counts)):
        n_total = int(totals[row])
        n_take = min(int(depth), n_total)
        if n_take == n_total:
            thinned[row] = counts[row]
        elif n_take:
            colors = [int(v) for v in counts[row]] + [int(other[row])]
            draw = rng.multivariate_hypergeometric(colors, n_take, method='marginals')
            if int(draw.sum()) != n_take or np.any(draw > colors):
                raise AssertionError('invalid multivariate hypergeometric draw')
            thinned[row] = draw[:-1]
    return thinned


def assert_patient_disjoint(train_patients, test_patients) -> None:
    overlap = set(train_patients) & set(test_patients)
    if overlap:
        raise ValueError(f'patient leakage across split: {sorted(overlap)}')
