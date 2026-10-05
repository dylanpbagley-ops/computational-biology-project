from __future__ import annotations
import csv, hashlib, json, math, statistics, sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import numpy as np

STAGE = Path(__file__).resolve().parents[1]
ROOT = STAGE.parent
sys.path.insert(0, str(STAGE / 'dependencies'))
from load_paired_counts import TARGETS, PAIR_ORDER, PAIR_NAMES, joint_thin_counts, load_sample, summarize_pair

PAIRS = {
    'PM01': {'pre': ('GSM5971019', 'P1393'), 'post': ('GSM5971023', 'P1555')},
    'PM07': {'pre': ('GSM5971021', 'P1389'), 'post': ('GSM5971025', 'P1479')},
    'PM10': {'pre': ('GSM5971020', 'P1217'), 'post': ('GSM5971024', 'P1438')},
}
SAMPLE_META = {gsm: {'patient': patient, 'visit': visit, 'specimen': specimen}
               for patient, visits in PAIRS.items() for visit, (gsm, specimen) in visits.items()}
SEEDS = list(range(20261003, 20261103))


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def csv_write(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f'Refusing empty output table: {path}')
    with path.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n', extrasaction='raise')
        writer.writeheader()
        writer.writerows(rows)


def finite_cell(value):
    if value is None:
        return 'NA'
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        value = float(value)
        if not math.isfinite(value):
            return 'NA'
    return value


def metric_value(metrics: dict, key: str):
    return metrics[key]


def main() -> int:
    started_local = datetime.now().astimezone()
    started_utc = datetime.now(timezone.utc)
    analysis_dir = STAGE / 'analysis'
    data_dir = STAGE / 'data'
    outdir = analysis_dir / 'results'
    outdir.mkdir(parents=True, exist_ok=True)
    protocol = analysis_dir / 'pilot_protocol.md'
    lock_path = analysis_dir / 'pilot_protocol.lock.json'
    lock = json.loads(lock_path.read_text(encoding='utf-8'))
    protocol_sha = digest(protocol)
    if protocol_sha != lock['protocol_sha256']:
        raise RuntimeError('Locked protocol hash differs from the saved protocol')
    required_pairs = {
        'PM01': (('GSM5971019','P1393'),('GSM5971023','P1555')),
        'PM07': (('GSM5971021','P1389'),('GSM5971025','P1479')),
        'PM10': (('GSM5971020','P1217'),('GSM5971024','P1438')),
    }
    for patient, (pre, post) in required_pairs.items():
        exact = f"- {patient}: {pre[0]} / {pre[1]} pre → {post[0]} / {post[1]} post."
        if exact not in protocol.read_text(encoding='utf-8'):
            raise RuntimeError(f'Locked specimen map mismatch for {patient}')
    manifest = json.loads((data_dir / 'download_manifest.json').read_text(encoding='utf-8'))
    expected_gsms = {gsm for v in PAIRS.values() for gsm, _ in v.values()}
    if {r['gsm'] for r in manifest['files']} != expected_gsms or manifest['total_bytes'] > 500 * 1024 * 1024:
        raise RuntimeError('Download manifest differs from locked six-sample set or exceeds total cap')
    qc_path = data_dir / 'input_structure_qc.tsv'
    with qc_path.open(encoding='utf-8-sig', newline='') as f:
        qc = {r['gsm']: r for r in csv.DictReader(f, delimiter='\t')}
    if set(qc) != expected_gsms or any(not r['status'].startswith('PASS') or r['duplicate_feature_row_coordinates'] != 'False' for r in qc.values()):
        raise RuntimeError('Input structure QC has a missing/failed sample or duplicate coordinate')
    file_rows = {r['gsm']: r for r in manifest['files']}
    loaded = {}
    for gsm in lock['gsm_accessions']:
        source = file_rows[gsm]
        path = STAGE / source['relative_path']
        if path.stat().st_size != int(source['bytes']) or digest(path) != source['sha256']:
            raise RuntimeError(f'Downloaded input hash/size mismatch: {gsm}')
        loaded[gsm] = load_sample(path, gsm)
        if len(loaded[gsm]['all_feature_types']) != 1 or loaded[gsm]['all_feature_types'][0] != 'Gene Expression':
            raise RuntimeError(f'{gsm}: unexpected modality appeared after structure QC')
        if loaded[gsm]['n_gene_expression_features'] != 33538:
            raise RuntimeError(f'{gsm}: changed feature count after structure QC')

    # Primary specimen-level values and a specimen-qualified raw source table.
    cell_rows = []
    sample_depth = {}
    primary_rows = []
    primary_lookup = {}
    for gsm in lock['gsm_accessions']:
        meta = SAMPLE_META[gsm]
        item = loaded[gsm]
        totals = item['total_umi']
        counts = item['target_counts']
        if not len(totals) or len(counts) != len(totals):
            raise RuntimeError(f'{gsm}: no complete barcode/count data')
        quantiles = np.quantile(totals, [0, .25, .5, .75, 1], method='linear')
        sample_depth[gsm] = {
            'n_cells': len(totals), 'minimum': int(quantiles[0]), 'p25': float(quantiles[1]),
            'median': float(quantiles[2]), 'p75': float(quantiles[3]), 'maximum': int(quantiles[4]),
            'mean': float(np.mean(totals)),
        }
        for i, barcode in enumerate(item['barcodes']):
            cell_rows.append({
                'cell_key': item['cell_keys'][i], 'gsm': gsm, 'participant': meta['patient'],
                'venetoclax_visit': meta['visit'], 'specimen_key': meta['specimen'], 'barcode': barcode,
                'total_gene_expression_umi': int(totals[i]),
                'TNFRSF17_umi': int(counts[i,0]), 'GPRC5D_umi': int(counts[i,1]), 'SLAMF7_umi': int(counts[i,2]),
                'cell_scope': 'CD138-selected filtered-matrix barcode; malignant state not assigned'
            })
        for pair_index, (ai, bi) in enumerate(PAIR_ORDER):
            metrics = summarize_pair(counts, ai, bi)
            row = {
                'participant': meta['patient'], 'venetoclax_visit': meta['visit'], 'specimen_key': meta['specimen'],
                'gsm': gsm, 'marker_A': PAIR_NAMES[pair_index][0], 'marker_B': PAIR_NAMES[pair_index][1],
                'n_cells': len(totals), 'raw_count_detection_rule': 'integer UMI count > 0',
            }
            row.update(metrics)
            primary_rows.append({k: finite_cell(v) for k, v in row.items()})
            primary_lookup[(gsm, pair_index)] = metrics
    csv_write(outdir / 'per_barcode_source_counts.tsv', cell_rows)
    csv_write(outdir / 'sample_depth_summary.tsv', [
        {'gsm': gsm, 'participant': SAMPLE_META[gsm]['patient'], 'visit': SAMPLE_META[gsm]['visit'],
         'specimen_key': SAMPLE_META[gsm]['specimen'], **sample_depth[gsm],
         'feature_scope': 'all 33,538 Gene Expression features', 'barcode_scope': 'all filtered matrix barcodes'}
        for gsm in lock['gsm_accessions']
    ])
    csv_write(outdir / 'primary_joint_detection.tsv', primary_rows)

    # Paired specimen changes are differences in proportions/conditionals, never raw-cell count subtraction.
    delta_rows = []
    figure_rows = []
    delta_metrics = [
        ('n11_fraction','both detected fraction'), ('n10_fraction','A only fraction'),
        ('n01_fraction','B only fraction'), ('n00_fraction','neither detected fraction'),
        ('a_detected_fraction','A marginal detection fraction'), ('b_detected_fraction','B marginal detection fraction'),
        ('or_coverage_fraction','OR coverage fraction'), ('double_undetected_fraction','double-undetected fraction'),
        ('b_given_a_undetected','B detection among A-undetected'), ('a_given_b_undetected','A detection among B-undetected'),
    ]
    for patient, visits in PAIRS.items():
        pre_gsm, pre_specimen = visits['pre']
        post_gsm, post_specimen = visits['post']
        for pair_index, pair in enumerate(PAIR_NAMES):
            pre = primary_lookup[(pre_gsm, pair_index)]
            post = primary_lookup[(post_gsm, pair_index)]
            for key, label in delta_metrics:
                pv, qv = pre[key], post[key]
                delta = (qv - pv) if pv is not None and qv is not None else None
                delta_rows.append({
                    'participant': patient, 'pre_gsm': pre_gsm, 'pre_specimen': pre_specimen,
                    'post_gsm': post_gsm, 'post_specimen': post_specimen,
                    'marker_A': pair[0], 'marker_B': pair[1], 'metric': key, 'metric_label': label,
                    'pre_value': finite_cell(pv), 'post_value': finite_cell(qv), 'post_minus_pre': finite_cell(delta),
                    'pre_n_cells': sample_depth[pre_gsm]['n_cells'], 'post_n_cells': sample_depth[post_gsm]['n_cells'],
                    'unit': 'fraction; post minus pre; descriptive participant-paired change'
                })
            for category_key, category in [('n11_fraction','both'), ('n10_fraction','A_only'), ('n01_fraction','B_only'), ('n00_fraction','neither')]:
                figure_rows.append({
                    'participant': patient, 'marker_A': pair[0], 'marker_B': pair[1], 'category': category,
                    'post_minus_pre_fraction': finite_cell(post[category_key] - pre[category_key]),
                    'pre_gsm': pre_gsm, 'post_gsm': post_gsm,
                })
    csv_write(outdir / 'participant_paired_deltas.tsv', delta_rows)
    csv_write(outdir / 'figure_joint_category_deltas.tsv', figure_rows)

    specimen_medians = [sample_depth[gsm]['median'] for gsm in lock['gsm_accessions']]
    depth_D = math.floor(min(specimen_medians))
    threshold_rows = []
    for gsm in lock['gsm_accessions']:
        below = int(np.count_nonzero(loaded[gsm]['total_umi'] < depth_D))
        n_cells = sample_depth[gsm]['n_cells']
        threshold_rows.append({
            'gsm': gsm, 'participant': SAMPLE_META[gsm]['patient'], 'visit': SAMPLE_META[gsm]['visit'],
            'specimen_key': SAMPLE_META[gsm]['specimen'], 'D': depth_D,
            'n_cells': n_cells, 'n_cells_below_D': below, 'fraction_below_D': below/n_cells,
            'median_total_gene_expression_umi': sample_depth[gsm]['median']
        })
    csv_write(outdir / 'depth_threshold_by_sample.tsv', threshold_rows)
    rng_probe = np.random.default_rng(20261003).multivariate_hypergeometric([3,4,2,30], 5, method='marginals')
    if int(rng_probe.sum()) != 5 or np.any(rng_probe > np.array([3,4,2,30])):
        raise RuntimeError('NumPy multivariate hypergeometric smoke check failed')
    depth_rows, paired_sens_rows = [], []
    if depth_D < 1:
        depth_status = 'UNAVAILABLE: D < 1; unresolved depth confounding'
    else:
        depth_status = 'PASS: joint multivariate hypergeometric sensitivity'
        for seed in SEEDS:
            child_seeds = np.random.SeedSequence(seed).generate_state(len(lock['gsm_accessions']), dtype=np.uint32)
            sens_lookup = {}
            for sample_pos, gsm in enumerate(lock['gsm_accessions']):
                item = loaded[gsm]
                counts = joint_thin_counts(item['target_counts'], item['total_umi'], depth_D, int(child_seeds[sample_pos]))
                below = int(np.count_nonzero(item['total_umi'] < depth_D))
                for pair_index, (ai, bi) in enumerate(PAIR_ORDER):
                    metrics = summarize_pair(counts, ai, bi)
                    sens_lookup[(gsm, pair_index)] = metrics
                    depth_rows.append({
                        'replicate_seed': seed, 'sample_rng_seed': int(child_seeds[sample_pos]),
                        'participant': SAMPLE_META[gsm]['patient'], 'venetoclax_visit': SAMPLE_META[gsm]['visit'],
                        'gsm': gsm, 'marker_A': PAIR_NAMES[pair_index][0], 'marker_B': PAIR_NAMES[pair_index][1],
                        'D': depth_D, 'n_cells': metrics['n_cells'], 'n_cells_below_D': below,
                        'fraction_below_D': below / metrics['n_cells'],
                        **{k: finite_cell(v) for k, v in metrics.items() if k != 'n_cells'}
                    })
            for patient, visits in PAIRS.items():
                pre_gsm, _ = visits['pre']; post_gsm, _ = visits['post']
                for pair_index, pair in enumerate(PAIR_NAMES):
                    pre = sens_lookup[(pre_gsm,pair_index)]; post = sens_lookup[(post_gsm,pair_index)]
                    for key, label in delta_metrics:
                        pv, qv = pre[key], post[key]
                        paired_sens_rows.append({
                            'replicate_seed': seed, 'participant': patient,
                            'marker_A': pair[0], 'marker_B': pair[1], 'metric': key,
                            'pre_gsm': pre_gsm, 'post_gsm': post_gsm,
                            'post_minus_pre': finite_cell(qv-pv if pv is not None and qv is not None else None),
                            'D': depth_D
                        })
        csv_write(outdir / 'depth_sensitivity_replicates.tsv', depth_rows)
        csv_write(outdir / 'depth_sensitivity_paired_deltas.tsv', paired_sens_rows)
        summary_rows = []
        group = defaultdict(list)
        for r in depth_rows:
            for key, _ in delta_metrics:
                if r[key] != 'NA':
                    group[(r['gsm'], r['marker_A'], r['marker_B'], key)].append(float(r[key]))
        for (gsm, marker_a, marker_b, key), vals in group.items():
            q = np.quantile(vals, [.05,.5,.95])
            summary_rows.append({
                'gsm': gsm, 'participant': SAMPLE_META[gsm]['patient'], 'visit': SAMPLE_META[gsm]['visit'],
                'marker_A': marker_a, 'marker_B': marker_b, 'metric': key, 'replicates': len(vals),
                'mean': float(np.mean(vals)), 'p05': float(q[0]), 'median': float(q[1]), 'p95': float(q[2]),
                'minimum': float(min(vals)), 'maximum': float(max(vals)), 'D': depth_D,
                'interpretation': 'within-specimen sampling sensitivity; not additional participants'
            })
        csv_write(outdir / 'depth_sensitivity_summary.tsv', summary_rows)
        paired_summary = []
        pgroup = defaultdict(list)
        for r in paired_sens_rows:
            if r['post_minus_pre'] != 'NA':
                pgroup[(r['participant'],r['marker_A'],r['marker_B'],r['metric'])].append(float(r['post_minus_pre']))
        for (patient, marker_a, marker_b, key), vals in pgroup.items():
            primary = next(r for r in delta_rows if r['participant']==patient and r['marker_A']==marker_a and r['marker_B']==marker_b and r['metric']==key)
            primary_delta = None if primary['post_minus_pre']=='NA' else float(primary['post_minus_pre'])
            same_sign = None if primary_delta is None or primary_delta == 0 else sum((v > 0) == (primary_delta > 0) and v != 0 for v in vals) / len(vals)
            q = np.quantile(vals,[.05,.5,.95])
            paired_summary.append({
                'participant': patient, 'marker_A': marker_a, 'marker_B': marker_b, 'metric': key,
                'replicates': len(vals), 'primary_post_minus_pre': finite_cell(primary_delta),
                'depth_mean_post_minus_pre': float(np.mean(vals)), 'depth_p05': float(q[0]),
                'depth_median': float(q[1]), 'depth_p95': float(q[2]), 'depth_min': float(min(vals)),
                'depth_max': float(max(vals)), 'fraction_replicates_same_nonzero_sign_as_primary': finite_cell(same_sign),
                'D': depth_D, 'interpretation': 'descriptive within-specimen joint-thinning sensitivity'
            })
        csv_write(outdir / 'depth_sensitivity_paired_summary.tsv', paired_summary)
    if depth_D >= 1 and not depth_rows:
        raise RuntimeError('Depth sensitivity unexpectedly produced no replicate rows')

    cell_scope = 'CD138-selected filtered-matrix barcodes; no per-cell malignant labels; specimens are not lineage-matched'
    result_manifest = {
        'analysis': 'fixed three-marker paired descriptive transcript-detection pilot',
        'started_at_local': started_local.isoformat(timespec='seconds'),
        'started_at_utc': started_utc.isoformat(timespec='seconds').replace('+00:00','Z'),
        'protocol_sha256': protocol_sha, 'input_download_manifest_sha256': digest(data_dir / 'download_manifest.json'),
        'n_participants': len(PAIRS), 'n_specimens': len(SAMPLE_META),
        'participant_specimen_pairs': {p:{v:{'gsm':x[0],'specimen':x[1]} for v,x in visits.items()} for p,visits in PAIRS.items()},
        'cell_scope': cell_scope, 'markers': list(TARGETS), 'marker_pairs': [list(x) for x in PAIR_NAMES],
        'primary_rule': 'deposited nonnegative integer RNA UMI count > 0',
        'depth_sensitivity': {'status': depth_status, 'D': depth_D, 'D_rule': 'floor(minimum of six specimen medians over Gene Expression UMIs)',
                              'seeds': SEEDS if depth_rows else [], 'replicates_per_specimen': len(SEEDS) if depth_rows else 0,
                              'sampling': 'multivariate hypergeometric over three target counts and all other Gene Expression counts combined'},
        'limitations': [
          'Only three participants; descriptive feasibility, no population-level inference or validation.',
          'No per-barcode malignant/clone annotation; do not describe barcodes as confirmed tumor cells.',
          'No protein measurement; RNA detection is not surface-antigen availability or functional escape.',
          'Separate biopsies are not lineage tracked; no barcodes are joined across specimens.',
          'Pre/post labels are relative to venetoclax; no causal treatment effect is identified.',
          'Chemistry remains unresolved: six sample records say CG00052 v2 while series overview says v3.',
          'Depth sensitivity is a sampling check, not a causal correction; cells remain subsamples.'
        ],
        'outcomes': {'samples': len(loaded), 'sample_pair_rows': len(primary_rows),
                     'barcode_rows': len(cell_rows), 'depth_replicate_rows': len(depth_rows)},
        'ended_at_local': datetime.now().astimezone().isoformat(timespec='seconds'),
        'ended_at_utc': datetime.now(timezone.utc).isoformat(timespec='seconds').replace('+00:00','Z')
    }
    (outdir / 'pilot_run_manifest.json').write_text(json.dumps(result_manifest, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(f'Pilot complete: {len(cell_rows)} barcode rows; {len(primary_rows)} specimen/pair summaries; D={depth_D}; {len(depth_rows)} depth replicate rows; {depth_status}')
    for patient, visits in PAIRS.items():
        print(patient, visits['pre'][0], '->', visits['post'][0])
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
