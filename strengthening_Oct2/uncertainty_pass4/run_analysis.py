"""Run specified post-discovery sensitivity scenarios, preserving source artifacts."""
import ast
import csv
import hashlib
import itertools
import json
from pathlib import Path

from uncertainty import auc, flip_bounds, zero_event_upper

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
BASE = PARENT.parent
SRC = BASE / 'ev_capture_model/src/aim13_scoping'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_tsv(path, rows):
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    specification_hash = digest(HERE / 'analysis_specification.md')
    hash_path = HERE / 'analysis_specification.sha256'
    if hash_path.exists():
        assert hash_path.read_text().strip() == specification_hash
    else:
        hash_path.write_text(specification_hash + '\n')
    module = ast.parse((SRC / 'target_labels.py').read_text(encoding='utf-8-sig'))
    targets = next(ast.literal_eval(node.value) for node in module.body if isinstance(node, ast.Assign)
                   and any(isinstance(t, ast.Name) and t.id == 'TARGETS' for t in node.targets))
    labels = {r[0]: r[3] == 'safe' for r in targets if r[3] in ['safe', 'toxic']}
    raw = json.loads((SRC / 'wholebody_v2.json').read_text())['reference_A_cell_type']
    flags = json.loads((SRC / 'flag_cleanup.json').read_text())['targets']
    genes30 = list(labels)
    assert len(genes30) == 30 and sum(labels.values()) == 19
    data = {g: {'confinement': raw[g]['top_share'], 'blood_flag': float(flags[g]['M2_flag'])} for g in genes30}
    for row in data.values():
        row['frozen_rule'] = row['confinement'] + row['blood_flag']
    names = ['frozen_rule', 'confinement', 'blood_flag']

    def evaluate(genes, y):
        return {name: auc([data[g][name] for g in genes], y) for name in names}

    genes = [g for g in genes30 if g != 'FAP']
    y = [labels[g] for g in genes]
    scenarios = []
    original = evaluate(genes30, [labels[g] for g in genes30])
    assert abs(original['frozen_rule'] - 191/209) < 1e-14
    scenarios.append({'scenario': 'original_30', 'n': 30, 'favorable': 19, 'unfavorable': 11, **original})
    for name, changed in [('exclude_FAP', []), ('exclude_FAP_change_GPRC5D', ['GPRC5D']),
                          ('exclude_FAP_change_NECTIN4', ['NECTIN4']),
                          ('exclude_FAP_change_both', ['GPRC5D', 'NECTIN4'])]:
        yy = [False if g in changed else labels[g] for g in genes]
        scenarios.append({'scenario': name, 'n': len(genes), 'favorable': sum(yy),
                          'unfavorable': len(genes)-sum(yy), **evaluate(genes, yy)})
    write_tsv(HERE / 'scenario_results.tsv', scenarios)

    uncertain = ['FOLR1', 'EGFR', 'EPCAM', 'NCAM1', 'F3', 'ERBB2', 'CD47', 'GPRC5D', 'NECTIN4', 'CLDN18']
    assert len(set(uncertain)) == 10 and set(uncertain) <= set(genes)
    assignments = []
    for bits in itertools.product([False, True], repeat=len(uncertain)):
        overrides = dict(zip(uncertain, bits))
        yy = [overrides.get(g, labels[g]) for g in genes]
        metrics = evaluate(genes, yy)
        assignments.append({'assignment_id': len(assignments), 'uncertain_assigned_favorable': ';'.join(g for g, bit in overrides.items() if bit),
                            'favorable': sum(yy), 'unfavorable': len(genes)-sum(yy), **metrics,
                            'rule_minus_confinement': metrics['frozen_rule']-metrics['confinement'],
                            'rule_minus_blood': metrics['frozen_rule']-metrics['blood_flag']})
    assert len(assignments) == 1024
    write_tsv(HERE / 'all_1024_assignments.tsv', assignments)
    metrics = names + ['rule_minus_confinement', 'rule_minus_blood']
    envelope = {name: {'min': min(assignments, key=lambda r: r[name]),
                       'max': max(assignments, key=lambda r: r[name])} for name in metrics}
    curves = []
    for name in names:
        for k in range(6):
            bound = flip_bounds([data[g][name] for g in genes], y, k)
            curves.append({'model': name, 'exact_label_changes': k,
                           'minimum_auc': bound['min']['auc'], 'maximum_auc': bound['max']['auc'],
                           'minimum_witness_flipped_genes': ';'.join(genes[i] for i in bound['min']['flipped_indices']),
                           'maximum_witness_flipped_genes': ';'.join(genes[i] for i in bound['max']['flipped_indices'])})
    write_tsv(HERE / 'arbitrary_flip_bounds.tsv', curves)
    inputs = [{'gene': g, 'historical_favorable': labels[g], **data[g],
               'excluded_human_scope': g == 'FAP', 'in_conditional_uncertainty_set': g in uncertain} for g in genes30]
    write_tsv(HERE / 'analysis_inputs.tsv', inputs)
    result = {
        'status': 'Hypothetical post-discovery sensitivity; no final amended clinical outcomes',
        'scenarios': scenarios, 'uncertainty_genes': uncertain,
        'fixed_favorable_genes': [g for g in genes if g not in uncertain and labels[g]],
        'fixed_unfavorable_genes': [g for g in genes if g not in uncertain and not labels[g]],
        'assignment_count': len(assignments), 'conditional_envelopes': envelope,
        'arbitrary_flip_bounds': curves,
        'zero_event_illustrations': [{'n': n, 'one_sided_95_upper': zero_event_upper(n)} for n in [4, 19, 27, 40, 100]],
        'limits': 'Ranges are conditional, not confidence intervals. Assignments have no probabilities. Fixed labels remain unadjudicated. No fitted model or independent validation.'
    }
    (HERE / 'results.json').write_text(json.dumps(result, indent=2) + '\n')

    manifest = json.loads((PARENT / 'before_hashes.json').read_text())
    modified = {r['file'].replace('\\', '/') for r in json.loads((PARENT / 'document_changes.json').read_text())}
    preserved = [rel for rel in manifest if rel.replace('\\', '/') not in modified]
    for rel in preserved:
        assert digest(BASE / rel) == manifest[rel], rel
    previous_count = 0
    for folder in ['clinical_review_pass2', 'clinical_review_pass3']:
        for rel, value in json.loads((PARENT / folder / 'artifact_hashes.json').read_text()).items():
            assert digest(PARENT / folder / rel) == value, rel
            previous_count += 1
    checks = {'original_auc_reproduced': original['frozen_rule'], 'historical_files_unchanged': len(preserved),
              'prior_pass_files_unchanged': previous_count, 'all_assignments_exhaustively_evaluated': len(assignments),
              'specification_sha256': specification_hash, 'source_hashes': {p.name: digest(p) for p in [SRC/'target_labels.py', SRC/'wholebody_v2.json', SRC/'flag_cleanup.json']}}
    (HERE / 'checks.json').write_text(json.dumps(checks, indent=2) + '\n')
    print(json.dumps({'scenarios': scenarios,
                      'envelopes': {name: [envelope[name]['min'][name], envelope[name]['max'][name]] for name in metrics},
                      'arbitrary_rule_bounds': [r for r in curves if r['model'] == 'frozen_rule'],
                      'zero_event_illustrations': result['zero_event_illustrations'], 'checks': checks}, indent=2))


if __name__ == '__main__':
    main()
