"""Assemble a score-withheld review template and verify evidence-table integrity."""
import ast
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
BASE = PARENT.parent

def read_tsv(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    assert all(None not in r and None not in r.values() for r in rows), path
    return rows

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

tree = ast.parse((BASE/'ev_capture_model/src/aim13_scoping/target_labels.py').read_text(encoding='utf-8'))
targets = next(ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign)
               and any(isinstance(t, ast.Name) and t.id == 'TARGETS' for t in n.targets))
favorable = {r[0] for r in targets if r[3] == 'safe'}
panel = {r[0] for r in targets if r[3] in ('safe', 'toxic')}
triage = read_tsv(HERE/'favorable_case_triage.tsv')
assert len(triage) == 19 and {r['gene'] for r in triage} == favorable
notes = read_tsv(HERE/'source_notes.tsv')
assert len(notes) == 6 and len({r['gene'] for r in notes}) == 6
assert {r['gene'] for r in notes} <= favorable
assert all(r['source'].startswith('https://') for r in notes)

# Copy source identifiers only. No prior labels, scores or assistant judgments.
old = read_tsv(PARENT/'clinical_review_pass3/source_only_review_form.tsv')
assert len(old) == 30 and {r['gene'] for r in old} == panel
identity = ['gene', 'agent_or_agents', 'primary_sources']
fields = identity + ['agent_cohort_event_id', 'human_eligibility', 'antigen_isoform',
    'endpoint_name', 'event_count', 'endpoint_evaluable_n', 'denominator_definition',
    'grade', 'injury_related_cessation', 'death', 'protocol_DLT',
    'target_attribution', 'attribution_evidence_location', 'normal_cell_depletion',
    'clinical_consequences', 'supportive_care', 'recovery',
    'exposure_activity_evidence', 'ascertainment', 'observation_window',
    'censoring_missingness', 'combination_confounders', 'overlapping_cohorts',
    'observational_state', 'adequacy_limitations', 'reviewer_name_date']
with (HERE/'review_form_score_withheld.tsv').open('w', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
    writer.writeheader()
    for row in sorted(old, key=lambda r: r['gene']):
        clean = {k: row[k] for k in identity}
        for note in notes:
            if note['gene'] == row['gene'] and note['source'] not in clean['primary_sources']:
                clean['primary_sources'] += '; ' + note['source']
        writer.writerow(clean)
form = read_tsv(HERE/'review_form_score_withheld.tsv')
assert all(not r[k] for r in form for k in fields if k not in identity)
assert not any('score' in k or 'historical' in k for k in fields)

manifest = json.loads((PARENT/'before_hashes.json').read_text())
modified = {r['file'].replace('\\', '/') for r in json.loads((PARENT/'document_changes.json').read_text())}
preserved = [rel for rel in manifest if rel.replace('\\', '/') not in modified]
for rel in preserved:
    assert digest(BASE/rel) == manifest[rel], rel
previous = 0
for folder in ['clinical_review_pass2', 'clinical_review_pass3', 'uncertainty_pass4']:
    for rel, value in json.loads((PARENT/folder/'artifact_hashes.json').read_text()).items():
        assert digest(PARENT/folder/rel) == value, (folder, rel)
        previous += 1
checks = {'favorable_cases_covered': len(triage), 'fresh_source_notes': len(notes),
          'review_template_targets': len(form), 'review_judgments_all_blank': True,
          'triage_counts': dict(Counter(r['working_status'] for r in triage)),
          'protected_historical_files_unchanged': len(preserved),
          'previous_manifest_files_unchanged': previous,
          'limits': 'Schema/coverage/integrity checks only; no independent clinical adjudication or new AUC.'}
(HERE/'checks.json').write_text(json.dumps(checks, indent=2)+'\n', encoding='utf-8')
files = sorted(p for p in HERE.iterdir() if p.is_file() and p.name != 'artifact_hashes.json')
(HERE/'artifact_hashes.json').write_text(json.dumps({p.name: digest(p) for p in files}, indent=2)+'\n')
print(json.dumps(checks, indent=2))
