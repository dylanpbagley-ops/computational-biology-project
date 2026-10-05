"""Filter an existing exhaustive table using disclosed, source-based constraints."""
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
with (HERE/'all_1024_assignments.tsv').open(encoding='utf-8', newline='') as stream:
    all_rows = list(csv.DictReader(stream, delimiter='\t'))
rows = [r for r in all_rows if not {'CD47', 'NECTIN4'} & set(r['uncertain_assigned_favorable'].split(';'))]
assert len(all_rows) == 1024 and len(rows) == 256
metrics = ['frozen_rule', 'confinement', 'blood_flag', 'rule_minus_confinement', 'rule_minus_blood']
for r in rows:
    for m in metrics:
        r[m] = float(r[m])
out = {'status': 'Post hoc source-constrained sensitivity, not adjudication',
       'forced_unfavorable': ['CD47', 'NECTIN4'], 'assignments': len(rows),
       'envelopes': {m: {'min': min(rows, key=lambda r: r[m]), 'max': max(rows, key=lambda r: r[m])} for m in metrics}}
(HERE/'constrained_results.json').write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps({m: [out['envelopes'][m]['min'][m], out['envelopes'][m]['max'][m]] for m in metrics}, indent=2))
