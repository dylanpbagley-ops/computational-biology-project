"""
AN INDEPENDENCE CHECK USING SOMEBODY ELSE'S TABLE.

The standing weakness of the outcome tests is that the safe/toxic labels are literature
judgements made by one person (§23.8 item 3). This script uses a label that person did not
make: the trial registry compiled by De Schepper et al. (Nature Communications 2025,
Supplementary Data 1) -- 623 registered clinical trials against 74 antibody-drug-conjugate
targets, with the phase of each.

It is important to be clear about what this table can and cannot test.

  IT CANNOT be used as a safety label. "Reached phase 3" measures development progress --
  efficacy, funding, commercial interest -- not on-target toxicity. HER2 has 36 phase-3 trials
  and also caused a fatal on-target CAR-T event.

  IT CAN be used as a SPECIFICITY CONTROL, and this is stated as a prediction before the run:
  confinement should NOT predict which targets reached phase 3. If it did, that would suggest
  confinement is really tracking "targets people invest in" rather than "targets that do not
  poison normal tissue", and the whole outcome result would be suspect.

  PREDICTION: separation near 0.5, not significant.
"""
import os, csv, json, collections
import numpy as np
import openpyxl
import target_labels as TL
import wholebody_v2 as WB

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'ev_capture_data', 'adc_autopsy'))
ALIAS = {'HER2': 'ERBB2', 'HER3': 'ERBB3', 'TROP2': 'TACSTD2', 'PDL1': 'CD274', 'CD73': 'NT5E'}

wb = openpyxl.load_workbook(os.path.join(D, '41467_2025_67840_MOESM3_ESM.xlsx'), read_only=True)
rows = list(wb['Sup data 1'].iter_rows(values_only=True))
wb.close()
hdr = rows[1]
gi, ph = hdr.index('Gene_Symbol'), hdr.index('Phase')
late = collections.defaultdict(bool)
ntrials = collections.Counter()
for r in rows[2:]:
    if not r[0]:
        continue
    g = ALIAS.get(str(r[gi]).strip(), str(r[gi]).strip())
    p = str(r[ph]).strip()
    ntrials[g] += 1
    if p in ('3', '2/3', '4'):
        late[g] = True

W = json.load(open(os.path.join(HERE, 'wholebody_v2.json')))
conf = {g: d['top_share'] for g, d in W['reference_A_cell_type'].items()}
label = {t[0]: t[3] for t in TL.TARGETS}

shared = sorted(g for g in ntrials if g in conf)
print('targets with both a confinement value and registered trials: %d' % len(shared))
print('  reached phase 3+: %d   did not: %d'
      % (sum(1 for g in shared if late[g]), sum(1 for g in shared if not late[g])))

lab2 = {g: ('safe' if late[g] else 'toxic') for g in shared}   # reusing the two-group machinery
s, p, ex, na, nb = WB.run_test({g: conf[g] for g in shared}, lab2, True)
print('\n  confinement vs "reached phase 3"   separation %.3f  %s p = %.4f  (n=%d vs %d)'
      % (s, 'exact' if ex else 'approx', p, na, nb))
print('  PREDICTED: near 0.5 and not significant  ->  %s'
      % ('AS PREDICTED' if (abs(s - 0.5) < 0.15 and p > 0.05) else 'NOT AS PREDICTED'))

print('\n  for comparison, on the targets in this same registry that also carry a SAFETY label:')
use = {g: conf[g] for g in shared if label.get(g) in ('safe', 'toxic')}
if len({label[g] for g in use}) == 2:
    s2, p2, ex2, na2, nb2 = WB.run_test(use, label, True)
    print('    confinement vs on-target toxicity  separation %.3f  %s p = %.4f  (n=%d vs %d)'
          % (s2, 'exact' if ex2 else 'approx', p2, na2, nb2))

print("""
  READING. The two rows answer different questions on overlapping targets. Confinement is
  uninformative about whether a programme advanced, and informative about whether it poisoned
  normal tissue. That is the behaviour a toxicity measure should have, and it is checked here
  against a table this project did not compile.""")

json.dump({'n_targets': len(shared), 'sep_phase3': s, 'p_phase3': p,
           'late_phase': {g: bool(late[g]) for g in shared}},
          open(os.path.join(HERE, 'trial_registry_control.json'), 'w'), indent=1)
