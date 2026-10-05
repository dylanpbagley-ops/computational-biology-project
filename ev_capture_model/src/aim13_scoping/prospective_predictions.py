"""
LOCKED PROSPECTIVE PREDICTIONS.

Everything in this project so far is retrospective: the clinical outcome was already known when
the measure was computed. The strongest possible answer to "did you fit this after the fact?" is
not another argument -- it is a dated, locked list of predictions about targets whose verdict
does not exist yet, published before the answer is known.

That is what this script produces.

THE PANEL. Targets that are (a) in registered clinical trials, from the De Schepper et al.
trial registry of 623 trials against 74 ADC targets, and (b) carry NO safety verdict in
target_labels.py -- that is, no severe on-target toxicity is documented and no approval has been
granted. These are the genuinely undecided cases.

THE RULE, fixed and identical to results section 23.2, with no free parameter:
    score = top_share + 1.0 if the home cell type is a haematopoietic lineage
    PREDICT TOLERATED   if score >= 0.25
    PREDICT ON-TARGET TOXICITY if score < 0.10
    NO CALL in between -- the framework is not claiming to resolve the middle.
The two thresholds come from where the 19 safe and 10 toxic targets actually sit (section 23.1)
and are written here once, before the undecided targets are scored.

WHAT WOULD FALSIFY IT. If, over the next several years, the targets predicted TOLERATED show
severe on-target toxicity at the same rate as those predicted TOXIC, the framework is wrong.
That is a real risk, taken deliberately.
"""
import os, json, collections, datetime
import numpy as np
import openpyxl
import target_labels as TL

HERE = os.path.dirname(os.path.abspath(__file__))
ADC = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'ev_capture_data', 'adc_autopsy'))
ALIAS = {'HER2': 'ERBB2', 'HER3': 'ERBB3', 'TROP2': 'TACSTD2', 'PDL1': 'CD274', 'CD73': 'NT5E'}

WEB = json.load(open(os.path.join(HERE, 'confinement_web.json'), encoding='utf-8'))
CT, HAEM = WEB['cell_types'], WEB['haem']
G = {r[0]: r for r in WEB['genes']}
LABEL = {t[0]: t[3] for t in TL.TARGETS}

SAFE_CUT, TOX_CUT = 0.25, 0.10

wb = openpyxl.load_workbook(os.path.join(ADC, '41467_2025_67840_MOESM3_ESM.xlsx'), read_only=True)
rows = list(wb['Sup data 1'].iter_rows(values_only=True))
wb.close()
hdr = rows[1]
gi, ph, dr, si = (hdr.index('Gene_Symbol'), hdr.index('Phase'),
                  hdr.index('Drug'), hdr.index('Tumor_Site'))
info = collections.defaultdict(lambda: {'phases': set(), 'drugs': set(), 'sites': set(), 'n': 0})
for r in rows[2:]:
    if not r[0]:
        continue
    g = ALIAS.get(str(r[gi]).strip(), str(r[gi]).strip())
    info[g]['phases'].add(str(r[ph]).strip())
    info[g]['drugs'].add(str(r[dr]).strip())
    info[g]['sites'].add(str(r[si]).strip())
    info[g]['n'] += 1

undecided = []
for g, d in info.items():
    if LABEL.get(g) in ('safe', 'toxic'):
        continue          # verdict already known -- not a prediction
    rec = G.get(g)
    if rec is None:
        continue
    share, ntypes, j = rec[1], rec[2], rec[3]
    score = share + (1.0 if HAEM[j] else 0.0)
    call = ('TOLERATED' if score >= SAFE_CUT else
            'ON-TARGET TOXICITY' if score < TOX_CUT else 'no call')
    undecided.append((score, g, share, ntypes, CT[j], bool(HAEM[j]), call, d))

undecided.sort(reverse=True)
today = datetime.date.today().isoformat()
print('LOCKED PREDICTIONS  --  generated %s' % today)
print('rule: score = confinement + 1.0 if the home is a blood lineage;'
      ' >= %.2f tolerated, < %.2f toxic' % (SAFE_CUT, TOX_CUT))
print('panel: targets in registered trials with NO safety verdict yet\n')
print('%-9s %7s %7s %-30s %-19s %5s %s'
      % ('GENE', 'score', 'n_types', 'home cell type', 'PREDICTION', 'trials', 'highest phase'))
for score, g, share, ntypes, home, haem, call, d in undecided:
    hp = max([p for p in d['phases'] if p and p != 'None'] or ['?'],
             key=lambda x: (len(x), x))
    print('%-9s %7.3f %7d %-30s %-19s %5d %s'
          % (g, score, ntypes, home[:30], call, d['n'], hp))

n_tol = sum(1 for u in undecided if u[6] == 'TOLERATED')
n_tox = sum(1 for u in undecided if u[6] == 'ON-TARGET TOXICITY')
n_no = sum(1 for u in undecided if u[6] == 'no call')
print('\n%d undecided targets: %d predicted tolerated, %d predicted toxic, %d no call'
      % (len(undecided), n_tol, n_tox, n_no))
print("""
HOW TO USE THIS. Print it, date it, and put it in the write-up unchanged. Do not revise a
prediction after a trial reports. If the list is wrong, that is the result, and it is worth
more to the field than another retrospective fit.""")

json.dump({'generated': today, 'rule': {'safe_cut': SAFE_CUT, 'toxic_cut': TOX_CUT,
                                        'formula': 'top_share + 1.0 if haematopoietic home'},
           'predictions': [{'gene': g, 'score': score, 'confinement': share,
                            'n_cell_types': ntypes, 'home': home, 'expendable_home': haem,
                            'prediction': call, 'n_trials': d['n'],
                            'tumour_sites': sorted(d['sites'])[:6]}
                           for score, g, share, ntypes, home, haem, call, d in undecided]},
          open(os.path.join(HERE, 'prospective_predictions.json'), 'w'), indent=1)
print('\nwrote prospective_predictions.json')
