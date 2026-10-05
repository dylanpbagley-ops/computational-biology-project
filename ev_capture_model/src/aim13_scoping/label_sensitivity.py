"""Does the headline result depend on Dylan's own safety-label judgements?

WHY THIS EXISTS
---------------
The 44 clinical safety labels in `target_labels.py` are the softest input in the project:
they are literature judgements made by one person, and every downstream statistic rests on
them. The standard defence is that the labels are frozen before computation and carry
citations. That is necessary but not sufficient -- it shows the labels were not fitted to
the answer, but it does not show the answer would survive if a reviewer disagreed with some
of them.

This script tests that directly, by re-running the section 23 whole-body test under
adversarial relabelings -- the specific reassignments a hostile reviewer is most likely to
propose, each defensible on its own terms:

  * NECTIN4 -> toxic     enfortumab vedotin carries real on-target skin effects; a reviewer
                         could argue they belong in the toxic class rather than the safe one.
  * CD33 -> toxic        gemtuzumab causes myelosuppression; called safe here because marrow
                         regenerates, but that is a judgement.
  * TACSTD2 -> toxic     sacituzumab's neutropenia is called payload-driven rather than
                         on-target; a reviewer could disagree.
  * EGFR -> safe         cetuximab is an approved drug; a reviewer could argue its skin
                         toxicity is manageable rather than disqualifying.

These are perturbations WITHIN the existing taxonomy (moving individual targets between
classes). They do NOT test the taxonomy itself -- whether "dose-limiting on-target toxicity"
is the right axis to classify on. That remains a genuine open limitation and should be
stated as such in the write-up rather than implied to be covered here.

RESULT SUMMARY (see print_summary for live numbers): the separation stays in a narrow band
across every relabeling tested and remains significant throughout. The finding does not
depend on these label calls.
"""
import itertools
import json
import os

import numpy as np

import target_labels as TL
from wholebody_v2 import sep, exact_p

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'wholebody_v2.json')
OUT = os.path.join(HERE, 'label_sensitivity.json')

# Each scenario: name -> {gene: new_label}. Empty dict = the frozen labels, unchanged.
SCENARIOS = {
    'as_frozen': {},
    'NECTIN4_toxic': {'NECTIN4': 'toxic'},
    'NECTIN4_CD33_toxic': {'NECTIN4': 'toxic', 'CD33': 'toxic'},
    'EGFR_safe': {'EGFR': 'safe'},
    'NECTIN4_CD33_TACSTD2_toxic': {'NECTIN4': 'toxic', 'CD33': 'toxic', 'TACSTD2': 'toxic'},
}


def _scores(A, two_factor):
    return {g: A[g]['top_share'] + (1.0 if (two_factor and TL.haematopoietic(A[g]['home'])) else 0.0)
            for g in A}


def run_scenario(A, base_labels, changes, two_factor):
    labels = dict(base_labels)
    labels.update(changes)
    scores = _scores(A, two_factor)
    a = [scores[g] for g in scores if labels.get(g) == 'safe']
    b = [scores[g] for g in scores if labels.get(g) == 'toxic']
    s = sep(a, b)
    p, exact = exact_p(a, b)
    return {'sep': float(s), 'p': float(p), 'p_is_exact': bool(exact),
            'n_safe': len(a), 'n_toxic': len(b)}


def compute_all():
    d = json.load(open(SRC))
    A = d['reference_A_cell_type']
    base = d['labels']
    out = {}
    for name, changes in SCENARIOS.items():
        out[name] = {
            'changes': changes,
            'one_factor': run_scenario(A, base, changes, two_factor=False),
            'two_factor': run_scenario(A, base, changes, two_factor=True),
        }
    seps = [v[k]['sep'] for v in out.values() for k in ('one_factor', 'two_factor')]
    ps = [v[k]['p'] for v in out.values() for k in ('one_factor', 'two_factor')]
    out['_summary'] = {
        'sep_min': float(min(seps)), 'sep_max': float(max(seps)),
        'p_max': float(max(ps)),
        'n_scenarios': len(SCENARIOS),
        'caveat': ('Perturbs labels WITHIN the existing taxonomy. Does not test the '
                   'taxonomy itself (whether dose-limiting on-target toxicity is the '
                   'right classification axis) -- that remains an open limitation.'),
    }
    return out


def print_summary():
    r = compute_all()
    print('=' * 78)
    print('LABEL SENSITIVITY: does the headline result depend on the safety-label calls?')
    print('=' * 78)
    print(f"{'scenario':30s} {'n':>7s}  {'1-factor':>18s}  {'2-factor':>18s}")
    for name in SCENARIOS:
        v = r[name]
        o, t = v['one_factor'], v['two_factor']
        n = f"{o['n_safe']}v{o['n_toxic']}"
        print(f"{name:30s} {n:>7s}  {o['sep']:.3f} (p={o['p']:.5f})  "
              f"{t['sep']:.3f} (p={t['p']:.5f})")
    s = r['_summary']
    print()
    print(f"separation range across all scenarios: {s['sep_min']:.3f} to {s['sep_max']:.3f}")
    print(f"worst p-value across all scenarios:    {s['p_max']:.5f}")
    print()
    print('CAVEAT:', s['caveat'])
    return r


if __name__ == '__main__':
    res = print_summary()
    with open(OUT, 'w') as f:
        json.dump(res, f, indent=1)
    print(f'\nwrote {OUT}')
