"""Is the 'ceiling in clinical units' measuring antigen biology, or measurement noise?

THE QUESTION
------------
ceiling_in_clinical_units.py reports, at harm budgets of 1% and 5%, the maximum fraction of
tumour cells killable. The numbers are small (0.6-1.9% at a 1% budget) and approved targets
(BCMA 1.0%, CD38 0.6%, SLAMF7 0.8%) do not visibly outrank the panel's control channels. That
pattern has a specific explanation to test: a harm budget of b sits at the (1-b) quantile of the
normal population. If that quantile lies in a contaminated tail (doublets, ambient antibody,
sticky background), every marker's threshold is set by noise, and 'kill' collapses to roughly
the budget itself no matter what the antigen is. Then the metric carries no antigen information
at that budget.

THE TEST (pre-stated)
---------------------
  1. Control channels. The panel includes FITC and PE, which are fluorophore/isotype-type tags
     with no target biology. Their 'ceiling' is an empirical noise floor.
  2. Permutation null. For every marker, tumour-vs-compartment labels are shuffled within each
     sample and the same statistic recomputed. Under the null the expected kill equals the
     budget; the spread of the null gives an envelope.
  3. Excess over chance = (kill - budget) / (1 - budget); 0 means indistinguishable from a
     coin flip at that budget, 1 means the whole tumour is killable.
  4. Budgets 1%, 5%, 10%, 20%, 35%, 50% are all reported, so the budget at which biology
     emerges from noise -- if any -- is visible rather than assumed.

A marker is called 'above floor' at a budget only if its ceiling exceeds BOTH the best control
channel AND the 95th percentile of its own permutation null.
"""
import json
import os
import numpy as np

import ceiling_in_clinical_units as C

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'clinical_units_noise_floor.json')

BUDGETS = (0.01, 0.05, 0.10, 0.20, 0.35, 0.50)
CONTROLS = ('FITC', 'PE')
N_PERM = 60
SEED = 20260928


def _ceiling(samples, ix, name, budget, perm_rng=None):
    """Median-over-samples kill per compartment, then the minimum over compartments (the binding
    one). With perm_rng, tumour/compartment labels are shuffled within each sample first."""
    per_comp = {}
    for lab in ('T', 'MY', 'NK', 'NB'):
        vals = []
        for s in samples:
            if lab not in s['comps']:
                continue
            t = s['tumor'][:, ix[name]]
            n = s['comps'][lab][:, ix[name]]
            if perm_rng is not None:
                pool = np.concatenate([t, n])
                perm_rng.shuffle(pool)
                t, n = pool[:len(t)], pool[len(t):]
            vals.append(C.max_kill_at_harm(t, n, budget))
        if vals:
            per_comp[lab] = float(np.median(vals))
    return min(per_comp.values()) if per_comp else None


def compute_all(n_perm=N_PERM):
    names, ix, samples = C.load_populations()
    alt = {m: C.load_populations(exclude_gate_marker=m)[2] for m in C.GATE_MARKERS}
    rng = np.random.default_rng(SEED)
    res = {}
    for name in names:
        use = alt[name] if name in C.GATE_MARKERS else samples
        row = {}
        for b in BUDGETS:
            obs = _ceiling(use, ix, name, b)
            null = np.array([_ceiling(use, ix, name, b, rng) for _ in range(n_perm)])
            row[f'{b:.2f}'] = {'ceiling': obs, 'excess_over_chance': (obs - b) / (1 - b),
                               'null_mean': float(null.mean()),
                               'null_p95': float(np.percentile(null, 95))}
        res[name] = row
    summary = {}
    for b in BUDGETS:
        k = f'{b:.2f}'
        ctrl = max(res[c][k]['ceiling'] for c in CONTROLS)
        above = sorted(n for n in names if n not in CONTROLS and
                       res[n][k]['ceiling'] > ctrl and
                       res[n][k]['ceiling'] > res[n][k]['null_p95'])
        summary[k] = {'best_control_ceiling': ctrl,
                      'median_marker_ceiling': float(np.median([res[n][k]['ceiling'] for n in names])),
                      'n_markers': len(names) - len(CONTROLS),
                      'n_above_floor': len(above), 'above_floor': above}
    return {'per_marker': res, 'summary': summary,
            'meta': {'budgets': list(BUDGETS), 'controls': list(CONTROLS),
                     'n_perm': n_perm, 'seed': SEED}}


def print_summary(r):
    print('=' * 92)
    print('NOISE FLOOR OF THE CLINICAL-UNITS CEILING  (control channels FITC, PE + permutation null)')
    print('=' * 92)
    print(f"{'budget':>7s} {'best control':>13s} {'median marker':>14s} {'#above floor':>13s}   examples")
    for k, v in r['summary'].items():
        print(f"{k:>7s} {v['best_control_ceiling']:13.3f} {v['median_marker_ceiling']:14.3f} "
              f"{v['n_above_floor']:6d}/{v['n_markers']:<6d}   {', '.join(v['above_floor'][:8])}")
    print('\nKey targets (ceiling / excess-over-chance):')
    for m in ('CD269', 'CD38', 'CD319', 'CD307', 'CD138'):
        cells = '  '.join(f"{k}: {r['per_marker'][m][k]['ceiling']:.3f}/"
                          f"{r['per_marker'][m][k]['excess_over_chance']:+.3f}"
                          for k in r['summary'])
        print(f'  {m:6s} {cells}')


if __name__ == '__main__':
    r = compute_all()
    print_summary(r)
    json.dump(r, open(OUT, 'w'), indent=1)
    print('\nwrote', OUT)
