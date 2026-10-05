"""Does the section 23 result survive replacing the binary flag with continuous weights?

Runs the whole-body safe-vs-toxic test under several weighting schemes and reports the
result under each. The deliverable is the SPREAD, not any single number.

Schemes:
  binary             the current rule (+1.0 for a haematopoietic home) -- baseline
  continuous         per-cell-type weights from expendability_weights.WEIGHTS
  regeneration_only  a coarser variant using only "does this tissue renew?"
  shuffled_*         null controls: the SAME continuous weights, randomly reassigned
                     across cell types. These destroy the biological meaning while keeping
                     the numeric distribution identical, so whatever separation survives
                     is what the scheme's arithmetic buys on its own.

Decision recorded BEFORE running (per the plan in Expansion_Plan_Five_Angles_Sept26.md):
the binary result stays primary unless the continuous scheme performs at least as well,
and any decrease is reported rather than buried.
"""
import json
import os

import numpy as np

import expendability_weights as W
from wholebody_v2 import sep, exact_p

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'wholebody_v2.json')
OUT = os.path.join(HERE, 'continuous_weight_sensitivity.json')

N_SHUFFLE = 200
SHUFFLE_SEED = 20260927


def run(A, labels, scheme, weight_fn=None, with_p=True):
    if weight_fn is None:
        scores = {g: W.score(A[g]['top_share'], A[g]['home'], scheme) for g in A}
    else:
        scores = {g: A[g]['top_share'] + (1.0 - weight_fn(A[g]['home'])) for g in A}
    a = [scores[g] for g in scores if labels.get(g) == 'safe']
    b = [scores[g] for g in scores if labels.get(g) == 'toxic']
    s = sep(a, b)
    if not with_p:
        # the shuffled null needs only the separation; computing a permutation p-value
        # inside each of N_SHUFFLE scrambles is ~200x wasted work
        return {'sep': float(s), 'n_safe': len(a), 'n_toxic': len(b)}
    p, exact = exact_p(a, b)
    return {'sep': float(s), 'p': float(p), 'p_is_exact': bool(exact),
            'n_safe': len(a), 'n_toxic': len(b)}


def compute_all():
    d = json.load(open(SRC))
    A = d['reference_A_cell_type']
    labels = d['labels']
    out = {}

    for scheme in ('binary', 'continuous', 'regeneration_only'):
        out[scheme] = run(A, labels, scheme)

    # --- null: same weight VALUES, scrambled across cell types -------------------------
    rng = np.random.default_rng(SHUFFLE_SEED)
    types = sorted({A[g]['home'] for g in A})
    vals = [W.weight(t, 'continuous') for t in types]
    seps = []
    for _ in range(N_SHUFFLE):
        perm = rng.permutation(vals)
        m = dict(zip(types, perm))
        seps.append(run(A, labels, 'continuous', weight_fn=lambda h: m[h],
                        with_p=False)['sep'])
    seps = np.array(seps)
    obs = out['continuous']['sep']
    out['shuffled_null'] = {
        'n_shuffles': N_SHUFFLE,
        'sep_mean': float(seps.mean()),
        'sep_p95': float(np.percentile(seps, 95)),
        'sep_max': float(seps.max()),
        'observed_continuous': obs,
        'p_vs_null': float((np.sum(seps >= obs - 1e-12) + 1) / (N_SHUFFLE + 1)),
    }

    real = [out[s]['sep'] for s in ('binary', 'continuous', 'regeneration_only')]
    out['_summary'] = {
        'sep_min': float(min(real)), 'sep_max': float(max(real)),
        'p_max': float(max(out[s]['p'] for s in ('binary', 'continuous', 'regeneration_only'))),
        'continuous_minus_binary': float(out['continuous']['sep'] - out['binary']['sep']),
    }
    return out


def print_summary():
    r = compute_all()
    print('=' * 78)
    print('CONTINUOUS EXPENDABILITY WEIGHTS vs THE BINARY FLAG (section 23 test)')
    print('=' * 78)
    print(f"{'scheme':22s} {'separation':>11s} {'p':>11s}")
    for s in ('binary', 'continuous', 'regeneration_only'):
        v = r[s]
        print(f"{s:22s} {v['sep']:11.3f} {v['p']:11.5f}")
    n = r['shuffled_null']
    print()
    print(f"shuffled null ({n['n_shuffles']} scrambles of the same weight values):")
    print(f"  mean {n['sep_mean']:.3f}   95th pct {n['sep_p95']:.3f}   max {n['sep_max']:.3f}")
    print(f"  observed continuous {n['observed_continuous']:.3f}, p vs null = {n['p_vs_null']:.4f}")
    s = r['_summary']
    print()
    print(f"continuous - binary = {s['continuous_minus_binary']:+.3f}")
    return r


if __name__ == '__main__':
    res = print_summary()
    with open(OUT, 'w') as f:
        json.dump(res, f, indent=1)
    print(f'\nwrote {OUT}')
