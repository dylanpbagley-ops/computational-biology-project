"""Bootstrap confidence intervals for every headline separation (AUC) statistic.

WHY THIS EXISTS
---------------
Every headline number in sections 23/26/27 is reported as a separation value plus a
permutation p-value. A p-value asserts that an effect is unlikely under the null; it says
nothing about how precisely the effect size itself is pinned down. With 19 safe vs 11 toxic
targets, that precision is genuinely limited, and stating the limit explicitly is more
defensible than reporting a bare point estimate.

These intervals are expected to be WIDE. That is the honest picture and it is the point:
the aim is calibrated awareness of precision, not the appearance of precision.

METHOD
------
Stratified nonparametric bootstrap over TARGETS, which is the sampling unit (not cell types,
not cells). The safe group and the toxic group are resampled independently with replacement,
preserving each group's size, and the separation is recomputed on each resample. The reported
interval is the percentile interval of that bootstrap distribution.

Percentile intervals were chosen over BCa deliberately: with n this small the bias-correction
and acceleration terms are themselves estimated from very few points, so BCa buys accuracy
that the sample size cannot actually support, at the cost of an interval that is much harder
to explain and defend. A percentile interval is transparent and conservative here.

A degenerate resample (every drawn target identical, so separation is exactly 0.5) is a real
possible outcome at this sample size and is KEPT, not discarded — dropping it would narrow
the interval by removing exactly the cases that demonstrate how little data there is.

`sep` is imported from wholebody_v2 rather than reimplemented, so there is one and only one
definition of separation in the codebase.
"""
import json
import os
import numpy as np

import target_labels as TL
from wholebody_v2 import sep

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'wholebody_v2.json')
OUT = os.path.join(HERE, 'auc_confidence_intervals.json')

N_BOOT = 10000
SEED = 0
LEVEL = 95.0


def bootstrap_ci(a, b, n_boot=N_BOOT, seed=SEED, level=LEVEL):
    """Percentile bootstrap CI for sep(a, b), resampling each group independently.

    Returns (point, lo, hi, n_degenerate) where n_degenerate counts resamples in which
    every value drawn was identical across both groups (separation exactly 0.5).
    """
    a = np.asarray(list(a), float)
    b = np.asarray(list(b), float)
    rng = np.random.default_rng(seed)
    point = sep(a, b)
    stats = np.empty(n_boot, float)
    degenerate = 0
    for i in range(n_boot):
        aa = rng.choice(a, size=len(a), replace=True)
        bb = rng.choice(b, size=len(b), replace=True)
        stats[i] = sep(aa, bb)
        if np.all(aa == aa[0]) and np.all(bb == bb[0]) and aa[0] == bb[0]:
            degenerate += 1
    tail = (100.0 - level) / 2.0
    lo = float(np.percentile(stats, tail))
    hi = float(np.percentile(stats, 100.0 - tail))
    return float(point), lo, hi, degenerate


def _scores(A, H, two_factor):
    """Reproduce the section 23 scoring rule exactly: top_share, plus 1.0 if the home
    cell type is haematopoietic (two-factor only)."""
    return {g: A[g]['top_share'] + (1.0 if (two_factor and H.get(g)) else 0.0) for g in A}


def _split(scores, labels, pos='safe', neg='toxic'):
    a = [scores[g] for g in scores if labels.get(g) == pos]
    b = [scores[g] for g in scores if labels.get(g) == neg]
    return a, b


def compute_all():
    d = json.load(open(SRC))
    A = d['reference_A_cell_type']
    H = d['haem_home_A']
    labels = d['labels']
    results = {}

    for name, two in (('A_one_factor_top_share', False), ('A_two_factor', True)):
        a, b = _split(_scores(A, H, two), labels)
        pt, lo, hi, deg = bootstrap_ci(a, b)
        results[name] = {
            'point': pt, 'ci_lo': lo, 'ci_hi': hi,
            'n_safe': len(a), 'n_toxic': len(b),
            'n_degenerate_resamples': deg,
        }

    # Reference B is whole-TISSUE, so its home is a tissue, not a cell type. The
    # haematopoietic flag must be evaluated against B's own homes -- borrowing reference
    # A's cell-type flags here silently inflates the result (0.866 vs the correct 0.751),
    # because A's flag fires on far more targets than B's tissue-level home does.
    B = d['reference_B_tissue']
    Hb = {g: bool(TL.haematopoietic(B[g]['home'])) for g in B}
    a, b = _split(_scores(B, Hb, True), labels)
    pt, lo, hi, deg = bootstrap_ci(a, b)
    results['B_two_factor'] = {
        'point': pt, 'ci_lo': lo, 'ci_hi': hi,
        'n_safe': len(a), 'n_toxic': len(b),
        'n_degenerate_resamples': deg,
    }

    # --- solid tumours only (blood-cancer targets removed) --------------------------
    # The smallest and most fragile outcome test (8 v 10), so the interval matters more
    # here than anywhere else.
    #
    # NOTE: the blood set is the EXPLICIT, hand-curated list in wholebody_v2.py, not a
    # computed one. An earlier version of this code filtered on "non-haematopoietic home"
    # instead, which is a different criterion and gave 0.851 at n=11v11 rather than the
    # documented 0.8125 at n=8v10. Reconstructing a published subset by inference rather
    # than by reusing its actual definition is exactly how a number drifts.
    BLOOD = {'CD19', 'MS4A1', 'CD22', 'CD79B', 'TNFRSF17', 'SLAMF7', 'CD38', 'TNFRSF8',
             'CD33', 'GPRC5D', 'CD52', 'PTPRC', 'IL3RA', 'FCER2', 'CD40', 'CD80',
             'HLA-DRA', 'SDC1', 'CD74', 'CD47'}
    solid = {g: A[g] for g in A if g not in BLOOD}
    a, b = _split(_scores(solid, H, False), labels)
    pt, lo, hi, deg = bootstrap_ci(a, b)
    results['solid_only_one_factor'] = {
        'point': pt, 'ci_lo': lo, 'ci_hi': hi,
        'n_safe': len(a), 'n_toxic': len(b),
        'n_degenerate_resamples': deg,
        'documented_point': (d.get('solid_only') or {}).get('sep'),
    }

    results['_meta'] = {
        'n_boot': N_BOOT, 'seed': SEED, 'level': LEVEL,
        'method': 'stratified percentile bootstrap over targets',
        'note': ('Intervals are wide by construction at this sample size. Degenerate '
                 'resamples are retained, not discarded.'),
    }
    return results


def print_summary():
    r = compute_all()
    print('=' * 78)
    print('BOOTSTRAP CONFIDENCE INTERVALS ON HEADLINE SEPARATIONS')
    print(f"{r['_meta']['n_boot']} resamples, {r['_meta']['level']:.0f}% percentile, "
          f"stratified over targets, seed={r['_meta']['seed']}")
    print('=' * 78)
    for k, v in r.items():
        if k.startswith('_'):
            continue
        width = v['ci_hi'] - v['ci_lo']
        print(f"{k:28s} {v['point']:.3f}  95% CI [{v['ci_lo']:.3f}, {v['ci_hi']:.3f}]  "
              f"width={width:.3f}  (n={v['n_safe']}v{v['n_toxic']})")
    print()
    print('Read these as: the effect is real (see permutation p-values) but its MAGNITUDE')
    print('is pinned down only to within the interval shown. Report both.')
    return r


if __name__ == '__main__':
    res = print_summary()
    with open(OUT, 'w') as f:
        json.dump(res, f, indent=1)
    print(f'\nwrote {OUT}')
