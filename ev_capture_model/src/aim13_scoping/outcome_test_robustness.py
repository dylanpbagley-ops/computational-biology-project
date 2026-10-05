"""Strengthening pass on the outcome tests: three questions the point estimates cannot answer.

1. HOW PRECISE are the two independent-data results (autopsy replication, registry control)?
   -> bootstrap CIs, same method as auc_confidence_intervals.py.

2. COULD THE SPECIFICITY CONTROL HAVE FAILED? The registry control is reported as "passes:
   separation 0.495, p = 0.53". A control only means something if it was capable of detecting the
   effect it is ruling out. With 9 vs 11 targets it may not be. This computes its statistical
   power at several true effect sizes (binormal model, one-sided Mann-Whitney at alpha = 0.05),
   so the write-up can state what the control does and does NOT exclude.

3. DOES THE EXPENDABILITY TERM ADD ANYTHING BEYOND CONFINEMENT? The headline claim is the
   two-factor rule (confinement x expendability), but the one-factor number alone is already
   0.880. The claim that the second factor helps is a claim about a DIFFERENCE (0.914 - 0.880),
   and a difference needs its own uncertainty. A PAIRED bootstrap is used (the same resampled
   targets score under both rules), because the two scores are computed on the same targets and
   an unpaired comparison of two intervals would throw that information away.

Everything is reported as-is. Nothing here is tuned; wide or unflattering results are kept.
"""
import json
import os
import numpy as np
from scipy import stats

import target_labels as TL
from wholebody_v2 import sep
from auc_confidence_intervals import bootstrap_ci, _scores, _split

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'outcome_test_robustness.json')
N_SIM = 4000
SEED = 20260928


# ---------------------------------------------------------------- 1. replication groups
def replication_groups():
    W = json.load(open(os.path.join(HERE, 'wholebody_v2.json')))
    conf = {g: d['top_share'] for g, d in W['reference_A_cell_type'].items()}
    label = {t[0]: t[3] for t in TL.TARGETS}

    # The documented 17-target replication (section 27.2, adc_autopsy_ceiling.head_to_head) is
    # defined on the WITHIN-PATIENT estimates, which require >= 5 patients per gene. The
    # 'pooled' key admits one extra labelled target (18) and gives 0.778, not the documented
    # 0.750. Reusing the real definition, not a plausible-looking one, is what stops drift.
    ceil = json.load(open(os.path.join(HERE, 'adc_autopsy_ceiling.json')))['within']
    shared = [g for g in ceil if g in conf and label.get(g) in ('safe', 'toxic')]
    autopsy_a = [conf[g] for g in shared if label[g] == 'safe']
    autopsy_b = [conf[g] for g in shared if label[g] == 'toxic']

    reg = json.load(open(os.path.join(HERE, 'trial_registry_control.json')))
    late = reg['late_phase']
    reg_a = [conf[g] for g in late if late[g] and g in conf]
    reg_b = [conf[g] for g in late if (not late[g]) and g in conf]
    return (autopsy_a, autopsy_b), (reg_a, reg_b)


# ---------------------------------------------------------------- 2. power of a control
def power_at_auc(true_auc, n_a, n_b, n_sim=N_SIM, seed=SEED, alpha=0.05):
    """Probability that a one-sided Mann-Whitney test at `alpha` detects a true separation of
    `true_auc` with group sizes n_a vs n_b. Binormal model: a ~ N(d, 1), b ~ N(0, 1) with
    d = sqrt(2) * Phi^-1(true_auc), which gives P(a > b) = true_auc exactly."""
    d = np.sqrt(2.0) * stats.norm.ppf(true_auc)
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(n_sim):
        a = rng.normal(d, 1.0, n_a)
        b = rng.normal(0.0, 1.0, n_b)
        p = stats.mannwhitneyu(a, b, alternative='greater', method='exact').pvalue
        hits += p < alpha
    return hits / n_sim


def min_detectable_auc(n_a, n_b, target_power=0.80, grid=None):
    grid = grid or [round(x, 2) for x in np.arange(0.55, 0.99, 0.02)]
    for t in grid:
        if power_at_auc(t, n_a, n_b, n_sim=1500) >= target_power:
            return t
    return None


# ---------------------------------------------------------------- 3. paired difference
def paired_difference_ci(n_boot=10000, seed=SEED, level=95.0):
    d = json.load(open(os.path.join(HERE, 'wholebody_v2.json')))
    A, H, labels = d['reference_A_cell_type'], d['haem_home_A'], d['labels']
    s1, s2 = _scores(A, H, False), _scores(A, H, True)
    safe = [g for g in A if labels.get(g) == 'safe']
    tox = [g for g in A if labels.get(g) == 'toxic']
    point = sep([s2[g] for g in safe], [s2[g] for g in tox]) - \
        sep([s1[g] for g in safe], [s1[g] for g in tox])
    rng = np.random.default_rng(seed)
    diffs = np.empty(n_boot)
    for i in range(n_boot):
        sa = rng.choice(safe, len(safe), replace=True)
        sb = rng.choice(tox, len(tox), replace=True)
        diffs[i] = (sep([s2[g] for g in sa], [s2[g] for g in sb]) -
                    sep([s1[g] for g in sa], [s1[g] for g in sb]))
    tail = (100.0 - level) / 2.0
    return {
        'point_difference': float(point),
        'ci_lo': float(np.percentile(diffs, tail)),
        'ci_hi': float(np.percentile(diffs, 100.0 - tail)),
        'p_boot_difference_le_0': float(np.mean(diffs <= 0)),
        'p_boot_difference_lt_0': float(np.mean(diffs < 0)),
        'frac_exactly_zero': float(np.mean(diffs == 0)),
        'n_safe': len(safe), 'n_toxic': len(tox),
    }


def expendability_tests(n_perm=20000, seed=SEED):
    """Direct tests of the EXPENDABILITY claim, separate from the AUC increment.

    The AUC increment (+0.033) is a weak lens on the second factor because confinement already
    ranks most safe targets above most toxic ones; the flag can only rescue the few pairs still
    mis-ordered. The sharper question is whether a haematopoietic home is associated with
    tolerability at all:
      (a) Fisher exact test on the 2x2 table (haem home vs safe/toxic);
      (b) a permutation null for the AUC increment: shuffle which targets carry the flag
          (keeping its count fixed) and ask how often the increment is as large as observed.
    """
    d = json.load(open(os.path.join(HERE, 'wholebody_v2.json')))
    A, H, labels = d['reference_A_cell_type'], d['haem_home_A'], d['labels']
    genes = [g for g in A if labels.get(g) in ('safe', 'toxic')]
    flag = np.array([bool(H.get(g)) for g in genes])
    safe = np.array([labels[g] == 'safe' for g in genes])
    conf = np.array([A[g]['top_share'] for g in genes])
    table = [[int((flag & safe).sum()), int((~flag & safe).sum())],
             [int((flag & ~safe).sum()), int((~flag & ~safe).sum())]]
    fisher_p = stats.fisher_exact(table, alternative='greater')[1]

    def inc(f):
        s2 = conf + f.astype(float)
        return sep(s2[safe], s2[~safe]) - sep(conf[safe], conf[~safe])
    obs = inc(flag)
    rng = np.random.default_rng(seed)
    null = np.array([inc(rng.permutation(flag)) for _ in range(n_perm)])
    rescued = sorted(g for g, f, s_, c in zip(genes, flag, safe, conf) if f and s_ and
                     any(c <= cc for cc, ss in zip(conf, safe) if not ss))
    return {'table_[[haem&safe, nonhaem&safe],[haem&toxic, nonhaem&toxic]]': table,
            'fisher_exact_one_sided_p': float(fisher_p),
            'observed_auc_increment': float(obs),
            'perm_p_increment_at_least_observed': float(np.mean(null >= obs - 1e-12)),
            'flagged_safe_targets_that_confinement_alone_ranks_below_some_toxic_target': rescued}


def compute_all():
    (aa, ab), (ra, rb) = replication_groups()
    res = {}
    for name, (a, b) in (('autopsy_replication', (aa, ab)), ('registry_control_phase3', (ra, rb))):
        pt, lo, hi, deg = bootstrap_ci(a, b)
        res[name] = {'point': pt, 'ci_lo': lo, 'ci_hi': hi, 'n_a': len(a), 'n_b': len(b),
                     'n_degenerate_resamples': deg}

    pw = {}
    for name, (na, nb) in (('registry_control_phase3', (len(ra), len(rb))),
                           ('autopsy_replication', (len(aa), len(ab)))):
        pw[name] = {'n_a': na, 'n_b': nb,
                    'power_at_auc': {str(t): power_at_auc(t, na, nb)
                                     for t in (0.50, 0.60, 0.65, 0.75, 0.88)},
                    'min_detectable_auc_80pct': min_detectable_auc(na, nb)}
    res['power'] = pw
    res['paired_two_factor_minus_one_factor'] = paired_difference_ci()
    res['expendability_direct_tests'] = expendability_tests()
    res['_meta'] = {'n_sim': N_SIM, 'seed': SEED,
                    'power_model': 'binormal, one-sided exact Mann-Whitney, alpha=0.05'}
    return res


def print_summary(r):
    print('=' * 78)
    print('OUTCOME-TEST ROBUSTNESS')
    print('=' * 78)
    for k in ('autopsy_replication', 'registry_control_phase3'):
        v = r[k]
        print(f"{k:26s} {v['point']:.3f}  95% CI [{v['ci_lo']:.3f}, {v['ci_hi']:.3f}]  "
              f"(n={v['n_a']}v{v['n_b']})")
    print('\nPOWER (probability the test detects a true separation of this size):')
    for k, v in r['power'].items():
        row = '  '.join(f"AUC {t}: {p:.2f}" for t, p in v['power_at_auc'].items())
        print(f"  {k:26s} n={v['n_a']}v{v['n_b']}   {row}   "
              f"| min detectable @80%: {v['min_detectable_auc_80pct']}")
    p = r['paired_two_factor_minus_one_factor']
    print(f"\nPAIRED two-factor minus one-factor: {p['point_difference']:+.3f}  "
          f"95% CI [{p['ci_lo']:+.3f}, {p['ci_hi']:+.3f}]  "
          f"P(diff<=0)={p['p_boot_difference_le_0']:.3f}  "
          f"P(diff exactly 0)={p['frac_exactly_zero']:.3f}")


def print_expendability(r):
    e = r['expendability_direct_tests']
    print('\nEXPENDABILITY DIRECT TESTS')
    for k, v in e.items():
        print(f'  {k}: {v}')


if __name__ == '__main__':
    r = compute_all()
    print_summary(r)
    print_expendability(r)
    json.dump(r, open(OUT, 'w'), indent=1)
    print('\nwrote', OUT)
