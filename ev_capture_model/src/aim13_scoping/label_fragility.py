"""Label fragility: how many labels must be wrong to erase the frozen rule's result? Precommitted in
label_fragility_precommitment.md. Rule and scores unchanged; only labels are swapped."""
import itertools, json, math, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'label_fragility.json')
WEAK_TO_SAFE = ['EGFR', 'EPCAM', 'NCAM1', 'FAP', 'CD47', 'ERBB2', 'CEACAM5']
WEAK_TO_TOX = ['NECTIN4']
K_MAX = 6


def load():
    d = json.load(open(os.path.join(HERE, 'flag_cleanup.json')))['targets']
    A = json.load(open(os.path.join(HERE, 'wholebody_v2.json')))['reference_A_cell_type']  # unrounded top_share
    genes = sorted(d)
    score = np.array([A[g]['top_share'] + (1.0 if d[g]['M2_flag'] else 0.0) for g in genes])
    lab = np.array([1 if d[g]['label'] == 'safe' else 0 for g in genes])
    return genes, score, lab


def rank_stats(score, is_safe):
    """Return separation and MW normal-approx two-sided-in-direction (one-sided safe>toxic) p, tie-corrected."""
    n1, n0 = int(is_safe.sum()), int((~is_safe).sum())
    if n1 == 0 or n0 == 0:
        return float('nan'), 1.0
    order = score.argsort(kind='mergesort')
    s = score[order]
    ranks = np.empty(len(s)); i = 0; tie_term = 0.0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        ranks[i:j + 1] = (i + j) / 2 + 1
        t = j - i + 1; tie_term += t ** 3 - t
        i = j + 1
    r = np.empty(len(s)); r[order] = ranks
    R1 = r[is_safe].sum()
    U = R1 - n1 * (n1 + 1) / 2
    auc = U / (n1 * n0)
    N = n1 + n0
    var = n1 * n0 / 12.0 * ((N + 1) - tie_term / (N * (N - 1)))
    z = (U - n1 * n0 / 2 - 0.5) / math.sqrt(var)
    p = 0.5 * math.erfc(z / math.sqrt(2))
    if tie_term == 0:
        p = exact_p(auc, n1, n0)
    return auc, p


_TAB = {}
def exact_table(n1, n0):
    """P(U >= u) for the Mann-Whitney U, no ties: counts of arrangements via DP on partitions."""
    key = (n1, n0)
    if key in _TAB:
        return _TAB[key]
    # f[m][n] as polynomial in U: use recurrence c(n1,n0,u) = c(n1-1,n0,u-n0) + c(n1,n0-1,u)
    from functools import lru_cache
    import sys; sys.setrecursionlimit(10000)
    maxu = n1 * n0
    # iterative DP over (a,b)
    dp = {}
    for a in range(n1 + 1):
        for b in range(n0 + 1):
            arr = np.zeros(a * b + 1)
            if a == 0 or b == 0:
                arr[0] = 1
            else:
                x = dp[(a - 1, b)]; y = dp[(a, b - 1)]
                arr[:len(x)] += 0
                arr_x = np.zeros(a * b + 1); arr_x[b:b + len(x)] = x
                arr_y = np.zeros(a * b + 1); arr_y[:len(y)] = y
                arr = arr_x + arr_y
            dp[(a, b)] = arr
    c = dp[(n1, n0)]
    tail = np.cumsum(c[::-1])[::-1] / c.sum()
    _TAB[key] = tail
    return tail


def exact_p(auc, n1, n0):
    u = int(round(auc * n1 * n0))
    return float(exact_table(n1, n0)[u])


def mc_p(score, is_safe, n=20000, seed=7):
    rng = np.random.default_rng(seed)
    a0, _ = rank_stats(score, is_safe)
    n1 = int(is_safe.sum()); cnt = 0
    for _ in range(n):
        perm = np.zeros(len(score), bool); perm[rng.choice(len(score), n1, replace=False)] = True
        a, _ = rank_stats(score, perm)
        cnt += a >= a0 - 1e-12
    return (cnt + 1) / (n + 1)


def main():
    genes, score, lab = load()
    is_safe = lab == 1
    auc0, p0 = rank_stats(score, is_safe)
    res = {'frozen': {'sep': auc0, 'p_normal_approx': p0, 'p_montecarlo': mc_p(score, is_safe)}}
    # validation on random flip sets
    rng = np.random.default_rng(11); worst_ratio = 1.0; val = []
    for _ in range(30):
        k = int(rng.integers(1, 5)); idx = rng.choice(30, k, replace=False)
        l2 = is_safe.copy(); l2[idx] = ~l2[idx]
        a, p = rank_stats(score, l2); pm = mc_p(score, l2, 5000, 3)
        val.append((p, pm)); worst_ratio = max(worst_ratio, max(p, pm) / max(min(p, pm), 1e-12))
    res['approx_validation'] = {'max_ratio_p_normal_vs_mc': worst_ratio, 'n': 30,
                                'note': 'MC has floor ~1/5001; ratios near the floor are not informative'}
    # F1 exhaustive
    F1 = {}
    for k in range(1, K_MAX + 1):
        seps, ps = [], []; worst = None
        for idx in itertools.combinations(range(30), k):
            l2 = is_safe.copy(); l2[list(idx)] = ~l2[list(idx)]
            a, p = rank_stats(score, l2)
            seps.append(a); ps.append(p)
            if worst is None or p > worst[0]:
                worst = (p, a, [genes[i] for i in idx])
        seps, ps = np.array(seps), np.array(ps)
        F1[k] = {'n_sets': len(seps), 'sep_min': float(seps.min()), 'sep_median': float(np.median(seps)),
                 'sep_max': float(seps.max()), 'frac_p_lt_0.05': float((ps < 0.05).mean()),
                 'frac_p_lt_0.001': float((ps < 0.001).mean()),
                 'worst_p': worst[0], 'worst_sep': worst[1], 'worst_set': worst[2]}
        print(k, F1[k], flush=True)
    res['F1'] = F1
    # which 4-swap sets remove significance: do they all move an approved blood target?
    blood_safe = [g for i, g in enumerate(genes) if is_safe[i] and score[i] >= 1.0]
    bad = []
    for idx in itertools.combinations(range(30), 4):
        l2 = is_safe.copy(); l2[list(idx)] = ~l2[list(idx)]
        if rank_stats(score, l2)[1] >= 0.05:
            bad.append([genes[i] for i in idx])
    res['k4_sets_losing_significance'] = {'n': len(bad), 'blood_safe_targets': blood_safe,
        'n_containing_a_blood_safe_target': sum(any(g in blood_safe for g in b) for b in bad)}
    fi = next((k for k in range(1, K_MAX + 1) if F1[k]['frac_p_lt_0.05'] < 1.0), None)
    res['fragility_index'] = fi if fi else f'> {K_MAX}'
    # F2 weak set
    weak = [(g, 'safe') for g in WEAK_TO_SAFE] + [(g, 'toxic') for g in WEAK_TO_TOX]
    gi = {g: i for i, g in enumerate(genes)}
    rows = []
    for mask in range(256):
        l2 = is_safe.copy(); ch = []
        for b, (g, new) in enumerate(weak):
            if mask >> b & 1:
                l2[gi[g]] = (new == 'safe'); ch.append(g)
        a, p = rank_stats(score, l2)
        rows.append({'flipped': ch, 'n_flipped': len(ch), 'sep': a, 'p': p, 'n_safe': int(l2.sum())})
    res['F2_all_subsets'] = rows
    res['F2_summary'] = {'sep_min': min(r['sep'] for r in rows), 'sep_max': max(r['sep'] for r in rows),
                         'frac_p_lt_0.05': float(np.mean([r['p'] < 0.05 for r in rows])),
                         'all_eight': next(r for r in rows if r['n_flipped'] == 8),
                         'by_n_flipped_min_sep': {n: min(r['sep'] for r in rows if r['n_flipped'] == n) for n in range(9)},
                         'by_n_flipped_max_p': {n: max(r['p'] for r in rows if r['n_flipped'] == n) for n in range(9)}}
    json.dump(res, open(OUT, 'w'), indent=1, default=str)
    print('frozen', res['frozen']); print('fragility_index', res['fragility_index']); print('F2', res['F2_summary'])
    print('approx validation', res['approx_validation'])


if __name__ == '__main__':
    main()
