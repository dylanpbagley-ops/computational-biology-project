"""Lemma linking top_share to exposure breadth (Sept 29, 2026). A proof plus a numerical check.

Setup. A gene has non-negative expression x_c in each of N cell types; peak = max_c x_c; top_share t = peak / sum_c x_c.
For a relative exposure level r in (0,1], call a cell type EXPOSED if x_c >= r * peak (the top type always is).
Let E(r) = number of exposed types.

LEMMA (upper bound).  E(r) <= 1 + (1/t - 1) / r.
Proof. Total expression S = sum x_c >= peak + (E-1) * r * peak (the top type plus E-1 other exposed types, each at
least r*peak; all other terms are >= 0). Divide by peak: 1/t = S/peak >= 1 + (E-1) r. Rearrange. QED.

LEMMA (lower bound).  t >= 1 / (E + (N - E) r).
Proof. S <= peak * E + (N - E) * r * peak, because exposed types are each <= peak and unexposed each < r*peak. QED.

TIGHTNESS. For given t, r, N the bound on E is attained (up to the integer floor) by putting E-1 types exactly at
r*peak and spreading the remaining mass below r*peak over the other types, provided N is large enough. So if the
only thing known about a gene is t, the largest number of exposed cell types that can be ruled in is exactly this
bound: top_share is the sufficient statistic for the WORST-CASE breadth.

CONSEQUENCE for the threshold form of the theorem. With weights w_c in [0,1] and exposure indicators, the summed
harm H(r) = sum_c w_c * 1[x_c >= r*peak] <= E(r) <= 1 + (1/t - 1)/r, which decreases as t rises. Ordering genes by
top_share is therefore ordering them by a guaranteed ceiling on threshold-form harm. It is one-directional: a high
top_share guarantees limited exposure; a low top_share does not guarantee wide exposure (the second lemma bounds it
from the other side only when E is known).
"""
import json, os, math
import numpy as np
import continuous_weights_v2 as C2
import target_labels as TL

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'breadth_bound_lemma.json')
RS = (0.03, 0.1, 0.3)


def bounds(t, r, N, E):
    return 1 + (1 / t - 1) / r, 1.0 / (E + (N - E) * r)


def worst_case_profile(t, r, N):
    """Construct a profile with top_share t that maximises exposed types at level r; return (E, actual top_share)."""
    peak = 1.0
    S = peak / t
    k = min(int(math.floor((S - peak) / (r * peak) + 1e-12)), N - 1)
    x = [peak] + [r * peak] * k
    rest = S - sum(x)
    others = N - 1 - k
    if others > 0 and rest > 1e-12:
        each = min(rest / others, r * peak * (1 - 1e-9))
        x += [each] * others
    x = np.array(x); x = x / x.sum()
    E = int((x >= r * x.max() - 1e-12).sum())
    return E, float(x.max())


def main():
    shares = C2.shares()
    labels = {g: v for g, v in json.load(open(os.path.join(HERE, 'flag_cleanup.json')))['targets'].items()}
    rows, viol = [], 0
    for g, (types, s) in shares.items():
        s = np.asarray(s, float); N = len(s); t = float(s.max())
        for r in RS:
            E = int((s >= r * s.max() - 1e-12).sum())
            ub, lb = bounds(t, r, N, E)
            ok = (E <= ub + 1e-9) and (t >= lb - 1e-9)
            viol += (not ok)
            rows.append({'gene': g, 'label': labels.get(g, {}).get('label'), 'r': r, 't': t, 'N': N, 'E': E,
                         'E_upper_bound': ub, 't_lower_bound': lb, 'holds': bool(ok), 'slack_ratio': E / ub})
    # tightness check by construction
    rng = np.random.default_rng(3); tight = []
    for _ in range(2000):
        t = float(rng.uniform(0.02, 0.9)); r = float(rng.choice(RS)); N = 154
        E, t2 = worst_case_profile(t, r, N)
        ub = 1 + (1 / t2 - 1) / r
        tight.append((ub - E) < 1.0 + 1e-6 and E <= ub + 1e-9)
    # random-vector violation search
    rv = 0
    for _ in range(20000):
        N = int(rng.integers(2, 200)); x = rng.pareto(1.0, N) * rng.uniform(0.1, 5); r = float(rng.uniform(0.01, 1))
        t = x.max() / x.sum(); E = int((x >= r * x.max()).sum())
        rv += not (E <= 1 + (1 / t - 1) / r + 1e-9 and t >= 1.0 / (E + (N - E) * r) - 1e-12)
    out = {'n_genes': len(shares), 'n_checks_panel': len(rows), 'violations_panel': int(viol),
           'random_vectors_checked': 20000, 'violations_random': int(rv),
           'tightness_constructions': len(tight), 'tightness_all_within_1': bool(all(tight)),
           'median_slack_ratio_by_r': {str(r): float(np.median([x['slack_ratio'] for x in rows if x['r'] == r])) for r in RS},
           'safe_vs_toxic_median_bound_r0.1': {lab: float(np.median([x['E_upper_bound'] for x in rows if x['r'] == 0.1 and x['label'] == lab]))
                                                  for lab in ('safe', 'toxic')},
           'rows': rows}
    from scipy.stats import spearmanr
    for r in RS:
        x = [a for a in rows if a['r'] == r]; pan = [a for a in x if a['label'] in ('safe', 'toxic')]
        out.setdefault('spearman_top_share_vs_E', {})[str(r)] = float(spearmanr([a['t'] for a in x], [a['E'] for a in x])[0])
        out.setdefault('nonvacuous_bound_all', {})[str(r)] = int(sum(a['E_upper_bound'] < a['N'] for a in x))
        out.setdefault('nonvacuous_bound_panel30', {})[str(r)] = int(sum(a['E_upper_bound'] < a['N'] for a in pan))
    json.dump(out, open(OUT, 'w'), indent=1)
    print({k: v for k, v in out.items() if k != 'rows'})
    for g in ('CD20', 'MS4A1', 'CEACAM5', 'ERBB2', 'EGFR', 'CD47', 'CD19'):
        for x in rows:
            if x['gene'] == g and x['r'] == 0.1:
                print(g, x['label'], round(x['t'], 3), 'E', x['E'], 'bound', round(x['E_upper_bound'], 1))


if __name__ == '__main__':
    main()
