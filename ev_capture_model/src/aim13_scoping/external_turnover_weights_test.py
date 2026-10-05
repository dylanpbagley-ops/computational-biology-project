"""Third continuous-expendability attempt: measured cell-turnover weights (Sender & Milo 2021).
Precommitted in external_turnover_precommitment.md (schemes X1-X6, predictions, decision rule) before running.
"""
import json, os
import numpy as np
import external_turnover_weights as X
import continuous_weights_v2 as C2
import wholebody_v2 as WB

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'external_turnover_weights_test.json')
N_SHUF, N_BOOT, SEED = 1000, 5000, 20260930
SCHEMES = ('X1', 'X2', 'X3', 'X4', 'X5', 'X6')


def score_all(S, A, sch, wfun, ffun):
    out = {}
    for g, (types, s) in S.items():
        top = float(s.max())
        home = A[g]['home']
        harm = float(sum(si * wfun(c) for c, si in zip(types, s)))
        out[g] = {'X1': 1 - harm, 'X2': 1 - wfun(home), 'X3': top + 1 - harm,
                  'X4': top + 1 - wfun(home), 'X5': top + ffun(30)(home), 'X6': top + ffun(365)(home)}[sch]
    return out


def compute_all():
    d = json.load(open(os.path.join(HERE, 'wholebody_v2.json')))
    labels, A = d['labels'], d['reference_A_cell_type']
    S = {g: v for g, v in C2.shares().items() if g in A}
    safe = [g for g in S if labels.get(g) == 'safe']
    tox = [g for g in S if labels.get(g) == 'toxic']
    types = sorted({c for t, _ in S.values() for c in t})
    wf = X.weight_fun(types)
    wvec = np.array([wf(c) for c in types])
    fvec = {30: np.array([X.flag_fun(30)(c) for c in types]), 365: np.array([X.flag_fun(365)(c) for c in types])}
    rng = np.random.default_rng(SEED)
    hh = d['haem_home_A']
    binary = {g: A[g]['top_share'] + (1.0 if hh[g] else 0.0) for g in A}
    one = {g: A[g]['top_share'] for g in A}
    base = {'one_factor': float(C2.sep_of(one, labels, safe, tox)[0]),
            'binary_two_factor': float(C2.sep_of(binary, labels, safe, tox)[0])}
    mapped_share = float(np.mean([sum(si for c, si in zip(*S[g]) if c in X.MAP) for g in S]))
    home_mapped = {g: A[g]['home'] in X.MAP for g in safe + tox}
    res = {'baselines': base, 'n_safe': len(safe), 'n_toxic': len(tox), 'n_types': len(types),
           'n_types_mapped': sum(c in X.MAP for c in types), 'neutral_weight': X.neutral_weight(types),
           'mean_expression_share_on_mapped_types': mapped_share,
           'targets_with_mapped_home': int(sum(home_mapped.values())), 'targets_total': len(home_mapped),
           'unmapped_home_targets': sorted(g for g, v in home_mapped.items() if not v), 'schemes': {}}
    for sch in SCHEMES:
        sc = score_all(S, A, sch, wf, X.flag_fun)
        s, a, b = C2.sep_of(sc, labels, safe, tox)
        p, exact = WB.exact_p(a, b)
        nul = np.empty(N_SHUF)
        for i in range(N_SHUF):
            if sch in ('X5', 'X6'):
                v = fvec[30 if sch == 'X5' else 365]; m = dict(zip(types, rng.permutation(v)))
                sci = score_all(S, A, sch, wf, lambda t, m=m: (lambda c: m[c]))
            else:
                m = dict(zip(types, rng.permutation(wvec)))
                sci = score_all(S, A, sch, lambda c, m=m: m[c], X.flag_fun)
            nul[i] = C2.sep_of(sci, labels, safe, tox)[0]
        boots = np.empty(N_BOOT); diffs = np.empty(N_BOOT)
        for i in range(N_BOOT):
            sa = list(rng.choice(safe, len(safe))); sb = list(rng.choice(tox, len(tox)))
            x = WB.sep([sc[g] for g in sa], [sc[g] for g in sb]); y = WB.sep([binary[g] for g in sa], [binary[g] for g in sb])
            boots[i] = x; diffs[i] = x - y
        loo = [C2.sep_of(sc, labels, [x for x in safe if x != g], [x for x in tox if x != g])[0] for g in safe + tox]
        res['schemes'][sch] = {
            'sep': float(s), 'p_label_permutation': float(p), 'p_exact': bool(exact),
            'ci': [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))],
            'diff_vs_binary': float(s - base['binary_two_factor']),
            'diff_ci': [float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))],
            'null_mean': float(nul.mean()), 'null_p95': float(np.percentile(nul, 95)),
            'p_vs_shuffled_weights': float((np.sum(nul >= s - 1e-12) + 1) / (N_SHUF + 1)),
            'loo_min': float(min(loo)), 'loo_max': float(max(loo))}
    b = base['binary_two_factor']
    ok = [k for k, v in res['schemes'].items() if v['sep'] >= b - 1e-12 and v['p_vs_shuffled_weights'] < 0.05]
    res['decision'] = {'rule': 'adopt only if sep >= binary AND shuffled-weight null p < 0.05',
                       'schemes_meeting_rule': ok, 'binary_stays_primary': len(ok) == 0}
    return res


def main():
    r = compute_all()
    print('baselines', {k: round(v, 3) for k, v in r['baselines'].items()})
    print('types mapped %d/%d; expression share on mapped types %.2f; targets with mapped home %d/%d'
          % (r['n_types_mapped'], r['n_types'], r['mean_expression_share_on_mapped_types'], r['targets_with_mapped_home'], r['targets_total']))
    print('%-3s %6s %8s %-14s %7s %-16s %7s %6s' % ('id', 'sep', 'p_perm', '95% CI', 'vs bin', 'paired CI', 'p_null', 'LOOmin'))
    for k, v in r['schemes'].items():
        print('%-3s %6.3f %8.5f [%.3f,%.3f] %+7.3f [%+.3f,%+.3f] %7.3f %6.3f' % (
            k, v['sep'], v['p_label_permutation'], *v['ci'], v['diff_vs_binary'], *v['diff_ci'], v['p_vs_shuffled_weights'], v['loo_min']))
    print('decision', r['decision'])
    json.dump(r, open(OUT, 'w'), indent=1)


if __name__ == '__main__':
    main()
