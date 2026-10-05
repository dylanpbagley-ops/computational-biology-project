"""Second attempt at continuous expendability: use the WHOLE expression distribution, not just the home cell type.

Precommitted in continuous_weights_v2_precommitment.md (schemes E1-E4, predictions, decision rule) before running.
The theorem (section 29.3) says harm = sum_c w_c * s_c, where s_c is the share of a gene's normal expression in
cell type c. The first continuous scheme used only the top cell type's weight, which is arbitrary for broadly
expressed genes (EGFR -> placental tissue). These schemes use every s_c.

  E1  1 - sum s_c w_c                    hand weights, no confinement term
  E2  sum s_c h_c                        fraction of expression in blood lineages
  E3  top_share + (1 - sum s_c w_c)      hand weights, additive with confinement
  E4  top_share + sum s_c h_c            smooth version of the section 23 binary rule
"""
import json, os
import numpy as np
import expendability_weights as W
import target_labels as TL
import wholebody_v2 as WB

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'continuous_weights_v2.json')
N_SHUF = 1000
N_BOOT = 5000
SEED = 20260929
SCHEMES = ('E1', 'E2', 'E3', 'E4')


def shares():
    vals, _ = WB.load('rna_single_cell_type.tsv.zip', ['nCPM', 'nTPM', 'pTPM'], ['Cell type'])
    S = {}
    for g, d in vals.items():
        v = np.array(list(d.values()), float)
        if v.sum() <= 0:
            continue
        S[g] = (list(d.keys()), v / v.sum())
    return S


def scores(S, scheme, wfun=None, hfun=None):
    wfun = wfun or (lambda c: W.weight(c, 'continuous'))
    hfun = hfun or (lambda c: 1.0 if TL.haematopoietic(c) else 0.0)
    out = {}
    for g, (types, s) in S.items():
        top = float(s.max())
        harm = float(sum(si * wfun(c) for c, si in zip(types, s)))
        haem = float(sum(si * hfun(c) for c, si in zip(types, s)))
        out[g] = {'E1': 1 - harm, 'E2': haem, 'E3': top + 1 - harm, 'E4': top + haem}[scheme]
    return out


def sep_of(sc, labels, safe=None, tox=None):
    a = [sc[g] for g in (safe or [g for g in sc if labels.get(g) == 'safe'])]
    b = [sc[g] for g in (tox or [g for g in sc if labels.get(g) == 'toxic'])]
    return WB.sep(a, b), a, b


def compute_all():
    d = json.load(open(os.path.join(HERE, 'wholebody_v2.json')))
    labels = d['labels']
    A = d['reference_A_cell_type']
    S = shares()
    S = {g: v for g, v in S.items() if g in A}
    safe = [g for g in S if labels.get(g) == 'safe']
    tox = [g for g in S if labels.get(g) == 'toxic']
    rng = np.random.default_rng(SEED)
    types = sorted({c for t, _ in S.values() for c in t})
    hv = {c: (1.0 if TL.haematopoietic(c) else 0.0) for c in types}
    wv = np.array([W.weight(c, 'continuous') for c in types])
    hvec = np.array([hv[c] for c in types])

    # baselines, recomputed from the frozen JSON
    base = {}
    hh = d['haem_home_A']
    one = {g: A[g]['top_share'] for g in A}
    two = {g: A[g]['top_share'] + (1.0 if hh[g] else 0.0) for g in A}
    old = {g: W.score(A[g]['top_share'], A[g]['home'], 'continuous') for g in A}
    for k, sc in (('one_factor', one), ('binary_two_factor', two), ('old_continuous_home', old)):
        base[k] = float(sep_of(sc, labels, safe, tox)[0])

    res = {'baselines': base, 'n_safe': len(safe), 'n_toxic': len(tox), 'schemes': {}}
    bin_sc = two
    for sch in SCHEMES:
        sc = scores(S, sch)
        s, a, b = sep_of(sc, labels, safe, tox)
        p, exact = WB.exact_p(a, b)
        # coverage of hand weights (expression share on types with an explicit weight)
        cov = float(np.mean([sum(si for c, si in zip(*S[g]) if c in W.WEIGHTS) for g in S]))
        # shuffled-weight null
        nul = np.empty(N_SHUF)
        for i in range(N_SHUF):
            if sch in ('E1', 'E3'):
                m = dict(zip(types, rng.permutation(wv)))
                sci = scores(S, sch, wfun=lambda c: m[c])
            else:
                m = dict(zip(types, rng.permutation(hvec)))
                sci = scores(S, sch, hfun=lambda c: m[c])
            nul[i] = sep_of(sci, labels, safe, tox)[0]
        # bootstrap CI and paired difference vs binary
        boots = np.empty(N_BOOT); diffs = np.empty(N_BOOT)
        for i in range(N_BOOT):
            sa = list(rng.choice(safe, len(safe))); sb = list(rng.choice(tox, len(tox)))
            x = WB.sep([sc[g] for g in sa], [sc[g] for g in sb])
            y = WB.sep([bin_sc[g] for g in sa], [bin_sc[g] for g in sb])
            boots[i] = x; diffs[i] = x - y
        # leave-one-target-out
        loo = []
        for g in safe + tox:
            s2 = [x for x in safe if x != g]; t2 = [x for x in tox if x != g]
            loo.append(sep_of(sc, labels, s2, t2)[0])
        res['schemes'][sch] = {
            'sep': float(s), 'p_label_permutation': float(p), 'p_exact': bool(exact),
            'ci': [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))],
            'diff_vs_binary': float(s - base['binary_two_factor']),
            'diff_ci': [float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))],
            'null_mean': float(nul.mean()), 'null_p95': float(np.percentile(nul, 95)),
            'p_vs_shuffled_weights': float((np.sum(nul >= s - 1e-12) + 1) / (N_SHUF + 1)),
            'loo_min': float(min(loo)), 'loo_max': float(max(loo)),
            'hand_weight_coverage': cov,
        }
    b = base['binary_two_factor']
    ok = [k for k, v in res['schemes'].items() if v['sep'] >= b - 1e-12 and v['p_vs_shuffled_weights'] < 0.05]
    res['decision'] = {'rule': 'adopt only if sep >= binary AND shuffled-weight null p < 0.05',
                       'schemes_meeting_rule': ok, 'binary_stays_primary': len(ok) == 0}
    return res


def main():
    r = compute_all()
    print('baselines', {k: round(v, 3) for k, v in r['baselines'].items()})
    print('%-4s %6s %8s %-16s %8s %-17s %8s %7s %6s' % ('id', 'sep', 'p_perm', '95% CI', 'vs bin', 'paired CI', 'p_null', 'LOOmin', 'cov'))
    for k, v in r['schemes'].items():
        print('%-4s %6.3f %8.5f [%.3f,%.3f] %+8.3f [%+.3f,%+.3f] %8.3f %7.3f %6.2f' % (
            k, v['sep'], v['p_label_permutation'], *v['ci'], v['diff_vs_binary'], *v['diff_ci'],
            v['p_vs_shuffled_weights'], v['loo_min'], v['hand_weight_coverage']))
    print('decision', r['decision'])
    json.dump(r, open(OUT, 'w'), indent=1)


if __name__ == '__main__':
    main()
