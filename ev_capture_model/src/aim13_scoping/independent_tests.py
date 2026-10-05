"""Test 1 (naming) and Test 2 (vital-share candidate rules). Precommitted in independent_tests_precommitment.md."""
import json, os
from math import comb
import numpy as np
from scipy.stats import rankdata, wilcoxon
import target_labels as TL
import target_labels_validation as V
import lineage_map as LM
import flag_cleanup as FC
import expendability_weights as W
import external_turnover_weights as X
import damage_tissue_map as DM
import wholebody_v2 as WB
import continuous_weights_v2 as C2

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'independent_tests.json')
WB.GENES.update(t[0] for t in V.VALIDATION)
EXPENDABLE = {'blood', 'blood_other', 'germline', 'placental'}
N_BOOT, SEED = 5000, 20261004


def sign_test_p(better, worse):
    n = better + worse
    if n == 0:
        return 1.0
    k = max(better, worse)
    return min(1.0, 2 * sum(comb(n, i) for i in range(k, n + 1)) / 2 ** n)


def raw_shares():
    vals, _ = WB.load('rna_single_cell_type.tsv.zip', ['nCPM', 'nTPM', 'pTPM'], ['Cell type'])
    out = {}
    for g, d in vals.items():
        ks = list(d.keys()); x = np.array([d[k] for k in ks], float)
        if x.sum() > 0:
            out[g] = (ks, x / x.sum())
    return out


def naming(S):
    types = S[next(iter(S))][0]
    tw = X.weight_fun(types)
    models = {
        'N1': lambda ks, s: s,
        'N2': lambda ks, s: s * np.array([LM.lineage(c) not in EXPENDABLE for c in ks], float),
        'N3': lambda ks, s: s * np.array([W.weight(c, 'continuous') for c in ks]),
        'N4': lambda ks, s: s * np.array([tw(c) for c in ks]),
    }
    res = {'per_target': {}, 'summary': {}}
    for g, dmg in DM.DAMAGED.items():
        ks, s = S[g]
        row = {}
        for m, f in models.items():
            r = rankdata(-f(ks, s), method='average')
            row[m] = float(min(r[ks.index(c)] for c in dmg))
        res['per_target'][g] = row
    for m in models:
        ranks = np.array([res['per_target'][g][m] for g in DM.DAMAGED])
        res['summary'][m] = {'median_best_rank': float(np.median(ranks)), 'mean_log2_rank': float(np.mean(np.log2(ranks))),
                             'hit_at_3': int((ranks <= 3).sum()), 'n': len(ranks)}
    for m in ('N2', 'N3', 'N4'):
        d = [res['per_target'][g][m] - res['per_target'][g]['N1'] for g in DM.DAMAGED]
        better = sum(x < 0 for x in d); worse = sum(x > 0 for x in d)
        res['summary'][m].update({'better_than_N1': better, 'worse_than_N1': worse, 'sign_test_p': sign_test_p(better, worse)})
    return res


def rules(S):
    d = json.load(open(os.path.join(HERE, 'wholebody_v2.json')))
    labels, A, hh = d['labels'], d['reference_A_cell_type'], d['haem_home_A']
    genes = [g for g in A if labels.get(g) in ('safe', 'toxic')]
    types = S[next(iter(S))][0]
    h0 = {c: 1.0 if TL.haematopoietic(c) else 0.0 for c in types}
    sub = {g: S[g] for g in genes}
    fl = FC.flags(sub, h0, LM.lineage)
    def R1(g):
        ks, s = S[g]
        return 1.0 - float((s * np.array([LM.lineage(c) not in EXPENDABLE for c in ks], float)).max())
    sc = {'M2': {g: A[g]['top_share'] + (1.0 if fl[g][1] else 0.0) for g in genes},
          'R1': {g: R1(g) for g in genes}, 'R2': {g: R1(g) + (1.0 if fl[g][1] else 0.0) for g in genes}}
    safe = [g for g in genes if labels[g] == 'safe']; tox = [g for g in genes if labels[g] == 'toxic']
    rng = np.random.default_rng(SEED)
    out = {}
    boots = {k: np.empty(N_BOOT) for k in sc}
    for i in range(N_BOOT):
        sa = list(rng.choice(safe, len(safe))); sb = list(rng.choice(tox, len(tox)))
        for k in sc:
            boots[k][i] = WB.sep([sc[k][g] for g in sa], [sc[k][g] for g in sb])
    for k in sc:
        s_, a_, b_ = C2.sep_of(sc[k], labels, safe, tox); p, _ = WB.exact_p(a_, b_)
        diff = boots[k] - boots['M2']
        out[k] = {'sep': float(s_), 'p': float(p), 'ci': [float(np.percentile(boots[k], 2.5)), float(np.percentile(boots[k], 97.5))],
                  'diff_vs_M2': float(s_ - WB.sep([sc['M2'][g] for g in safe], [sc['M2'][g] for g in tox])),
                  'diff_ci': [float(np.percentile(diff, 2.5)), float(np.percentile(diff, 97.5))]}
    # validation primary set (genes present in the reference)
    vg = [g for g in V.PRIMARY_SAFE + V.PRIMARY_TOXIC if g in S]
    vsub = {g: S[g] for g in vg}; vfl = FC.flags(vsub, h0, LM.lineage)
    vs = {'M2': {g: float(S[g][1].max()) + (1.0 if vfl[g][1] else 0.0) for g in vg},
          'R1': {g: R1(g) for g in vg}, 'R2': {g: R1(g) + (1.0 if vfl[g][1] else 0.0) for g in vg}}
    ps = [g for g in V.PRIMARY_SAFE if g in S]; pt = [g for g in V.PRIMARY_TOXIC if g in S]
    for k in vs:
        out[k]['validation_sep'] = float(WB.sep([vs[k][g] for g in ps], [vs[k][g] for g in pt]))
    return out


def main():
    S = raw_shares()
    nm = naming(S); rl = rules(S)
    a_ok = nm['summary']['N2']['sign_test_p'] < 0.05 and nm['summary']['N2']['better_than_N1'] > nm['summary']['N2']['worse_than_N1']
    dec = {}
    for k in ('R1', 'R2'):
        b_ok = rl[k]['sep'] >= 0.88 and rl[k]['diff_ci'][0] <= 0 <= rl[k]['diff_ci'][1]
        c_ok = rl[k]['validation_sep'] >= rl['M2']['validation_sep'] - 1e-12
        dec[k] = {'a_naming': bool(a_ok), 'b_panel': bool(b_ok), 'c_validation': bool(c_ok), 'adopt': bool(a_ok and b_ok and c_ok)}
    res = {'naming': nm, 'rules': rl, 'decision': dec}
    json.dump(res, open(OUT, 'w'), indent=1)
    print('NAMING best rank of a documented-damaged cell type (1 = perfect)')
    print('%-8s %6s %6s %6s %6s' % ('target', 'N1', 'N2', 'N3', 'N4'))
    for g, r in nm['per_target'].items():
        print('%-8s %6.1f %6.1f %6.1f %6.1f' % (g, r['N1'], r['N2'], r['N3'], r['N4']))
    for m, s in nm['summary'].items():
        print(m, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in s.items()})
    for k, v in rl.items():
        print(k, {a: (round(b, 3) if isinstance(b, float) else [round(z, 3) for z in b] if isinstance(b, list) else b) for a, b in v.items()})
    print('decision', dec)


if __name__ == '__main__':
    main()
