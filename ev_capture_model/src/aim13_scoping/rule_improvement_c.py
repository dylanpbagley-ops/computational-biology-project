"""Round 4C: vital-mass and absolute-threshold measures. Precommitted in rule_improvement_c_precommitment.md."""
import json, os
import numpy as np
from scipy.stats import rankdata
import target_labels as TL
import lineage_map as LM
import continuous_weights_v2 as C2
import flag_cleanup as FC
import rule_improvement as RI
import wholebody_v2 as WB

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'rule_improvement_c.json')
N_PERM, N_BOOT, SEED = 2000, 5000, 20261003
EXPENDABLE = {'blood', 'blood_other', 'germline', 'placental'}


def measures_c(vals, genes):
    m = {k: {} for k in ('C0', 'V1', 'V2', 'V3', 'V4', 'D1', 'D2', 'D3')}
    for g in genes:
        d = vals[g]; ks = list(d.keys()); x = np.array([d[k] for k in ks], float)
        s = x / x.sum()
        vital = np.array([LM.lineage(k) not in EXPENDABLE for k in ks])
        m['C0'][g] = float(s.max())
        m['V1'][g] = float(-np.log10(1 + x[vital].sum()))
        m['V2'][g] = float(-np.log10(1 + (x[vital].max() if vital.any() else 0.0)))
        m['V3'][g] = float(s[~vital].sum())
        m['D1'][g] = float(-(x >= 1).sum()); m['D2'][g] = float(-(x >= 10).sum()); m['D3'][g] = float(-(x >= 50).sum())
    gs = list(genes)
    r0 = rankdata([m['C0'][g] for g in gs]); r1 = rankdata([m['V1'][g] for g in gs])
    for g, a, b in zip(gs, r0, r1):
        m['V4'][g] = float((a + b) / 2)
    return m


def compute_all():
    d = json.load(open(os.path.join(HERE, 'wholebody_v2.json')))
    labels = d['labels']; vals = RI.raw_expression()
    S = {g: v for g, v in C2.shares().items() if g in d['reference_A_cell_type']}
    genes = [g for g in S if labels.get(g) in ('safe', 'toxic')]
    types = sorted({c for t, _ in S.values() for c in t})
    h0 = {c: 1.0 if TL.haematopoietic(c) else 0.0 for c in types}
    fl = FC.flags(S, h0, LM.lineage); flag = {g: 1.0 if fl[g][1] else 0.0 for g in genes}
    m = measures_c(vals, genes)
    sch = {}
    for k in ('C0', 'V1', 'V2', 'V3', 'V4', 'D1', 'D2', 'D3'):
        sch[k] = {g: m[k][g] for g in genes}
        sch[k + '+flag'] = {g: 1000.0 * flag[g] + m[k][g] for g in genes}
    ids = list(sch); is_safe = np.array([labels[g] == 'safe' for g in genes])
    mat = np.array([[sch[k][g] for g in genes] for k in ids])
    obs = np.array([RI.auc(mat[i], is_safe) for i in range(len(ids))])
    rng = np.random.default_rng(SEED)
    per = np.empty((N_PERM, len(ids)))
    for p in range(N_PERM):
        lab = rng.permutation(is_safe); per[p] = [RI.auc(mat[i], lab) for i in range(len(ids))]
    maxes = per.max(axis=1); bi = ids.index('C0+flag')
    si = np.where(is_safe)[0]; ti = np.where(~is_safe)[0]
    boots = np.empty((N_BOOT, len(ids)))
    for b in range(N_BOOT):
        sa = rng.choice(si, len(si)); sb = rng.choice(ti, len(ti))
        for i in range(len(ids)):
            boots[b, i] = WB.sep(mat[i][sa], mat[i][sb])
    res = {'schemes': {}}
    for i, k in enumerate(ids):
        loo = [RI.auc(mat[i][np.arange(len(genes)) != j], is_safe[np.arange(len(genes)) != j]) for j in range(len(genes))]
        diff = boots[:, i] - boots[:, bi]
        res['schemes'][k] = {'sep': float(obs[i]),
            'p_label_permutation': float((np.sum(per[:, i] >= obs[i] - 1e-12) + 1) / (N_PERM + 1)),
            'p_familywise': float((np.sum(maxes >= obs[i] - 1e-12) + 1) / (N_PERM + 1)),
            'ci': [float(np.percentile(boots[:, i], 2.5)), float(np.percentile(boots[:, i], 97.5))],
            'diff_vs_c0_flag': float(obs[i] - obs[bi]),
            'diff_ci': [float(np.percentile(diff, 2.5)), float(np.percentile(diff, 97.5))], 'loo_min': float(min(loo))}
    prom = [k for k, v in res['schemes'].items() if v['sep'] >= 0.944 and v['p_familywise'] < 0.05]
    res['decision'] = {'promoted': prom, 'null_max_p95': float(np.percentile(maxes, 95))}
    return res


def main():
    r = compute_all()
    print('null max-separation 95th pct %.3f' % r['decision']['null_max_p95'])
    for k, v in r['schemes'].items():
        print('%-9s %6.3f p_perm %.4f p_family %.3f CI [%.3f,%.3f] vs base %+.3f [%+.3f,%+.3f] LOOmin %.3f' % (
            k, v['sep'], v['p_label_permutation'], v['p_familywise'], *v['ci'], v['diff_vs_c0_flag'], *v['diff_ci'], v['loo_min']))
    print('decision', r['decision'])
    json.dump(r, open(OUT, 'w'), indent=1)


if __name__ == '__main__':
    main()
