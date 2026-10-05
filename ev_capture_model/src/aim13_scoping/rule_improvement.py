"""Round 4B: family of confinement / abundance measures. Precommitted in rule_improvement_precommitment.md."""
import json, os, collections
import numpy as np
from scipy.stats import rankdata
import target_labels as TL
import lineage_map as LM
import continuous_weights_v2 as C2
import flag_cleanup as FC
import wholebody_v2 as WB

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'rule_improvement.json')
N_PERM, N_BOOT, SEED = 2000, 5000, 20261002
MEASURES = ('C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8')


def raw_expression():
    vals, _ = WB.load('rna_single_cell_type.tsv.zip', ['nCPM', 'nTPM', 'pTPM'], ['Cell type'])
    return vals


def measures(vals, genes):
    m = {k: {} for k in MEASURES}
    for g in genes:
        d = vals[g]; ks = list(d.keys()); x = np.array([d[k] for k in ks], float)
        tot = x.sum(); s = x / tot; n = len(x)
        agg = collections.defaultdict(float)
        for c, si in zip(ks, s):
            agg[LM.lineage(c)] += si
        nz = s[s > 0]
        m['C0'][g] = float(s.max())
        m['C1'][g] = float(max(agg.values()))
        m['C2'][g] = float(1 - (-(nz * np.log(nz)).sum()) / np.log(n))
        m['C3'][g] = float(((1 - x / x.max()).sum()) / (n - 1))
        m['C4'][g] = float(np.sort(s)[-3:].sum())
        m['C5'][g] = float(-(x > 0.10 * x.max()).sum())
        m['C6'][g] = float(-np.log10(1 + tot))
        m['C7'][g] = float(-np.log10(1 + x.max()))
    gs = list(genes)
    r0 = rankdata([m['C0'][g] for g in gs]); r6 = rankdata([m['C6'][g] for g in gs])
    for g, a, b in zip(gs, r0, r6):
        m['C8'][g] = float((a + b) / 2)
    return m


def build_schemes(vals, genes, flag):
    m = measures(vals, genes)
    sch = {}
    for k in MEASURES:
        sch[k] = {g: m[k][g] for g in genes}
        sch[k + '+flag'] = {g: 1000.0 * flag[g] + m[k][g] for g in genes}
    return sch


def auc(scores, is_safe):
    r = rankdata(scores); ns = is_safe.sum(); nt = (~is_safe).sum()
    return (r[is_safe].sum() - ns * (ns + 1) / 2) / (ns * nt)


def compute_all():
    d = json.load(open(os.path.join(HERE, 'wholebody_v2.json')))
    labels = d['labels']
    vals = raw_expression()
    S = {g: v for g, v in C2.shares().items() if g in d['reference_A_cell_type']}
    genes = [g for g in S if labels.get(g) in ('safe', 'toxic')]
    types = sorted({c for t, _ in S.values() for c in t})
    h0 = {c: 1.0 if TL.haematopoietic(c) else 0.0 for c in types}
    fl = FC.flags(S, h0, LM.lineage)
    flag = {g: 1.0 if fl[g][1] else 0.0 for g in genes}
    sch = build_schemes(vals, genes, flag)
    ids = list(sch.keys())
    is_safe = np.array([labels[g] == 'safe' for g in genes])
    mat = np.array([[sch[k][g] for g in genes] for k in ids])          # schemes x genes
    obs = np.array([auc(mat[i], is_safe) for i in range(len(ids))])
    rng = np.random.default_rng(SEED)
    maxes = np.empty(N_PERM); per = np.empty((N_PERM, len(ids)))
    for p in range(N_PERM):
        lab = rng.permutation(is_safe)
        per[p] = [auc(mat[i], lab) for i in range(len(ids))]
        maxes[p] = per[p].max()
    base = 'C0+flag'; bi = ids.index(base)
    safe_idx = np.where(is_safe)[0]; tox_idx = np.where(~is_safe)[0]
    boots = np.empty((N_BOOT, len(ids)))
    for b in range(N_BOOT):
        sa = rng.choice(safe_idx, len(safe_idx)); sb = rng.choice(tox_idx, len(tox_idx))
        for i in range(len(ids)):
            boots[b, i] = WB.sep(mat[i][sa], mat[i][sb])
    res = {'n_safe': int(is_safe.sum()), 'n_toxic': int((~is_safe).sum()), 'schemes': {}}
    for i, k in enumerate(ids):
        loo = []
        for j in range(len(genes)):
            keep = np.array([q != j for q in range(len(genes))])
            loo.append(auc(mat[i][keep], is_safe[keep]))
        diff = boots[:, i] - boots[:, bi]
        res['schemes'][k] = {'sep': float(obs[i]),
            'p_label_permutation': float((np.sum(per[:, i] >= obs[i] - 1e-12) + 1) / (N_PERM + 1)),
            'p_familywise': float((np.sum(maxes >= obs[i] - 1e-12) + 1) / (N_PERM + 1)),
            'ci': [float(np.percentile(boots[:, i], 2.5)), float(np.percentile(boots[:, i], 97.5))],
            'diff_vs_c0_flag': float(obs[i] - obs[bi]),
            'diff_ci': [float(np.percentile(diff, 2.5)), float(np.percentile(diff, 97.5))],
            'loo_min': float(min(loo))}
    prom = [k for k, v in res['schemes'].items() if v['sep'] >= 0.944 and v['p_familywise'] < 0.05]
    res['decision'] = {'rule': 'promote only if sep >= 0.944 and familywise p < 0.05; must then be confirmed on new targets',
                       'promoted': prom, 'null_max_mean': float(maxes.mean()), 'null_max_p95': float(np.percentile(maxes, 95))}
    return res


def main():
    r = compute_all()
    print('null max-separation across 18 schemes: mean %.3f, 95th pct %.3f' % (r['decision']['null_max_mean'], r['decision']['null_max_p95']))
    print('%-9s %6s %8s %8s %-14s %7s %-16s %6s' % ('scheme', 'sep', 'p_perm', 'p_family', '95% CI', 'vs base', 'paired CI', 'LOOmin'))
    for k, v in r['schemes'].items():
        print('%-9s %6.3f %8.4f %8.3f [%.3f,%.3f] %+7.3f [%+.3f,%+.3f] %6.3f' % (k, v['sep'], v['p_label_permutation'], v['p_familywise'], *v['ci'], v['diff_vs_c0_flag'], *v['diff_ci'], v['loo_min']))
    print('decision', r['decision'])
    json.dump(r, open(OUT, 'w'), indent=1)


if __name__ == '__main__':
    main()
