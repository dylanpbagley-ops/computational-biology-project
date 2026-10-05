"""Round 5B: threshold-form harm. Precommitted in threshold_harm_precommitment.md."""
import json, os
import numpy as np
from scipy.stats import rankdata
import target_labels as TL
import target_labels_validation as V
import lineage_map as LM
import expendability_weights as W
import rule_improvement as RI
import flag_cleanup as FC
import continuous_weights_v2 as C2
import wholebody_v2 as WB

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'threshold_harm.json')
WB.GENES.update(t[0] for t in V.VALIDATION)
EXPENDABLE = {'blood', 'blood_other', 'germline', 'placental'}
N_PERM, N_BOOT, SEED = 2000, 5000, 20261005


def scheme_scores(vals, genes):
    out = {}
    for t in (1, 10, 50):
        out['TH1(%d)' % t] = {}; out['TH2(%d)' % t] = {}
        for g in genes:
            d = vals[g]; ks = list(d.keys()); x = np.array([d[k] for k in ks], float)
            above = x >= t
            vit = np.array([LM.lineage(k) not in EXPENDABLE for k in ks])
            wv = np.array([W.weight(k, 'continuous') for k in ks])
            out['TH1(%d)' % t][g] = float(-(above & vit).sum())
            out['TH2(%d)' % t][g] = float(-(wv * above).sum())
    return out


def main():
    d = json.load(open(os.path.join(HERE, 'wholebody_v2.json')))
    labels, A, hh = d['labels'], d['reference_A_cell_type'], d['haem_home_A']
    vals, _ = WB.load('rna_single_cell_type.tsv.zip', ['nCPM', 'nTPM', 'pTPM'], ['Cell type'])
    genes = [g for g in A if labels.get(g) in ('safe', 'toxic')]
    sch = scheme_scores(vals, genes)
    sch['frozen rule'] = {g: A[g]['top_share'] + (1.0 if hh[g] else 0.0) for g in genes}
    ids = list(sch); is_safe = np.array([labels[g] == 'safe' for g in genes])
    mat = np.array([[sch[k][g] for g in genes] for k in ids])
    obs = np.array([RI.auc(mat[i], is_safe) for i in range(len(ids))])
    rng = np.random.default_rng(SEED)
    per = np.array([[RI.auc(mat[i], rng.permutation(is_safe)) for i in range(len(ids))] for _ in range(N_PERM)])
    maxes = per.max(axis=1); bi = ids.index('frozen rule')
    si = np.where(is_safe)[0]; ti = np.where(~is_safe)[0]
    boots = np.empty((N_BOOT, len(ids)))
    for b in range(N_BOOT):
        sa = rng.choice(si, len(si)); sb = rng.choice(ti, len(ti))
        for i in range(len(ids)):
            boots[b, i] = WB.sep(mat[i][sa], mat[i][sb])
    # validation primary set
    S = {g: v for g, v in C2.shares().items() if g in V.PRIMARY_SAFE + V.PRIMARY_TOXIC}
    vg = [g for g in V.PRIMARY_SAFE + V.PRIMARY_TOXIC if g in S]
    vs = scheme_scores(vals, vg)
    types = sorted({c for t, _ in S.values() for c in t})
    h0 = {c: 1.0 if TL.haematopoietic(c) else 0.0 for c in types}
    vfl = FC.flags(S, h0, LM.lineage)
    vs['frozen rule'] = {g: float(S[g][1].max()) + (1.0 if vfl[g][1] else 0.0) for g in vg}
    ps = [g for g in V.PRIMARY_SAFE if g in S]; pt = [g for g in V.PRIMARY_TOXIC if g in S]
    res = {'schemes': {}}
    for i, k in enumerate(ids):
        diff = boots[:, i] - boots[:, bi]
        res['schemes'][k] = {'sep': float(obs[i]), 'p_familywise': float((np.sum(maxes >= obs[i] - 1e-12) + 1) / (N_PERM + 1)),
            'ci': [float(np.percentile(boots[:, i], 2.5)), float(np.percentile(boots[:, i], 97.5))],
            'diff_vs_frozen': float(obs[i] - obs[bi]), 'diff_ci': [float(np.percentile(diff, 2.5)), float(np.percentile(diff, 97.5))],
            'validation_sep': float(WB.sep([vs[k][g] for g in ps], [vs[k][g] for g in pt]))}
    res['decision'] = {'promoted': [k for k, v in res['schemes'].items() if k != 'frozen rule' and v['sep'] >= 0.944 and v['p_familywise'] < 0.05],
                       'null_max_p95': float(np.percentile(maxes, 95))}
    json.dump(res, open(OUT, 'w'), indent=1)
    print('null max 95th pct %.3f' % res['decision']['null_max_p95'])
    for k, v in res['schemes'].items():
        print('%-12s sep %.3f p_family %.3f CI [%.3f,%.3f] vs frozen %+.3f [%+.3f,%+.3f] validation %.3f' % (k, v['sep'], v['p_familywise'], *v['ci'], v['diff_vs_frozen'], *v['diff_ci'], v['validation_sep']))
    print(res['decision'])


if __name__ == '__main__':
    main()
