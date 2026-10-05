"""Round 6A: are the expendability weights the bottleneck? Precommitted in weights_diagnostic_precommitment.md."""
import json, os, collections
import numpy as np
import target_labels as TL
import lineage_map as LM
import flag_cleanup as FC
import continuous_weights_v2 as C2

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'weights_diagnostic.json')
SEED = 20260929
N_PERM = 200
FAMS = ['blood', 'blood_other', 'germline', 'placental', 'neural', 'endocrine', 'endothelial', 'stromal_muscle', 'other', 'epithelial']
GRID = np.linspace(0, 1.5, 7)


def auc(sa, sb):
    d = sa[:, None] - sb[None, :]
    return float(((d > 1e-12).sum() + 0.5 * (np.abs(d) <= 1e-12).sum()) / d.size)


def smooth(sa, sb, tau=0.04):
    d = (sa[:, None] - sb[None, :]) / tau
    return float(np.mean(1 / (1 + np.exp(-d))))


def fit(t, X, y, rng, starts=8, sweeps=4):
    """maximise smoothed AUC over a in [0,1.5]^K; y=1 for safe. returns a."""
    K = X.shape[1]
    best, besta = -1, None
    inits = [np.zeros(K)] + [rng.choice(GRID, K) for _ in range(starts)]
    for a in inits:
        a = a.copy()
        cur = smooth((t + X @ a)[y == 1], (t + X @ a)[y == 0])
        for _ in range(sweeps):
            moved = False
            for k in range(K):
                for v in GRID:
                    if v == a[k]:
                        continue
                    b = a.copy(); b[k] = v
                    s = (t + X @ b)
                    val = smooth(s[y == 1], s[y == 0])
                    if val > cur + 1e-9:
                        cur, a, moved = val, b, True
            if not moved:
                break
        if cur > best:
            best, besta = cur, a
    return besta


def main():
    import wholebody_v2 as WB
    d = json.load(open(os.path.join(HERE, 'wholebody_v2.json')))
    labels, A = d['labels'], d['reference_A_cell_type']
    S = {g: v for g, v in C2.shares().items() if g in A}
    safe = [g for g in S if labels.get(g) == 'safe']; tox = [g for g in S if labels.get(g) == 'toxic']
    genes = safe + tox
    y = np.array([1] * len(safe) + [0] * len(tox))
    types = sorted({c for tt, _ in S.values() for c in tt}); LM.check(types)
    h0 = {c: 1.0 if TL.haematopoietic(c) else 0.0 for c in types}
    fl = FC.flags(S, h0, LM.lineage)
    t = np.array([float(S[g][1].max()) for g in genes])
    F = np.zeros((len(genes), len(FAMS)))          # composition shares
    Fh = np.zeros((len(genes), len(FAMS)))         # one-hot largest family
    for i, g in enumerate(genes):
        types_g, s = S[g]
        for c, si in zip(types_g, s):
            F[i, FAMS.index(LM.lineage(c))] += si
        Fh[i, FAMS.index(fl[g][2])] = 1
    frozen = t + Fh[:, 0]
    res = {'n_safe': len(safe), 'n_toxic': len(tox), 'frozen_sep': auc(frozen[y == 1], frozen[y == 0])}
    nb = Fh[:, 0] == 0
    res['n_nonblood_safe'] = int((y[nb] == 1).sum()); res['n_nonblood_toxic'] = int((y[nb] == 0).sum())
    res['blood_flag_only_sep'] = auc(Fh[y == 1, 0], Fh[y == 0, 0])
    res['frozen_within_nonblood_sep'] = auc(t[nb & (y == 1)], t[nb & (y == 0)])
    res['blood_flagged'] = {'safe': int(Fh[y == 1, 0].sum()), 'toxic': int(Fh[y == 0, 0].sum())}
    rng = np.random.default_rng(SEED)
    for cls, X in (('F', Fh), ('G', F)):
        a = fit(t, X, y, rng)
        sc = t + X @ a
        r = {'a': dict(zip(FAMS, map(float, a))), 'oracle_sep': auc(sc[y == 1], sc[y == 0]),
             'oracle_within_nonblood': auc(sc[nb & (y == 1)], sc[nb & (y == 0)])}
        # nested LOTO
        pred = np.zeros(len(genes))
        for j in range(len(genes)):
            m = np.arange(len(genes)) != j
            aj = fit(t[m], X[m], y[m], rng, starts=4, sweeps=3)
            pred[j] = t[j] + X[j] @ aj
        r['nested_cv_sep'] = auc(pred[y == 1], pred[y == 0])
        # permutation null of the oracle
        nul = []
        for _ in range(N_PERM):
            yp = rng.permutation(y)
            ap = fit(t, X, yp, rng, starts=3, sweeps=3)
            sp = t + X @ ap
            nul.append(auc(sp[yp == 1], sp[yp == 0]))
        nul = np.array(nul)
        r['null_oracle_mean'] = float(nul.mean()); r['null_oracle_p95'] = float(np.percentile(nul, 95))
        r['p_oracle_vs_null'] = float((np.sum(nul >= r['oracle_sep'] - 1e-12) + 1) / (N_PERM + 1))
        res[cls] = r
    # which pairs stay misordered under the oracle F
    aF = np.array([res['F']['a'][k] for k in FAMS]); scF = t + Fh @ aF
    bad = [(genes[i], genes[j]) for i in range(len(genes)) if y[i] == 1 for j in range(len(genes)) if y[j] == 0 and scF[i] <= scF[j]]
    res['oracle_F_misordered_pairs'] = len(bad)
    res['oracle_F_misordered_toxic_counts'] = dict(collections.Counter(b for _, b in bad))
    res['oracle_F_misordered_safe_counts'] = dict(collections.Counter(a for a, _ in bad))
    res['largest_family'] = {g: fl[g][2] for g in genes}
    json.dump(res, open(OUT, 'w'), indent=1)
    return res


if __name__ == '__main__':
    r = main()
    for k in ('n_safe', 'n_toxic', 'frozen_sep', 'blood_flag_only_sep', 'blood_flagged', 'n_nonblood_safe', 'n_nonblood_toxic', 'frozen_within_nonblood_sep'):
        print(k, r[k])
    for cls in 'FG':
        x = r[cls]
        print(cls, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in x.items() if k != 'a'})
        print('  a:', {k: round(v, 2) for k, v in x['a'].items()})
    print('oracle F misordered', r['oracle_F_misordered_pairs'], r['oracle_F_misordered_toxic_counts'], r['oracle_F_misordered_safe_counts'])
