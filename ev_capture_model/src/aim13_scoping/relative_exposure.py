"""Round 6C: relative-exposure harm. Precommitted in relative_exposure_precommitment.md."""
import json, os
import numpy as np
import target_labels as TL
import target_labels_validation as V
import expendability_weights as W
import lineage_map as LM
import flag_cleanup as FC
import continuous_weights_v2 as C2
import wholebody_v2 as WB
import faers_severity as F
import faers_rules as FR
from faers_rules import auc, SYS, SI

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'relative_exposure.json')
RS = (0.03, 0.1, 0.3)
T0 = 10.0
N_PERM, N_BOOT, SEED = 1000, 5000, 20261007
VARS = [(wn, ag) for wn in ('binary', 'hand', 'death', 'severe') for ag in ('max', 'sum')] + [('death_ihc', 'max')]


def rel_exposure(vals, genes, types, t0=T0):
    ti = {c: i for i, c in enumerate(types)}
    X = np.zeros((len(genes), len(types)))
    for gi, g in enumerate(genes):
        for c, x in vals[g].items():
            X[gi, ti[c]] = x
    ref = np.maximum(X.max(1), t0)[:, None]
    E = np.mean([(X >= r * ref) for r in RS], axis=0).astype(float)
    return E


def sys_max(Ec, sysof):
    out = np.zeros((Ec.shape[0], len(SYS)))
    for j, s in enumerate(sysof):
        out[:, s] = np.maximum(out[:, s], Ec[:, j])
    return out


def score(wn, ag, Ec, Es, Eihc, w_cell, ws):
    if wn == 'death_ihc':
        rel = np.where(np.isfinite(Eihc), Eihc / np.maximum(np.nanmax(Eihc, axis=1, keepdims=True), 1e-12), np.nan)
        return -(((Es + rel) / 2) * ws['death']).max(1)
    H = Ec * w_cell
    return -(H.max(1) if ag == 'max' else H.sum(1))


def fit_sys(Es, y, form, rng, starts=4, sweeps=3):
    """diagnostic: nonnegative system-level weights maximising smoothed AUC (in-sample)."""
    grid = np.linspace(0, 1, 5)
    K = Es.shape[1] - 1
    def val(u):
        s = -((Es[:, :K] * u).max(1) if form == 'max' else (Es[:, :K] * u).sum(1))
        a, b = s[y == 1], s[y == 0]
        d = (a[:, None] - b[None, :]) / 0.05
        return float(np.mean(1 / (1 + np.exp(-np.clip(d, -50, 50)))))
    best, bu = -1, None
    for st in range(starts):
        u = rng.choice(grid, K) if st else np.full(K, 0.5)
        cur = val(u)
        for _ in range(sweeps):
            moved = False
            for k in range(K):
                for v in grid:
                    if v != u[k]:
                        b = u.copy(); b[k] = v; c = val(b)
                        if c > cur + 1e-9:
                            cur, u, moved = c, b, True
            if not moved: break
        if cur > best: best, bu = cur, u
    return bu


def main():
    d = json.load(open(os.path.join(HERE, 'wholebody_v2.json')))
    labels, A = d['labels'], d['reference_A_cell_type']
    panel = [g for g in A if labels.get(g) in ('safe', 'toxic')]
    y = np.array([1 if labels[g] == 'safe' else 0 for g in panel])
    vgenes = list(V.PRIMARY_SAFE + V.PRIMARY_TOXIC + V.HEDGED_TOXIC)
    vals, _ = WB.load('rna_single_cell_type.tsv.zip', ['nCPM', 'nTPM', 'pTPM'], ['Cell type'])
    S_all = C2.shares(); S = {g: v for g, v in S_all.items() if g in A or g in vgenes}
    allg = panel + [g for g in vgenes if g in vals and g not in panel]
    npan = len(panel)
    M = FR.build(allg, vals, S)
    types, sysof = M['types'], M['sysof']
    ws = FR.weight_sets()
    hand_w = np.array([W.weight(c, 'continuous') for c in types])
    bin_w = np.array([0.0 if TL.haematopoietic(c) else 1.0 for c in types])
    Ec = rel_exposure(vals, allg, types); Es = sys_max(Ec, sysof)
    def wcell(wn, sw=None):
        if wn == 'binary': return bin_w
        if wn == 'hand': return hand_w
        return (sw if sw is not None else ws[wn.replace('_ihc', '')])[sysof]
    types_S = sorted({c for tt, _ in S.values() for c in tt}); LM.check(types_S)
    h0 = {c: 1.0 if TL.haematopoietic(c) else 0.0 for c in types_S}
    fl = FC.flags(S, h0, LM.lineage)
    frozen = np.array([float(S[g][1].max()) + (1.0 if fl[g][1] else 0.0) if g in S else np.nan for g in allg])
    nonblood = np.array([not fl[g][1] if g in fl else True for g in allg])
    ihc_have = np.isfinite(M['Eihc']).all(1)
    sc = {v: score(v[0], v[1], Ec, Es, M['Eihc'], wcell(v[0]), ws) for v in VARS}
    rng0 = np.random.default_rng(SEED)
    null = {v: np.empty(N_PERM) for v in VARS}
    for i in range(N_PERM):
        pi = np.concatenate([rng0.permutation(len(F.SYSTEM_PTS)), [len(F.SYSTEM_PTS)]])
        pw = {'binary': rng0.permutation(bin_w), 'hand': rng0.permutation(hand_w)}
        wsp = {k: v[pi] for k, v in ws.items()}
        for v in VARS:
            wn, ag = v
            wc = pw[wn] if wn in pw else wsp[wn.replace('_ihc', '')][sysof]
            s = score(wn, ag, Ec[:npan], Es[:npan], M['Eihc'][:npan], wc, wsp)
            m = ihc_have[:npan] if wn == 'death_ihc' else np.ones(npan, bool)
            null[v][i] = FR.evaluate(s, y, m)[0]
    fam = np.max(np.array([null[v] for v in VARS]), axis=0)
    rng = np.random.default_rng(SEED + 1)
    si = np.where(y == 1)[0]; ti_ = np.where(y == 0)[0]
    bidx = [(rng.choice(si, len(si)), rng.choice(ti_, len(ti_))) for _ in range(N_BOOT)]
    res = {'n_panel': npan, 'frozen_sep': FR.evaluate(frozen[:npan], y, np.ones(npan, bool))[0], 'variants': {},
           'null_family_p95': float(np.percentile(fam, 95)), 'n_panel_with_ihc': int(ihc_have[:npan].sum())}
    ix = {g: i for i, g in enumerate(allg)}
    for v in VARS:
        s_all = sc[v]; m0 = ihc_have[:npan] if v[0] == 'death_ihc' else np.ones(npan, bool)
        s, m = FR.evaluate(s_all[:npan], y, m0); fz, _ = FR.evaluate(frozen[:npan], y, m)
        diffs, boots = [], []
        for a, b in bidx:
            a = a[m[a]]; b = b[m[b]]
            x = auc(s_all[a], s_all[b]); z = auc(frozen[a], frozen[b]); boots.append(x); diffs.append(x - z)
        loo = [FR.evaluate(s_all[:npan], y, m & (np.arange(npan) != j))[0] for j in range(npan) if m[j]]
        nb = nonblood[:npan] & m
        vs = [g for g in V.PRIMARY_SAFE if g in ix]; vt = [g for g in V.PRIMARY_TOXIC if g in ix]; vh = [g for g in V.HEDGED_TOXIC if g in ix]
        def vsep(a_, b_):
            a = [s_all[ix[g]] for g in a_ if np.isfinite(s_all[ix[g]]) and g in S]
            b = [s_all[ix[g]] for g in b_ if np.isfinite(s_all[ix[g]]) and g in S]
            return (auc(a, b), len(a), len(b)) if a and b else (None, len(a), len(b))
        res['variants']['%s_%s' % v] = {'sep': s, 'n': int(m.sum()), 'frozen_same_genes': fz, 'diff_vs_frozen': s - fz,
            'diff_ci': [float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))], 'ci': [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))],
            'p_shuffled_weights': float((np.sum(null[v] >= s - 1e-12) + 1) / (N_PERM + 1)), 'p_familywise': float((np.sum(fam >= s - 1e-12) + 1) / (N_PERM + 1)),
            'loo_min': float(min(loo)), 'nonblood_sep': FR.evaluate(s_all[:npan], y, nb)[0], 'nonblood_frozen': FR.evaluate(frozen[:npan], y, nb)[0], 'n_nonblood': int(nb.sum()),
            'validation_primary': vsep(vs, vt), 'validation_safe_vs_hedged': vsep(vs, vh)}
    rob = {}
    for t0 in (3.0, 10.0, 30.0):
        Ec2 = rel_exposure(vals, allg, types, t0)
        for wn in ('death', 'binary'):
            for ag in ('max', 'sum'):
                s2 = score(wn, ag, Ec2, sys_max(Ec2, sysof), M['Eihc'], wcell(wn), ws)
                rob['T0=%g_%s_%s' % (t0, wn, ag)] = FR.evaluate(s2[:npan], y, np.ones(npan, bool))[0]
    res['robustness_T0'] = rob
    # diagnostics: system-level oracle and nested CV
    rng2 = np.random.default_rng(SEED + 2); diag = {}
    E29 = Es[:npan]
    for form in ('sum', 'max'):
        u = fit_sys(E29, y, form, rng2)
        s = -((E29[:, :-1] * u).max(1) if form == 'max' else (E29[:, :-1] * u).sum(1))
        pred = np.zeros(npan)
        for j in range(npan):
            mk = np.arange(npan) != j
            uj = fit_sys(E29[mk], y[mk], form, rng2, starts=3, sweeps=2)
            pred[j] = -((E29[j, :-1] * uj).max() if form == 'max' else (E29[j, :-1] * uj).sum())
        nl = []
        for _ in range(60):
            yp = rng2.permutation(y); up = fit_sys(E29, yp, form, rng2, starts=2, sweeps=2)
            sp = -((E29[:, :-1] * up).max(1) if form == 'max' else (E29[:, :-1] * up).sum(1)); nl.append(auc(sp[yp == 1], sp[yp == 0]))
        diag[form] = {'oracle_sep': auc(s[y == 1], s[y == 0]), 'nested_cv_sep': auc(pred[y == 1], pred[y == 0]),
                      'null_oracle_mean': float(np.mean(nl)), 'null_oracle_p95': float(np.percentile(nl, 95)),
                      'weights': {SYS[i]: float(u[i]) for i in range(len(u))}}
    res['diagnostic_system_weight_oracle'] = diag
    res['decision'] = {'promoted': [k for k, v in res['variants'].items() if v['sep'] >= 0.944 and v['p_familywise'] < 0.05]}
    json.dump(res, open(OUT, 'w'), indent=1)
    return res


if __name__ == '__main__':
    r = main()
    print('frozen', round(r['frozen_sep'], 3), 'family null p95', round(r['null_family_p95'], 3), 'ihc n', r['n_panel_with_ihc'])
    for k, v in r['variants'].items():
        vp = v['validation_primary']; vh = v['validation_safe_vs_hedged']
        print('%-14s sep %.3f (n=%d, frozen %.3f) diff %+.3f [%+.3f,%+.3f] shufp %.3f famp %.3f LOOmin %.3f nonblood %.3f (frozen %.3f n=%d) val %s hedged %s' % (
            k, v['sep'], v['n'], v['frozen_same_genes'], v['diff_vs_frozen'], *v['diff_ci'], v['p_shuffled_weights'], v['p_familywise'], v['loo_min'],
            v['nonblood_sep'], v['nonblood_frozen'], v['n_nonblood'], None if vp[0] is None else (round(vp[0], 3), vp[1], vp[2]), None if vh[0] is None else round(vh[0], 3)))
    print('robustness', {k: round(v, 3) for k, v in r['robustness_T0'].items()})
    for f, x in r['diagnostic_system_weight_oracle'].items():
        print('oracle', f, {k: round(v, 3) for k, v in x.items() if k != 'weights'})
    print('promoted', r['decision']['promoted'])
