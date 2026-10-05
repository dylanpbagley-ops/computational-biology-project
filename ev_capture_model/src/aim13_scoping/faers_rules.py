"""Round 6B: measured-severity weights and multi-mechanism rules A1-A7. Precommitted in faers_precommitment.md (+ addendum)."""
import json, os, io, zipfile, collections
import numpy as np
import target_labels as TL
import target_labels_validation as V
import expendability_weights as W
import lineage_map as LM
import flag_cleanup as FC
import continuous_weights_v2 as C2
import wholebody_v2 as WB
import faers_severity as F

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'faers_rules.json')
WB.GENES.update(t[0] for t in V.VALIDATION)
TS = (1, 3, 10, 30, 100)
LEVEL = {'not detected': 0.0, 'low': 1 / 3, 'medium': 2 / 3, 'high': 1.0}
N_PERM, N_BOOT, SEED = 1000, 5000, 20261006
SYS = list(F.SYSTEM_PTS) + ['unmapped']
SI = {s: i for i, s in enumerate(SYS)}
VARIANTS = [(r, w) for r in ('A1', 'A2', 'A3', 'A4', 'A6', 'A7') for w in ('death', 'severe')] + [('A5', 'hand')]


def auc(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = a[:, None] - b[None, :]
    return float(((d > 1e-12).sum() + 0.5 * (np.abs(d) <= 1e-12).sum()) / d.size)


def weight_sets():
    c = json.load(open(F.COUNTS))['systems']
    out = {}
    for name, key in (('death', 'deaths'), ('severe', 'severe')):
        f = {s: c[s][key] / c[s]['reports'] for s in F.SYSTEM_PTS}
        mx = max(f.values())
        w = np.array([f[s] / mx for s in F.SYSTEM_PTS] + [float(np.median([f[s] / mx for s in F.SYSTEM_PTS]))])
        out[name] = w
    return out


def load_ihc(genes):
    z = zipfile.ZipFile(os.path.join(WB.D, 'normal_ihc_data.tsv.zip'))
    E = collections.defaultdict(lambda: np.zeros(len(SYS))); have = set()
    with z.open(z.namelist()[0]) as f:
        fh = io.TextIOWrapper(f, 'utf-8'); fh.readline()
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if p[1] not in genes or p[6] == 'Uncertain':
                continue
            lv = LEVEL.get(p[5].strip().lower())
            if lv is None:
                continue
            have.add(p[1])
            s = F.ihc_system(p[2], p[4]) or 'unmapped'
            i = SI[s]
            E[p[1]][i] = max(E[p[1]][i], lv)
    return E, have


def build(genes, vals, S):
    """matrices over `genes`: RNA exposure by system, IHC exposure by system, share by system, home one-hot, and cell-level exposure."""
    types = sorted({c for g in genes for c in vals[g]})
    ti = {c: i for i, c in enumerate(types)}
    sysof = np.array([SI[F.CELL_SYSTEM.get(c, 'unmapped')] for c in types])
    Erna = np.zeros((len(genes), len(SYS))); Ec = np.zeros((len(genes), len(types)))
    Sh = np.zeros((len(genes), len(SYS))); Hm = np.zeros((len(genes), len(SYS))); top = np.full(len(genes), np.nan)
    for gi, g in enumerate(genes):
        for c, x in vals[g].items():
            e = float(np.mean([x >= t for t in TS]))
            Ec[gi, ti[c]] = e
            Erna[gi, sysof[ti[c]]] = max(Erna[gi, sysof[ti[c]]], e)
        if g in S:
            tps, s = S[g]
            for c, si in zip(tps, s):
                Sh[gi, SI[F.CELL_SYSTEM.get(c, 'unmapped')]] += si
            Hm[gi, SI[F.CELL_SYSTEM.get(tps[int(np.argmax(s))], 'unmapped')]] = 1
            top[gi] = float(s.max())
    ihc, have = load_ihc(set(genes))
    Eihc = np.array([ihc[g] if g in have else np.full(len(SYS), np.nan) for g in genes])
    return dict(types=types, sysof=sysof, Erna=Erna, Ec=Ec, Sh=Sh, Hm=Hm, top=top, Eihc=Eihc)


def score(rule, w, M, hand=None):
    """w: weight vector over SYS (or over cell types for A5). returns score per gene (higher = predicted safer)."""
    if rule == 'A1': return M['top'] + 1 - (M['Hm'] * w).sum(1)
    if rule == 'A2': return M['top'] + 1 - (M['Sh'] * w).sum(1)
    if rule == 'A3': return -(M['Erna'] * w).max(1)
    if rule == 'A4': return -(M['Erna'] * w).sum(1)
    if rule == 'A6': return -(M['Eihc'] * w).max(1)
    if rule == 'A7': return -(((M['Erna'] + M['Eihc']) / 2) * w).max(1)
    if rule == 'A5': return -(M['Ec'] * w).max(1)
    raise ValueError(rule)


def evaluate(sc, y, mask):
    m = mask & np.isfinite(sc)
    return auc(sc[m & (y == 1)], sc[m & (y == 0)]), m


def main():
    d = json.load(open(os.path.join(HERE, 'wholebody_v2.json')))
    labels, A, hh = d['labels'], d['reference_A_cell_type'], d['haem_home_A']
    panel = [g for g in A if labels.get(g) in ('safe', 'toxic')]
    y = np.array([1 if labels[g] == 'safe' else 0 for g in panel])
    vgenes = [g for g in V.PRIMARY_SAFE + V.PRIMARY_TOXIC + V.HEDGED_TOXIC]
    vals, _ = WB.load('rna_single_cell_type.tsv.zip', ['nCPM', 'nTPM', 'pTPM'], ['Cell type'])
    S_all = C2.shares()
    S = {g: v for g, v in S_all.items() if g in A or g in vgenes}
    allg = panel + [g for g in vgenes if g in vals and g not in panel]
    types_all = sorted({c for g in allg for c in vals[g]})
    F.check_cell_types(sorted({c for t, _ in S_all.values() for c in t}))
    M = build(allg, vals, S)
    npan = len(panel)
    ws = weight_sets()
    hand_w = np.array([W.weight(c, 'continuous') for c in M['types']])
    # frozen rule on panel and validation genes
    types_S = sorted({c for tt, _ in S.values() for c in tt}); LM.check(types_S)
    h0 = {c: 1.0 if TL.haematopoietic(c) else 0.0 for c in types_S}
    fl = FC.flags(S, h0, LM.lineage)
    frozen = np.array([float(S[g][1].max()) + (1.0 if fl[g][1] else 0.0) if g in S else np.nan for g in allg])
    nonblood = np.array([not fl[g][1] if g in fl else True for g in allg])
    res = {'n_panel': npan, 'panel_genes': panel, 'variants': {}}
    scores = {}
    for rule, wn in VARIANTS:
        w = hand_w if rule == 'A5' else ws[wn]
        scores[(rule, wn)] = score(rule, w, M)
    ihc_have = np.isfinite(M['Eihc']).all(1)
    res['n_panel_with_ihc'] = int(ihc_have[:npan].sum())
    sym = np.random.default_rng(SEED)
    # shuffled-weight null (same permutation applied to both measured weightings; A5 permutes hand weights over cell types)
    null = {v: np.empty(N_PERM) for v in VARIANTS}
    for i in range(N_PERM):
        pi = np.concatenate([sym.permutation(len(F.SYSTEM_PTS)), [len(F.SYSTEM_PTS)]])
        ph = sym.permutation(hand_w)
        for rule, wn in VARIANTS:
            w = ph if rule == 'A5' else ws[wn][pi]
            sc = score(rule, w, M)[:npan]
            null[(rule, wn)][i] = evaluate(sc, y, np.ones(npan, bool) if rule not in ('A6', 'A7') else ihc_have[:npan])[0]
    fam = np.max(np.array([null[v] for v in VARIANTS]), axis=0)
    rng = np.random.default_rng(SEED + 1)
    si = np.where(y == 1)[0]; ti = np.where(y == 0)[0]
    bidx = [(rng.choice(si, len(si)), rng.choice(ti, len(ti))) for _ in range(N_BOOT)]
    for rule, wn in VARIANTS:
        sc = scores[(rule, wn)]
        need_ihc = rule in ('A6', 'A7')
        base_mask = ihc_have[:npan] if need_ihc else np.ones(npan, bool)
        s, m = evaluate(sc[:npan], y, base_mask)
        fz, _ = evaluate(frozen[:npan], y, m)
        p_shuf = float((np.sum(null[(rule, wn)] >= s - 1e-12) + 1) / (N_PERM + 1))
        p_fam = float((np.sum(fam >= s - 1e-12) + 1) / (N_PERM + 1))
        diffs = []; boots = []
        for sa, sb in bidx:
            sa = sa[m[sa]]; sb = sb[m[sb]]
            x = auc(sc[sa], sc[sb]); z = auc(frozen[sa], frozen[sb]); boots.append(x); diffs.append(x - z)
        loo = [evaluate(sc[:npan], y, m & (np.arange(npan) != j))[0] for j in range(npan) if m[j]]
        nb = nonblood[:npan] & m
        sn, _ = evaluate(sc[:npan], y, nb)
        fn, _ = evaluate(frozen[:npan], y, nb)
        # validation
        vs = [g for g in V.PRIMARY_SAFE if g in allg]; vt = [g for g in V.PRIMARY_TOXIC if g in allg]
        vh = [g for g in V.HEDGED_TOXIC if g in allg]
        ix = {g: i for i, g in enumerate(allg)}
        def vsep(sset, tset, need_S):
            a = [sc[ix[g]] for g in sset if np.isfinite(sc[ix[g]]) and (g in S or not need_S)]
            b = [sc[ix[g]] for g in tset if np.isfinite(sc[ix[g]]) and (g in S or not need_S)]
            return (auc(a, b), len(a), len(b)) if a and b else (None, len(a), len(b))
        res['variants']['%s_%s' % (rule, wn)] = {
            'sep': s, 'n': int(m.sum()), 'frozen_same_genes': fz, 'diff_vs_frozen': s - fz,
            'diff_ci': [float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))],
            'ci': [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))],
            'p_shuffled_weights': p_shuf, 'p_familywise': p_fam, 'loo_min': float(min(loo)),
            'nonblood_sep': sn, 'nonblood_frozen': fn, 'n_nonblood': int(nb.sum()),
            'validation_primary_S_only': vsep(vs, vt, True), 'validation_primary_all_scorable': vsep(vs, vt, False),
            'validation_safe_vs_hedged_toxic': vsep(vs, vh, False)}
    # robustness: A3 with rank weights, A3 with unmapped = min / max
    rob = {}
    for wn in ('death', 'severe'):
        w0 = ws[wn]; ranks = np.concatenate([np.argsort(np.argsort(w0[:-1])) + 1, [0]]) / 28.0
        ranks[-1] = np.median(ranks[:-1])
        for tag, w in (('rank', ranks), ('unmapped_min', np.concatenate([w0[:-1], [w0[:-1].min()]])), ('unmapped_max', np.concatenate([w0[:-1], [w0[:-1].max()]]))):
            rob['A3_%s_%s' % (wn, tag)] = evaluate(score('A3', w, M)[:npan], y, np.ones(npan, bool))[0]
    res['robustness_A3'] = rob
    res['frozen_sep'] = evaluate(frozen[:npan], y, np.ones(npan, bool))[0]
    res['null_family_p95'] = float(np.percentile(fam, 95))
    res['weights'] = {k: {s: float(v[i]) for i, s in enumerate(SYS)} for k, v in ws.items()}
    # per-target harm table for inspection (A3 death, A7 death)
    tab = {}
    for j, g in enumerate(panel):
        h3 = (M['Erna'][j] * ws['death'])
        k = int(np.argmax(h3))
        tab[g] = {'label': labels[g], 'A3_death_H': float(h3.max()), 'A3_worst_system': SYS[k],
                  'A7_death_H': float(scores[('A7', 'death')][j]) if ihc_have[j] else None}
    res['per_target'] = tab
    promoted = [k for k, v in res['variants'].items() if v['sep'] >= 0.944 and v['p_familywise'] < 0.05]
    res['decision'] = {'promoted': promoted, 'rule': 'sep >= 0.944 and family-wise p < 0.05 over all 13 variants'}
    json.dump(res, open(OUT, 'w'), indent=1)
    return res


if __name__ == '__main__':
    r = main()
    print('n panel', r['n_panel'], 'with IHC', r['n_panel_with_ihc'], 'frozen', round(r['frozen_sep'], 3), 'null fam p95', round(r['null_family_p95'], 3))
    for k, v in r['variants'].items():
        print('%-10s sep %.3f (n=%d frozen %.3f) diff %+.3f [%+.3f,%+.3f] shufp %.3f famp %.3f LOOmin %.3f nonblood %.3f (frozen %.3f, n=%d) val %s val_all %s hedged %s' % (
            k, v['sep'], v['n'], v['frozen_same_genes'], v['diff_vs_frozen'], *v['diff_ci'], v['p_shuffled_weights'], v['p_familywise'], v['loo_min'],
            v['nonblood_sep'], v['nonblood_frozen'], v['n_nonblood'],
            None if v['validation_primary_S_only'][0] is None else round(v['validation_primary_S_only'][0], 3),
            None if v['validation_primary_all_scorable'][0] is None else round(v['validation_primary_all_scorable'][0], 3),
            None if v['validation_safe_vs_hedged_toxic'][0] is None else round(v['validation_safe_vs_hedged_toxic'][0], 3)))
    print('robustness', {k: round(v, 3) for k, v in r['robustness_A3'].items()})
    print('promoted', r['decision']['promoted'])
