"""
WHOLE-BODY RESTRICTION, v2 -- expanded panel, two independent normal references, and the
two-factor test that turns the project's thesis into something falsifiable.

v1 (wholebody_restriction.py / wholebody_toxicity_test.py) showed that CONFINEMENT across
154 normal cell types separates safe from on-target-toxic targets. It also showed the one
case that breaks a pure-confinement account: CEACAM5 is confined and still caused colitis.

v2 asks the question that case forces:
    is 'confined AND confined to an expendable lineage' actually better than 'confined'?
The expendability factor is assigned from the IDENTITY of the home cell type alone
(haematopoietic = continuously regenerated from marrow with an established clinical
replacement), never from the drug outcome, so the test is not circular.

Labels live in target_labels.py and are never edited here.
References:
  A. HPA single-cell consensus  rna_single_cell_type.tsv.zip   (154 cell types)
  B. HPA bulk tissue consensus  rna_tissue_consensus.tsv.zip   (independent, tissue level)
"""
import os, io, json, zipfile, collections, itertools, numpy as np
import target_labels as TL

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'ev_capture_data', 'normal_tissue_atlas'))
GENES = {t[0] for t in TL.TARGETS}


def load(zip_name, unit_col_candidates, group_col_candidates):
    z = zipfile.ZipFile(os.path.join(D, zip_name))
    name = z.namelist()[0]
    vals = collections.defaultdict(dict)
    with z.open(name) as f:
        fh = io.TextIOWrapper(f, 'utf-8')
        hdr = fh.readline().rstrip('\n').split('\t')
        gi = hdr.index('Gene name')
        ci = next(hdr.index(c) for c in group_col_candidates if c in hdr)
        vi = next(hdr.index(c) for c in unit_col_candidates if c in hdr)
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if p[gi] in GENES:
                try:
                    vals[p[gi]][p[ci]] = float(p[vi])
                except (ValueError, IndexError):
                    pass
    return vals, hdr


def confinement(vals):
    out = {}
    for g, d in vals.items():
        v = np.array(list(d.values()), float)
        ks = list(d.keys())
        tot = v.sum()
        if tot <= 0:
            continue
        top = int(np.argmax(v))
        order = np.argsort(-v)[:3]
        out[g] = {
            'top_share': float(v[top] / tot),
            'n_types': int((v > 0.10 * v[top]).sum()),
            'home': ks[top],
            'top3': [(ks[int(i)], float(v[int(i)] / tot)) for i in order],
        }
    return out


# ---- rank-separation (= AUC) and its exact permutation p-value ---------------------------
def sep(a, b):
    a, b = list(a), list(b)
    gt = sum(1 for x in a for y in b if x > y)
    eq = sum(1 for x in a for y in b if x == y)
    return (gt + 0.5 * eq) / (len(a) * len(b))


def exact_p(a, b, max_exact=250000):
    """One-sided exact permutation p for sep(a,b) >= observed. Falls back to a large
    random permutation sample when the exact enumeration is too big; the fallback is
    reported so it is never mistaken for an exact result."""
    a, b = list(a), list(b)
    allv = a + b
    n, k = len(allv), len(a)
    obs = sep(a, b)
    from math import comb
    total = comb(n, k)
    if total <= max_exact:
        c = 0
        for idx in itertools.combinations(range(n), k):
            s = set(idx)
            aa = [allv[i] for i in idx]
            bb = [allv[i] for i in range(n) if i not in s]
            c += (sep(aa, bb) >= obs - 1e-12)
        return c / total, True
    rng = np.random.default_rng(0)
    arr = np.array(allv, float)
    c = 0
    R = 200000
    for _ in range(R):
        p = rng.permutation(n)
        c += (sep(arr[p[:k]], arr[p[k:]]) >= obs - 1e-12)
    return (c + 1) / (R + 1), False


def run_test(scores, label, higher_is_safer=True, pos='safe', neg='toxic'):
    A = [scores[g] for g in scores if label[g] == pos]
    B = [scores[g] for g in scores if label[g] == neg]
    if not higher_is_safer:
        A = [-x for x in A]
        B = [-x for x in B]
    s = sep(A, B)
    p, ex = exact_p(A, B)
    return s, p, ex, len(A), len(B)


def main():
    print('=' * 78)
    print('REFERENCE A: HPA single-cell consensus (cell-type level)')
    print('=' * 78)
    vA, _ = load('rna_single_cell_type.tsv.zip', ['nCPM', 'nTPM', 'pTPM'], ['Cell type'])
    cA = confinement(vA)
    missing = sorted(GENES - set(cA))
    print('genes resolved: %d / %d  %s' % (len(cA), len(GENES),
          ('MISSING: ' + ', '.join(missing)) if missing else ''))

    label = {t[0]: t[3] for t in TL.TARGETS}
    common = {t[0]: t[1] for t in TL.TARGETS}

    rows = sorted(cA.items(), key=lambda kv: -kv[1]['top_share'])
    print('\n%-9s %-13s %9s %8s %-30s %-6s %s' %
          ('GENE', 'common', 'top_share', 'n_types', 'home cell type', 'haem?', 'outcome'))
    for g, d in rows:
        print('%-9s %-13s %9.3f %8d %-30s %-6s %s' %
              (g, common[g], d['top_share'], d['n_types'], d['home'][:30],
               'YES' if TL.haematopoietic(d['home']) else 'no', label[g]))

    # ---- restrict to the outcome-bearing set ------------------------------------------
    use = {g: d for g, d in cA.items() if label[g] in ('safe', 'toxic')}
    print('\n' + '-' * 78)
    print('TEST 1 (one factor): does CONFINEMENT alone separate safe from toxic?')
    print('-' * 78)
    for key, hs, desc in (('top_share', True, 'top_share  (higher = more confined)'),
                          ('n_types', False, 'n_types    (fewer  = more confined)')):
        sc = {g: d[key] for g, d in use.items()}
        s, p, ex, na, nb = run_test(sc, label, hs)
        print('  %-38s separation %.3f  %s p = %.5f   (safe n=%d vs toxic n=%d)'
              % (desc, s, 'exact' if ex else 'approx', p, na, nb))

    print('\n' + '-' * 78)
    print('TEST 2 (two factors): confinement AND an expendable (haematopoietic) home')
    print('-' * 78)
    print('  Rule, fixed in advance: score = top_share, PLUS 1.0 if the home cell type is')
    print('  haematopoietic. The +1 cannot be beaten by any confinement value, so this is')
    print('  strictly "expendable first, confinement second" -- a lexicographic rule with')
    print('  no fitted parameter.')
    sc2 = {g: d['top_share'] + (1.0 if TL.haematopoietic(d['home']) else 0.0)
           for g, d in use.items()}
    s2, p2, ex2, na, nb = run_test(sc2, label, True)
    print('  two-factor score                       separation %.3f  %s p = %.5f'
          % (s2, 'exact' if ex2 else 'approx', p2))

    sc1 = {g: d['top_share'] for g, d in use.items()}
    s1, p1, _, _, _ = run_test(sc1, label, True)
    print('\n  one factor  %.3f   ->   two factors  %.3f     (change %+.3f)' % (s1, s2, s2 - s1))

    # which targets does each rule get wrong?
    def errors(sc):
        A = [(g, sc[g]) for g in sc if label[g] == 'safe']
        B = [(g, sc[g]) for g in sc if label[g] == 'toxic']
        bad = []
        for gt, vt in B:
            for gs, vs in A:
                if vt > vs:
                    bad.append((gt, gs))
        return bad
    e1, e2 = errors(sc1), errors(sc2)
    print('\n  mis-ordered (toxic ranked safer than a safe target):')
    print('    one factor : %d pairs, involving toxic targets %s'
          % (len(e1), sorted({common[a] for a, _ in e1})))
    print('    two factors: %d pairs, involving toxic targets %s'
          % (len(e2), sorted({common[a] for a, _ in e2})))

    # ---- leave-one-out robustness -----------------------------------------------------
    print('\n' + '-' * 78)
    print('TEST 3: leave-one-out. Drop each target in turn and re-run the two-factor test.')
    print('-' * 78)
    seps, ps = [], []
    for drop in list(use):
        sub = {g: v for g, v in sc2.items() if g != drop}
        if len({label[g] for g in sub}) < 2:
            continue
        s, p, _, _, _ = run_test(sub, label, True)
        seps.append(s); ps.append(p)
    print('  separation  min %.3f  median %.3f  max %.3f' % (min(seps), float(np.median(seps)), max(seps)))
    print('  p-value     min %.5f  median %.5f  max %.5f  |  p<0.05 in %d/%d drops'
          % (min(ps), float(np.median(ps)), max(ps), sum(1 for p in ps if p < 0.05), len(ps)))

    # ---- efficacy failures as a specificity control ------------------------------------
    print('\n' + '-' * 78)
    print('TEST 4 (specificity): targets that failed for EFFICACY, not toxicity, should NOT')
    print('        look toxic. If the measure flagged them too, it would be measuring')
    print('        "was ever tried and did not work", not on-target toxicity.')
    print('-' * 78)
    eff = {g: d for g, d in cA.items() if label[g] == 'efficacy_fail'}
    safe = {g: d for g, d in cA.items() if label[g] == 'safe'}
    tox = {g: d for g, d in cA.items() if label[g] == 'toxic'}
    def f(dd, k): return [v[k] for v in dd.values()]
    print('  mean top_share   safe %.3f | efficacy-fail %.3f | toxic %.3f'
          % (np.mean(f(safe, 'top_share')), np.mean(f(eff, 'top_share')), np.mean(f(tox, 'top_share'))))
    s_ef, p_ef, ex_ef, na, nb = run_test(
        {g: d['top_share'] for g, d in list(safe.items()) + list(eff.items())},
        label, True, pos='safe', neg='efficacy_fail')
    print('  safe vs efficacy-fail: separation %.3f  %s p = %.5f  (n=%d vs %d)'
          % (s_ef, 'exact' if ex_ef else 'approx', p_ef, na, nb))
    print('  (this one SHOULD be near 0.5 and non-significant.)')


    # ---- TEST 5: is the result an artifact of germ / placental cell types? --------------
    print('\n' + '-' * 78)
    print('TEST 5 (sensitivity): several targets have their single highest cell type in germ')
    print('        or placental tissue (CA9 -> spermatids, DLL3 -> oocytes, EGFR ->')
    print('        syncytiotrophoblasts). Those sites are immune-privileged and are NOT where')
    print('        the clinical toxicity fell. Drop every germ/placental/reproductive cell')
    print('        type from the reference and recompute from scratch.')
    print('-' * 78)
    DROP = ('spermat', 'oocyte', 'sperm', 'trophoblast', 'syncytio', 'placenta', 'decidual',
            'granulosa', 'leydig', 'sertoli', 'epididym', 'seminal', 'germ', 'ductus deferens',
            'extravillous')
    vA_clean = {g: {k: v for k, v in d.items() if not any(w in k.lower() for w in DROP)}
                for g, d in vA.items()}
    dropped = len(next(iter(vA.values()))) - len(next(iter(vA_clean.values())))
    cA_clean = confinement(vA_clean)
    print('  removed %d of %d cell types' % (dropped, len(next(iter(vA.values())))))
    useC = {g: d for g, d in cA_clean.items() if label[g] in ('safe', 'toxic')}
    s1c, p1c, ex1c, na, nb = run_test({g: d['top_share'] for g, d in useC.items()}, label, True)
    sc2c = {g: d['top_share'] + (1.0 if TL.haematopoietic(d['home']) else 0.0)
            for g, d in useC.items()}
    s2c, p2c, ex2c, _, _ = run_test(sc2c, label, True)
    print('  one factor  separation %.3f  p = %.5f   (was %.3f)' % (s1c, p1c, s1))
    print('  two factors separation %.3f  p = %.5f   (was %.3f)' % (s2c, p2c, s2))
    moved = [(g, cA[g]['home'], cA_clean[g]['home']) for g in cA_clean
             if cA[g]['home'] != cA_clean[g]['home']]
    print('  targets whose home cell type changed: %d' % len(moved))
    for g, a, b in sorted(moved):
        print('     %-9s %-30s -> %s' % (g, a[:30], b))

    # ---- TEST 6: does it hold for SOLID-TUMOUR targets alone? --------------------------
    print('\n' + '-' * 78)
    print('TEST 6 (the load-bearing one): the obvious objection is that this measure only')
    print('        works because blood targets are confined and blood targets are the safe')
    print('        ones. Remove EVERY blood-cancer target and re-run on solid-tumour')
    print('        targets only. If it survives, the finding is not a blood artifact.')
    print('-' * 78)
    BLOOD = {'CD19', 'MS4A1', 'CD22', 'CD79B', 'TNFRSF17', 'SLAMF7', 'CD38', 'TNFRSF8',
             'CD33', 'GPRC5D', 'CD52', 'PTPRC', 'IL3RA', 'FCER2', 'CD40', 'CD80',
             'HLA-DRA', 'SDC1', 'CD74', 'CD47'}
    useS = {g: d for g, d in cA.items() if label[g] in ('safe', 'toxic') and g not in BLOOD}
    solid_safe = sorted(g for g in useS if label[g] == 'safe')
    solid_tox = sorted(g for g in useS if label[g] == 'toxic')
    print('  solid safe  (n=%d): %s' % (len(solid_safe), ', '.join(common[g] for g in solid_safe)))
    print('  solid toxic (n=%d): %s' % (len(solid_tox), ', '.join(common[g] for g in solid_tox)))
    sS, pS, exS, na, nb = run_test({g: d['top_share'] for g, d in useS.items()}, label, True)
    sSn, pSn, _, _, _ = run_test({g: d['n_types'] for g, d in useS.items()}, label, False)
    print('  top_share separation %.3f  %s p = %.5f' % (sS, 'exact' if exS else 'approx', pS))
    print('  n_types   separation %.3f  p = %.5f' % (sSn, pSn))
    print('  (no haematopoietic target remains, so the two-factor rule reduces to the one-')
    print('   factor rule here by construction, and is not reported separately.)')

    # ---- reference B, independent -----------------------------------------------------
    print('\n' + '=' * 78)
    print('REFERENCE B: HPA bulk tissue consensus -- an INDEPENDENT normal reference')
    print('=' * 78)
    vB, hdrB = load('rna_tissue_consensus.tsv.zip', ['nTPM', 'nCPM', 'pTPM'], ['Tissue'])
    cB = confinement(vB)
    print('columns: %s' % hdrB)
    print('genes resolved: %d / %d ; tissues: %d'
          % (len(cB), len(GENES), len(next(iter(vB.values())))))
    useB = {g: d for g, d in cB.items() if label[g] in ('safe', 'toxic')}
    for key, hs, desc in (('top_share', True, 'top_share'), ('n_types', False, 'n_types')):
        sc = {g: d[key] for g, d in useB.items()}
        s, p, ex, na, nb = run_test(sc, label, hs)
        print('  one factor %-10s separation %.3f  %s p = %.5f  (n=%d vs %d)'
              % (desc, s, 'exact' if ex else 'approx', p, na, nb))
    scB2 = {g: d['top_share'] + (1.0 if TL.haematopoietic(d['home']) else 0.0)
            for g, d in useB.items()}
    n_haem_B = sum(1 for g, d in useB.items() if TL.haematopoietic(d['home']))
    sB2, pB2, exB2, na, nb = run_test(scB2, label, True)
    print('  two factors           separation %.3f  %s p = %.5f'
          % (sB2, 'exact' if exB2 else 'approx', pB2))
    print('  NOTE: reference B is indexed by TISSUE (bone marrow, lymph node), not by cell')
    print('  type, so the haematopoietic flag fires on %d of %d targets and the two-factor'
          % (n_haem_B, len(useB)))
    print('  rule collapses onto the one-factor rule. Reference B tests CONFINEMENT ONLY.')
    print('  It is also weaker by construction: 51 tissues cannot resolve a B cell from the')
    print('  lymph node it sits in, which is exactly why the cell-type reference is primary.')

    # ---- agreement between the two references -----------------------------------------
    shared = sorted(set(cA) & set(cB))
    xa = np.array([cA[g]['top_share'] for g in shared])
    xb = np.array([cB[g]['top_share'] for g in shared])
    from scipy.stats import spearmanr
    rho, prho = spearmanr(xa, xb)
    print('\n  agreement between the two independent references: Spearman rho = %.3f (p = %.2g, n = %d)'
          % (rho, prho, len(shared)))

    json.dump({'reference_A_cell_type': cA, 'reference_A_germ_removed': cA_clean,
               'solid_only': {'sep': sS, 'p': pS, 'n_safe': len(solid_safe), 'n_toxic': len(solid_tox)},
               'germ_removed': {'one_factor': s1c, 'two_factor': s2c, 'p_two': p2c}, 'reference_B_tissue': cB,
               'labels': {g: label[g] for g in label},
               'haem_home_A': {g: bool(TL.haematopoietic(d['home'])) for g, d in cA.items()},
               'tests': {'A_one_factor_top_share': s1, 'A_two_factor': s2, 'A_two_factor_p': p2,
                         'B_two_factor': sB2, 'B_two_factor_p': pB2,
                         'loo_sep_min': min(seps), 'loo_p_max': max(ps),
                         'spearman_A_vs_B': float(rho)}},
              open(os.path.join(HERE, 'wholebody_v2.json'), 'w'), indent=1)
    print('\nwrote wholebody_v2.json')


if __name__ == '__main__':
    main()
