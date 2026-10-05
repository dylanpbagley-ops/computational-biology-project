"""Score the frozen rule on the validation set (target_labels_validation.py) and run the label-flip sensitivity for the main panel.
Precommitted in validation_precommitment.md."""
import json, os
import numpy as np
import target_labels as TL
import target_labels_validation as V
import lineage_map as LM
import flag_cleanup as FC
import wholebody_v2 as WB

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'validation_test.json')
WB.GENES.update(t[0] for t in V.VALIDATION)          # the loader only reads genes it is told about


def shares_all():
    vals, _ = WB.load('rna_single_cell_type.tsv.zip', ['nCPM', 'nTPM', 'pTPM'], ['Cell type'])
    S = {}
    for g, d in vals.items():
        v = np.array(list(d.values()), float)
        if v.sum() > 0:
            S[g] = (list(d.keys()), v / v.sum())
    return S


def rule_scores(S, genes):
    types = sorted({c for g in genes if g in S for c in S[g][0]})
    h0 = {c: 1.0 if TL.haematopoietic(c) else 0.0 for c in types}
    sub = {g: S[g] for g in genes if g in S}
    fl = FC.flags(sub, h0, LM.lineage)
    return {g: {'top_share': float(sub[g][1].max()), 'blood_share': float(fl[g][0]), 'blood_largest': bool(fl[g][1]),
                'home_lineage': fl[g][2], 'home': sub[g][0][int(np.argmax(sub[g][1]))],
                'rule': float(sub[g][1].max()) + (1.0 if fl[g][1] else 0.0)} for g in sub}


def sep(a, b):
    return WB.sep(a, b)


def exact_p_2sided_perm(a, b):
    p, exact = WB.exact_p(a, b)
    return p, exact


def main():
    S = shares_all()
    genes = [t[0] for t in V.VALIDATION]
    missing = [g for g in genes if g not in S]
    sc = rule_scores(S, genes)
    res = {'missing_genes': missing, 'targets': {}}
    for t in V.VALIDATION:
        g = t[0]
        if g in sc:
            res['targets'][g] = dict(sc[g], outcome=t[4], tier=t[5])
    prim_s = [g for g in V.PRIMARY_SAFE if g in sc]; prim_t = [g for g in V.PRIMARY_TOXIC if g in sc]
    hed = [g for g in V.HEDGED_TOXIC if g in sc]
    rule = lambda gs: [sc[g]['rule'] for g in gs]; ts = lambda gs: [sc[g]['top_share'] for g in gs]
    p, ex = exact_p_2sided_perm(rule(prim_s), rule(prim_t))
    res['primary'] = {'n_safe': len(prim_s), 'n_toxic': len(prim_t), 'sep_rule': sep(rule(prim_s), rule(prim_t)),
                      'sep_top_share_only': sep(ts(prim_s), ts(prim_t)), 'p_exact_permutation': float(p),
                      'min_attainable_p': 1.0 / float(np.math.comb(len(prim_s) + len(prim_t), len(prim_t))) if hasattr(np, 'math') else None}
    med = float(np.median(rule(prim_s)))
    res['primary']['toxic_below_median_safe'] = {g: bool(sc[g]['rule'] < med) for g in prim_t}
    res['hedged'] = {'genes': hed, 'sep_hedged_toxic_vs_primary_safe': sep(rule(prim_s), rule(hed)) if hed else None,
                     'sep_all_toxic_vs_primary_safe': sep(rule(prim_s), rule(prim_t + hed))}
    # pooled (NOT independent)
    d = json.load(open(os.path.join(HERE, 'wholebody_v2.json')))
    A = d['reference_A_cell_type']; hh = d['haem_home_A']; labels = d['labels']
    main_s = [g for g in A if labels.get(g) == 'safe']; main_t = [g for g in A if labels.get(g) == 'toxic']
    sc_main = {g: A[g]['top_share'] + (1.0 if hh[g] else 0.0) for g in A}
    res['pooled_not_independent'] = {'n_safe': len(main_s) + len(prim_s), 'n_toxic': len(main_t) + len(prim_t),
        'sep': sep([sc_main[g] for g in main_s] + rule(prim_s), [sc_main[g] for g in main_t] + rule(prim_t))}
    # label-flip sensitivity on the main panel (rule = top_share + home-flag; identical flags to M2 on all 30)
    flips = {'F0 none': ({}), 'F1 NECTIN4->toxic': {'NECTIN4': 'toxic'}, 'F2 F3,FOLR1->safe': {'F3': 'safe', 'FOLR1': 'safe'},
             'F3 F1+F2': {'NECTIN4': 'toxic', 'F3': 'safe', 'FOLR1': 'safe'}, 'F4 TNFRSF17->toxic': {'TNFRSF17': 'toxic'}}
    res['flip_sensitivity'] = {}
    for name, fl in flips.items():
        lab = dict(labels); lab.update(fl)
        s_ = [g for g in A if lab.get(g) == 'safe']; t_ = [g for g in A if lab.get(g) == 'toxic']
        p_, ex_ = WB.exact_p([sc_main[g] for g in s_], [sc_main[g] for g in t_])
        res['flip_sensitivity'][name] = {'n_safe': len(s_), 'n_toxic': len(t_), 'sep': sep([sc_main[g] for g in s_], [sc_main[g] for g in t_]),
                                         'one_factor_sep': sep([A[g]['top_share'] for g in s_], [A[g]['top_share'] for g in t_]), 'p': float(p_)}
    json.dump(res, open(OUT, 'w'), indent=1)
    print('missing genes:', missing)
    for g, t in res['targets'].items():
        print('%-8s %-6s %-8s home=%-30s top=%.2f blood=%.2f lin=%-14s rule=%.2f' % (g, t['outcome'], t['tier'], t['home'][:30], t['top_share'], t['blood_share'], t['home_lineage'], t['rule']))
    print(json.dumps({k: res[k] for k in ('primary', 'hedged', 'pooled_not_independent')}, indent=1))
    for k, v in res['flip_sensitivity'].items():
        print(k, v)


if __name__ == '__main__':
    main()
