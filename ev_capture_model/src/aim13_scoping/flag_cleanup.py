"""Round 4A: majority / plurality blood flag. Precommitted in flag_cleanup_precommitment.md."""
import json, os, collections
import numpy as np
import target_labels as TL
import lineage_map as LM
import continuous_weights_v2 as C2
import wholebody_v2 as WB

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'flag_cleanup.json')
N_SHUF, N_BOOT, SEED = 1000, 5000, 20261001
SCHEMES = ('M1', 'M2')


def fisher_p(a, b, c, d):
    from math import comb
    n1, n2, k = a + b, c + d, a + c
    tot = comb(n1 + n2, k)
    def pr(x): return comb(n1, x) * comb(n2, k - x) / tot
    p0 = pr(a)
    return sum(pr(x) for x in range(max(0, k - n2), min(n1, k) + 1) if pr(x) <= p0 + 1e-12)


def flags(S, hmap, lin):
    """Return {gene: (B, blood_is_largest_lineage)} for a blood indicator map hmap and lineage function lin."""
    out = {}
    for g, (types, s) in S.items():
        B = float(sum(si * hmap[c] for c, si in zip(types, s)))
        agg = collections.defaultdict(float)
        for c, si in zip(types, s):
            agg['blood' if hmap[c] else lin(c)] += si
        top_lin = max(agg, key=agg.get)
        out[g] = (B, top_lin == 'blood', top_lin)
    return out


def score(S, hmap, lin, sch):
    fl = flags(S, hmap, lin)
    return {g: float(S[g][1].max()) + (1.0 if ((fl[g][0] >= 0.5) if sch == 'M1' else fl[g][1]) else 0.0) for g in S}


def compute_all():
    d = json.load(open(os.path.join(HERE, 'wholebody_v2.json')))
    labels, A, hh = d['labels'], d['reference_A_cell_type'], d['haem_home_A']
    S = {g: v for g, v in C2.shares().items() if g in A}
    safe = [g for g in S if labels.get(g) == 'safe']; tox = [g for g in S if labels.get(g) == 'toxic']
    types = sorted({c for t, _ in S.values() for c in t}); LM.check(types)
    h0 = {c: 1.0 if TL.haematopoietic(c) else 0.0 for c in types}
    rng = np.random.default_rng(SEED)
    binary = {g: A[g]['top_share'] + (1.0 if hh[g] else 0.0) for g in A}
    one = {g: A[g]['top_share'] for g in A}
    base = {'one_factor': float(C2.sep_of(one, labels, safe, tox)[0]),
            'binary_home_rule': float(C2.sep_of(binary, labels, safe, tox)[0])}
    fl = flags(S, h0, LM.lineage)
    table = {g: {'label': labels[g], 'home': A[g]['home'], 'home_flag': bool(hh[g]), 'top_share': round(A[g]['top_share'], 3),
                 'blood_share': round(fl[g][0], 3), 'M1_flag': bool(fl[g][0] >= 0.5), 'M2_flag': bool(fl[g][1]),
                 'home_lineage': fl[g][2]} for g in safe + tox}
    res = {'baselines': base, 'targets': table, 'schemes': {}, 'flag_counts': {}}
    for k, key in (('home_rule', 'home_flag'), ('M1', 'M1_flag'), ('M2', 'M2_flag')):
        ns = sum(table[g][key] for g in safe); nt = sum(table[g][key] for g in tox)
        res['flag_counts'][k] = {'safe_flagged': ns, 'toxic_flagged': nt, 'fisher_p': fisher_p(ns, len(safe) - ns, nt, len(tox) - nt)}
    for sch in SCHEMES:
        sc = score(S, h0, LM.lineage, sch)
        s, a, b = C2.sep_of(sc, labels, safe, tox); p, exact = WB.exact_p(a, b)
        nul = np.empty(N_SHUF)
        hv = np.array([h0[c] for c in types])
        lin_null = lambda c: 'shadow' if TL.haematopoietic(c) else LM.lineage(c)
        for i in range(N_SHUF):
            m = dict(zip(types, rng.permutation(hv)))
            # null lineage: originally-blood types must NOT stay in the blood family when their flag is shuffled away
            nul[i] = C2.sep_of(score(S, m, lin_null, sch), labels, safe, tox)[0]
        boots = np.empty(N_BOOT); diffs = np.empty(N_BOOT)
        for i in range(N_BOOT):
            sa = list(rng.choice(safe, len(safe))); sb = list(rng.choice(tox, len(tox)))
            x = WB.sep([sc[g] for g in sa], [sc[g] for g in sb]); y = WB.sep([binary[g] for g in sa], [binary[g] for g in sb])
            boots[i] = x; diffs[i] = x - y
        loo = [C2.sep_of(sc, labels, [x for x in safe if x != g], [x for x in tox if x != g])[0] for g in safe + tox]
        res['schemes'][sch] = {'sep': float(s), 'p_label_permutation': float(p),
            'ci': [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))],
            'diff_vs_home_rule': float(s - base['binary_home_rule']),
            'diff_ci': [float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))],
            'p_vs_shuffled_flags': float((np.sum(nul >= s - 1e-12) + 1) / (N_SHUF + 1)), 'loo_min': float(min(loo))}
    ok = [k for k, v in res['schemes'].items() if v['sep'] >= 0.900 and v['diff_ci'][0] <= 0 <= v['diff_ci'][1]]
    res['decision'] = {'rule': 'cleaned flag replaces home rule if sep >= 0.900 and paired CI includes 0',
                       'schemes_meeting_rule': ok, 'cleanup_adopted': len(ok) > 0}
    return res


def main():
    r = compute_all()
    print('baselines', {k: round(v, 3) for k, v in r['baselines'].items()})
    for k, v in r['flag_counts'].items():
        print('flag', k, v)
    for k, v in r['schemes'].items():
        print(k, 'sep %.3f p %.5f CI [%.3f,%.3f] vs home %+.3f [%+.3f,%+.3f] shufflep %.3f LOOmin %.3f' % (
            v['sep'], v['p_label_permutation'], *v['ci'], v['diff_vs_home_rule'], *v['diff_ci'], v['p_vs_shuffled_flags'], v['loo_min']))
    print('decision', r['decision'])
    for g, t in sorted(r['targets'].items(), key=lambda x: (x[1]['label'], -x[1]['blood_share'])):
        print('%-9s %-5s home=%-28s top=%.2f B=%.2f homeflag=%d M1=%d M2=%d lin=%s' % (g, t['label'], t['home'][:28], t['top_share'], t['blood_share'], t['home_flag'], t['M1_flag'], t['M2_flag'], t['home_lineage']))
    json.dump(r, open(OUT, 'w'), indent=1)


if __name__ == '__main__':
    main()
