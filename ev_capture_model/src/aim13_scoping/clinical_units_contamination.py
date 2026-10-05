"""Why the clinical-units ceiling collapses at low harm budgets, tested directly.

HYPOTHESIS (stated before the run)
----------------------------------
clinical_units_dlbcl.py shows CD19 kills 89.5% of tumour while sparing 99% of T cells, but only
~1% while sparing 99% of 'other leukocytes' (O) -- i.e. chance. The proposed mechanism: a harm
budget b sits at the (1-b) quantile of the normal compartment. If a few percent of that
compartment are actually B-lineage cells that slipped past the gate, the quantile lands inside
the tumour's own distribution and 'kill' collapses to ~b regardless of the antigen.

If true, the low-budget number is measuring GATING PURITY, not antigen biology.

TEST
----
1. Estimate B-lineage contamination of T and O using ONLY B-lineage markers other than the
   target (leave-target-out), so the estimate is not circular: a cell is 'B-like' if at least
   K_MIN of the other B-lineage markers exceed the tumour's PCT-th percentile.
2. Remove those cells and recompute kill-at-harm. If the hypothesis is right, the O-compartment
   numbers should jump toward the T-compartment numbers.
3. Sweep PCT and K_MIN so the conclusion is shown not to depend on one threshold.

Markers OUTSIDE the B-lineage set (CD40, CD80, HLA-DR ...) are decontaminated with all four
B-lineage markers. Normal B cells (NB) are never decontaminated: they are B cells legitimately,
they are the expendable compartment by design.

WHAT WOULD FALSIFY IT: if removing B-like cells leaves O-compartment kill at ~chance, the
collapse has some other cause and this explanation must be dropped.
"""
import json
import os
import numpy as np

import clinical_units_dlbcl as D

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'clinical_units_contamination.json')
B_LINEAGE = ('CD19', 'CD20', 'CD22', 'CD79B')
PCT_DEFAULT, K_DEFAULT = 25, 2
PCT_SWEEP = (10, 25, 50)
K_SWEEP = (1, 2, 3)
TARGETS = ('CD19', 'CD20', 'CD22', 'CD79B', 'CD40', 'CD80', 'HLA-DR')


def blike(pt, comp, ix, target, pct, k):
    others = [m for m in B_LINEAGE if m != target]
    cnt = np.zeros(len(comp), int)
    for m in others:
        cnt += comp[:, ix[m]] > np.percentile(pt['B'][:, ix[m]], pct)
    return cnt >= min(k, len(others))


def run(pct=PCT_DEFAULT, k=K_DEFAULT, budgets=D.BUDGETS):
    ant, pts = D.load()
    ix = {a: i for i, a in enumerate(ant)}
    res = {}
    for tgt in TARGETS:
        contam = {'T': [], 'O': []}
        before = {'T': [], 'O': [], 'NB': []}
        after = {'T': [], 'O': [], 'NB': []}
        for pt in pts:
            t = pt['B'][:, ix[tgt]]
            for lab in ('T', 'O', 'NB'):
                comp = pt['comp'][lab]
                before[lab].append(D.kills(t, comp[:, ix[tgt]], budgets))
                if lab == 'NB':
                    after[lab].append(before[lab][-1])
                    continue
                m = blike(pt, comp, ix, tgt, pct, k)
                contam[lab].append(float(m.mean()))
                keep = comp[~m]
                after[lab].append(D.kills(t, keep[:, ix[tgt]], budgets) if len(keep) > 200
                                  else np.full(len(budgets), np.nan))
        med = lambda L: np.nanmedian(np.array(L), axis=0)
        res[tgt] = {
            'contamination_frac_median': {l: float(np.median(v)) for l, v in contam.items()},
            'kill_before': {l: [float(x) for x in med(v)] for l, v in before.items()},
            'kill_after': {l: [float(x) for x in med(v)] for l, v in after.items()},
        }
        for tag in ('before', 'after'):
            kk = res[tgt][f'kill_{tag}']
            res[tgt][f'ceiling_nonexp_{tag}'] = [float(min(kk['T'][j], kk['O'][j]))
                                                 for j in range(len(budgets))]
    return res


def compute_all():
    main = run()
    sweep = {}
    for pct in PCT_SWEEP:
        for k in K_SWEEP:
            r = run(pct, k)
            sweep[f'pct{pct}_k{k}'] = {t: {'O_after_at_1pct': r[t]['kill_after']['O'][0],
                                           'O_contam': r[t]['contamination_frac_median']['O']}
                                       for t in ('CD19', 'CD20', 'CD22', 'CD79B')}
    return {'main': main, 'sweep': sweep, 'meta': {'budgets': list(D.BUDGETS),
            'pct': PCT_DEFAULT, 'k_min': K_DEFAULT, 'b_lineage': list(B_LINEAGE)}}


def print_summary(r):
    print('=' * 100)
    print('DECONTAMINATING THE NORMAL COMPARTMENTS (leave-target-out B-lineage consensus)')
    print('budgets:', [f'{b:.0%}' for b in D.BUDGETS])
    print('=' * 100)
    for t, v in r['main'].items():
        c = v['contamination_frac_median']
        print(f"\n{t:7s} B-like fraction  T: {c['T']:.3f}   O: {c['O']:.3f}")
        for lab in ('T', 'O'):
            print(f"   {lab}: kill before {[round(x,3) for x in v['kill_before'][lab]]}"
                  f"  after {[round(x,3) for x in v['kill_after'][lab]]}")
        print(f"   non-expendable ceiling  before {[round(x,3) for x in v['ceiling_nonexp_before']]}"
              f"  after {[round(x,3) for x in v['ceiling_nonexp_after']]}")
    print('\nROBUSTNESS -- O-compartment kill at 1% harm after decontamination, by threshold')
    for key, v in r['sweep'].items():
        print(f"  {key:9s} " + '  '.join(f"{t}: {x['O_after_at_1pct']:.2f} (contam {x['O_contam']:.3f})"
                                          for t, x in v.items()))


if __name__ == '__main__':
    r = compute_all()
    print_summary(r)
    json.dump(r, open(OUT, 'w'), indent=1)
    print('\nwrote', OUT)
