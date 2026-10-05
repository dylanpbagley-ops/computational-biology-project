"""The ceiling in clinical units, on data where the antigen actually separates the populations.

WHY THIS EXISTS
---------------
clinical_units_noise_floor.py showed that on the myeloma CITE-seq set the kill-at-fixed-harm
metric sits at the noise floor for essentially every marker, because marker-level tumour-vs-
normal separation there is weak (published AUCs: BCMA 0.61, CD38 0.29). A tail metric cannot be
informative when there is almost no separation to read. This repeats the analysis on the DLBCL
CyTOF cohort, where the B-lineage antigens separate B cells from other leukocytes strongly, and
asks the two questions the metric exists to answer:

  Q1. For each antigen, what is the most tumour that could be killed while sparing 99% / 95% of
      each normal compartment (T cells, other leukocytes, batch-matched normal B cells)?
  Q2. Does that number, computed with normal B cells treated as expendable, line up with what
      happened clinically -- approved (CD19, CD20, CD22, CD79B) vs. failed (CD40, CD80, HLA-DR)?
      This is DESCRIPTIVE (n = 4 vs 3); no p-value is claimed, none would be meaningful.

The 'tumour' population is the QC-gated B-cell population per patient, the same one used for
three_compartment_ceiling.py, so results are directly comparable to that published ceiling.
Null envelope: labels shuffled within each patient/compartment, 30 permutations.

TWO CEILINGS PER ANTIGEN
------------------------
  all-compartments   min over T, other leukocytes and normal B  (the strict bound)
  non-expendable     min over T and other leukocytes only (normal B cells treated as expendable,
                     as B-cell aplasia is managed clinically with immunoglobulin replacement)
The gap between them is the direct cost of sparing normal B cells -- the project's central
expendability claim expressed in units of tumour killed.
"""
import glob
import json
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'ev_capture_data', 'gated_cache'))
OUT = os.path.join(HERE, 'clinical_units_dlbcl.json')

BUDGETS = (0.01, 0.05, 0.10, 0.20)
SUB = 12000
SEED = 31          # same subsample seed family as three_compartment_ceiling.py
N_PERM = 30
NULL_SEED = 20260928
APPROVED = ('CD19', 'CD20', 'CD22', 'CD79B')
FAILED = ('CD40', 'CD80', 'HLA-DR')


def kills(t, n, budgets=BUDGETS):
    """Max tumour fraction killable at each harm budget (threshold at the 1-b quantile of normal)."""
    th = np.quantile(n, [1.0 - b for b in budgets])
    return np.array([(t > x).mean() for x in th])


def auc(x, y):
    n1, n2 = len(x), len(y)
    z = np.concatenate([x, y])
    order = np.argsort(z, kind='mergesort')
    zs = z[order]
    r = np.empty(len(z))
    i = 0
    while i < len(z):
        j = i
        while j + 1 < len(z) and zs[j + 1] == zs[i]:
            j += 1
        r[i:j + 1] = (i + j) / 2.0 + 1
        i = j + 1
    rr = np.empty(len(z))
    rr[order] = r
    return (rr[:n1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n2)


def load():
    files = {os.path.basename(p)[:-4]: p for p in glob.glob(os.path.join(CACHE, '*.npz'))}
    ant = [str(a) for a in np.load(next(iter(files.values())), allow_pickle=True)['antigens']]
    bat = lambda k: k.split('_')[0]
    tumors = sorted(k for k in files if '_DLB_' in k)
    rlns = sorted(k for k in files if '_RLN_' in k and 'pooled' not in k)
    normal_b = {b: np.concatenate([np.load(files[k], allow_pickle=True)['B'] for k in rlns
                                   if bat(k) == b], axis=0).astype(float)
                for b in sorted({bat(k) for k in rlns})}
    rng = np.random.default_rng(SEED)
    sub = lambda a: a[rng.choice(len(a), SUB, replace=False)] if len(a) > SUB else a
    pts = []
    for t in tumors:
        d = np.load(files[t], allow_pickle=True)
        B, T, O = d['B'].astype(float), d['T'].astype(float), d['O'].astype(float)
        if min(len(B), len(T), len(O)) < 500 or bat(t) not in normal_b:
            continue
        pts.append({'id': t, 'B': sub(B), 'comp': {'T': sub(T), 'O': sub(O), 'NB': sub(normal_b[bat(t)])}})
    return ant, pts


def per_marker(ant, pts, perm_rng=None):
    """Median-over-patients kill per (marker, compartment, budget); with perm_rng labels are
    shuffled within each patient/compartment first."""
    out = np.zeros((len(ant), 3, len(BUDGETS)))
    for c, lab in enumerate(('T', 'O', 'NB')):
        acc = np.zeros((len(pts), len(ant), len(BUDGETS)))
        for p, pt in enumerate(pts):
            for i in range(len(ant)):
                t, n = pt['B'][:, i], pt['comp'][lab][:, i]
                if perm_rng is not None:
                    pool = np.concatenate([t, n])
                    perm_rng.shuffle(pool)
                    t, n = pool[:len(t)], pool[len(t):]
                acc[p, i] = kills(t, n)
        out[:, c, :] = np.median(acc, axis=0)
    return out


def compute_all(n_perm=N_PERM):
    ant, pts = load()
    obs = per_marker(ant, pts)                     # marker x comp x budget
    rng = np.random.default_rng(NULL_SEED)
    null_all = np.zeros((n_perm, len(ant), len(BUDGETS)))
    for k in range(n_perm):
        pm = per_marker(ant, pts, rng)
        null_all[k] = pm.min(axis=1)
    p95 = np.percentile(null_all, 95, axis=0)      # marker x budget
    res = {}
    for i, a in enumerate(ant):
        row = {}
        for j, b in enumerate(BUDGETS):
            ce_all = float(obs[i, :, j].min())
            ce_ne = float(obs[i, :2, j].min())
            row[f'{b:.2f}'] = {
                'per_compartment': {'T': float(obs[i, 0, j]), 'O': float(obs[i, 1, j]),
                                    'NB': float(obs[i, 2, j])},
                'ceiling_all': ce_all, 'ceiling_nonexpendable': ce_ne,
                'cost_of_sparing_normal_B': ce_ne - ce_all,
                'null_p95_all': float(p95[i, j]), 'above_null_all': bool(ce_all > p95[i, j]),
            }
        res[a] = row
    summary = {}
    for b in BUDGETS:
        k = f'{b:.2f}'
        summary[k] = {
            'n_above_null': sum(res[a][k]['above_null_all'] for a in ant),
            'n_markers': len(ant),
            'approved_nonexp': {a: res[a][k]['ceiling_nonexpendable'] for a in APPROVED if a in res},
            'failed_nonexp': {a: res[a][k]['ceiling_nonexpendable'] for a in FAILED if a in res},
            'approved_all': {a: res[a][k]['ceiling_all'] for a in APPROVED if a in res},
            'failed_all': {a: res[a][k]['ceiling_all'] for a in FAILED if a in res},
        }
    return {'per_marker': res, 'summary': summary,
            'meta': {'n_patients': len(pts), 'budgets': list(BUDGETS), 'n_perm': n_perm,
                     'sub': SUB, 'seed': SEED, 'note': 'fractions of a population; calibration-invariant'}}


def crosscheck_against_published(ant, pts, tol=0.03):
    """Per-compartment median AUC from THESE populations vs three_compartment_results.json."""
    pub = json.load(open(os.path.join(HERE, 'three_compartment_results.json')))
    rows = {}
    for a in ('CD19', 'CD20', 'CD22', 'CD79B'):
        i = ant.index(a)
        mine = {lab: float(np.median([auc(pt['B'][:, i], pt['comp'][lab][:, i]) for pt in pts]))
                for lab in ('T', 'O', 'NB')}
        rows[a] = {'mine': mine, 'published': {k: pub[a][k] for k in ('T', 'O', 'NB')},
                   'max_abs_diff': max(abs(mine[k] - pub[a][k]) for k in ('T', 'O', 'NB'))}
    return rows


def print_summary(r):
    print('=' * 100)
    print(f"CEILING IN CLINICAL UNITS -- DLBCL CyTOF ({r['meta']['n_patients']} patients)")
    print('max % of tumour B cells killable while sparing (1-b) of each normal compartment')
    print('=' * 100)
    for k, v in r['summary'].items():
        print(f"\nbudget {float(k)*100:.0f}%   markers above null envelope: {v['n_above_null']}/{v['n_markers']}")
        print('   strict (all compartments)  approved:', {a: round(x, 3) for a, x in v['approved_all'].items()},
              ' failed:', {a: round(x, 3) for a, x in v['failed_all'].items()})
        print('   non-expendable only        approved:', {a: round(x, 3) for a, x in v['approved_nonexp'].items()},
              ' failed:', {a: round(x, 3) for a, x in v['failed_nonexp'].items()})


if __name__ == '__main__':
    r = compute_all()
    ant, pts = load()
    r['crosscheck_vs_published'] = crosscheck_against_published(ant, pts)
    print_summary(r)
    for a, v in r['crosscheck_vs_published'].items():
        print(f"  crosscheck {a}: max |mine - published| = {v['max_abs_diff']:.3f}")
    json.dump(r, open(OUT, 'w'), indent=1)
    print('\nwrote', OUT)
