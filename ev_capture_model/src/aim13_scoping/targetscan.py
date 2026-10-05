"""
targetscan -- command-line front end to this project's target-safety screen.

    python targetscan.py CD19 MS4A1 CD47 ERBB2 CLDN18 MSLN
    python targetscan.py --top 25            # the 25 most lineage-confined human genes
    python targetscan.py --ceiling 4.2 0.8 2.1 0.9    # tumour mu/sigma, normal mu/sigma (log10)

What it reports per gene:
  confinement   share of the gene's total normal expression in its single highest cell type
  cell types    how many normal cell types carry it above 10% of that maximum
  home          the highest-expressing normal cell type
  expendable    whether that home is a haematopoietic lineage (regenerated from marrow, with
                an established clinical replacement)
  verdict       the two-factor rule from section 23.2, and where the gene falls relative to
                the 19 clinically safe and 10 on-target-toxic targets in the frozen panel

The ceiling mode implements the ROC-ceiling theorem: given tumour and normal antigen
distributions, no kill function that is monotone in antigen can beat the ROC convex hull, so
the best achievable (tumour killed, normal spared) trade-off is computable in advance.
"""
import os, sys, json, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WEB = json.load(open(os.path.join(HERE, 'confinement_web.json'), encoding='utf-8'))
CT, HAEM = WEB['cell_types'], WEB['haem']
GENES = {r[0]: r for r in WEB['genes']}
RANK = {r[0]: i + 1 for i, r in enumerate(WEB['genes'])}
N = len(WEB['genes'])

_v2 = os.path.join(HERE, 'wholebody_v2.json')
PANEL = json.load(open(_v2, encoding='utf-8')) if os.path.exists(_v2) else None


def score(rec):
    """The two-factor rule: expendable home first, confinement second."""
    return rec[1] + (1.0 if HAEM[rec[3]] else 0.0)


def describe(name):
    rec = GENES.get(name) or GENES.get(name.upper())
    if rec is None:
        return '%-10s  not found in the reference (check the HGNC symbol)' % name
    g, share, ntypes, j = rec
    home = CT[j]
    exp = HAEM[j]
    out = ['%s' % g,
           '   confinement      %.3f      (rank %d of %d human genes)' % (share, RANK[g], N),
           '   cell types       %d' % ntypes,
           '   home             %s' % home,
           '   expendable home  %s' % ('YES - blood lineage, clinically replaceable' if exp
                                       else 'no  - solid tissue')]
    if PANEL:
        lab = PANEL['labels']
        cA = PANEL['reference_A_cell_type']
        safe = [cA[k]['top_share'] + (1.0 if PANEL['haem_home_A'][k] else 0.0)
                for k in cA if lab.get(k) == 'safe']
        tox = [cA[k]['top_share'] + (1.0 if PANEL['haem_home_A'][k] else 0.0)
               for k in cA if lab.get(k) == 'toxic']
        s = score(rec)
        pc_safe = 100.0 * sum(1 for v in safe if v < s) / len(safe)
        pc_tox = 100.0 * sum(1 for v in tox if v < s) / len(tox)
        out.append('   two-factor score %.3f' % s)
        out.append('   scores above %.0f%% of the 19 clinically safe targets and %.0f%% of the '
                   '10 on-target-toxic ones' % (pc_safe, pc_tox))
        if s >= 1.0:
            out.append('   READ: confined to an expendable lineage -- the pattern every approved '
                       'blood-cancer target follows.')
        elif share >= 0.25:
            out.append('   READ: confined, but to a tissue you cannot lose. This is the CEACAM5 '
                       'pattern: narrow is not the same as safe.')
        else:
            out.append('   READ: broadly expressed across normal tissue. This is the HER2 / CD47 '
                       'pattern, where on-target toxicity has repeatedly been dose-limiting.')
    return '\n'.join(out)


def ceiling(mu_t, sd_t, mu_n, sd_n, n=200000, seed=0):
    """ROC ceiling for lognormal tumour and normal antigen distributions (log10 mu/sigma)."""
    rng = np.random.default_rng(seed)
    t = rng.normal(mu_t, sd_t, n)
    m = rng.normal(mu_n, sd_n, n)
    thr = np.quantile(np.concatenate([t, m]), np.linspace(0, 1, 2001))
    tpr = np.array([(t >= x).mean() for x in thr])
    fpr = np.array([(m >= x).mean() for x in thr])
    J = tpr - fpr
    k = int(np.argmax(J))
    # AUC = P(tumour antigen > normal antigen), computed exactly by rank statistic
    order = np.argsort(np.concatenate([t, m]), kind='mergesort')
    ranks = np.empty(2 * n, float)
    ranks[order] = np.arange(1, 2 * n + 1)
    auc = float((ranks[:n].sum() - n * (n + 1) / 2.0) / (n * n))
    return {'best_threshold_log10': float(thr[k]), 'tumour_killed': float(tpr[k]),
            'normal_killed': float(fpr[k]), 'youden_J': float(J[k]), 'auc_approx': auc}


def main(argv):
    if not argv:
        print(__doc__)
        return
    if argv[0] == '--top':
        k = int(argv[1]) if len(argv) > 1 else 20
        print('%-10s %10s %8s  %-34s %s' % ('GENE', 'confine', 'n_types', 'home', 'expendable'))
        for g, share, ntypes, j in WEB['genes'][:k]:
            print('%-10s %10.3f %8d  %-34s %s'
                  % (g, share, ntypes, CT[j][:34], 'YES' if HAEM[j] else 'no'))
        return
    if argv[0] == '--ceiling':
        mu_t, sd_t, mu_n, sd_n = (float(x) for x in argv[1:5])
        r = ceiling(mu_t, sd_t, mu_n, sd_n)
        print('ROC ceiling for these two antigen distributions')
        print('  no monotone kill function can beat this trade-off:')
        print('    best single threshold   10^%.2f molecules/cell' % r['best_threshold_log10'])
        print('    tumour cells killed     %.1f%%' % (100 * r['tumour_killed']))
        print('    normal cells killed     %.1f%%' % (100 * r['normal_killed']))
        print("    Youden's J              %.3f" % r['youden_J'])
        print('    separation (AUC)        %.3f' % r['auc_approx'])
        return
    for a in argv:
        print(describe(a))
        print()


if __name__ == '__main__':
    main(sys.argv[1:])
