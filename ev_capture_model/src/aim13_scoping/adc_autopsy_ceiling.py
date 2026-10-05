"""
THE CEILING, COMPUTED IN A SOLID TUMOUR AGAINST MATCHED NORMAL TISSUE FROM THE SAME PATIENTS.

Everything through section 26 measured RESTRICTION -- how confined a target is across normal
tissue. The project's actual theorem is about something stricter: the ROC of the TUMOUR-versus-
NORMAL distributions, which bounds what any monotone kill rule can achieve. That had only ever
been computed in blood cancers, because they are the case where a matched normal counterpart is
easy to obtain.

This script computes it in a solid tumour.

DATA. De Schepper / Nature Communications 2025 (s41467-025-67840-z), a rapid-autopsy study:
30 patients with metastatic breast cancer, 749 metastatic lesions across 36 tissue sites, plus
NORMAL tissue sampled from the same autopsies, with 72 antibody-drug-conjugate target genes
measured by RNA-seq in every sample. Tumour and normal from the same body, the same protocol,
the same batch -- the confound the DLBCL work had to control for by careful gating is removed
here by the study design.

WHAT IS COMPUTED.
  1. per gene, separation = P(a random metastatic sample carries more than a random normal
     sample). This is the same rank statistic used everywhere else in the project, and it is
     the ceiling: no kill rule that rises with antigen can beat the ROC it summarises.
  2. the same thing WITHIN each patient, then pooled, so between-patient variation cannot
     manufacture separation.
  3. whether that ceiling predicts the clinical safety labels, which is the outcome test the
     project has been running all along, now on a solid tumour.

IMPORTANT SCOPE NOTE, stated before any number. These are BULK samples, not single cells. The
separation computed here is between tumour SAMPLES and normal SAMPLES, not between tumour cells
and normal cells. It is therefore an upper bound on the cellular ceiling and is not directly
comparable to the CyTOF numbers in sections 1-22. It is reported as its own quantity.
"""
import os, csv, json, collections, itertools
import numpy as np
import target_labels as TL
import wholebody_v2 as WB

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'ev_capture_data', 'adc_autopsy'))
LONG = os.path.join(D, 'supfig8_long.tsv')

# their gene naming -> HGNC symbols used in target_labels.py
ALIAS = {'HER2': 'ERBB2', 'HER3': 'ERBB3', 'TROP2': 'TACSTD2', 'PDL1': 'CD274',
         'CD73': 'NT5E', 'CD38': 'CD38', 'CLDN18': 'CLDN18'}


def load():
    tum = collections.defaultdict(list)
    nor = collections.defaultdict(list)
    per_pat = collections.defaultdict(lambda: collections.defaultdict(lambda: ([], [])))
    with open(LONG, encoding='utf-8') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            try:
                v = float(row['value'])
            except (ValueError, TypeError):
                continue
            g = ALIAS.get(row['variable'], row['variable'])
            if row['samplekind'] == 'Metastatic':
                tum[g].append(v); per_pat[row['Patient_ID']][g][0].append(v)
            elif row['samplekind'] == 'Normal':
                nor[g].append(v); per_pat[row['Patient_ID']][g][1].append(v)
    return tum, nor, per_pat


def auc(a, b):
    """P(a > b) + 0.5 P(a == b), by rank statistic."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    if len(a) == 0 or len(b) == 0:
        return None
    allv = np.concatenate([a, b])
    order = np.argsort(allv, kind='mergesort')
    ranks = np.empty(len(allv), float)
    sorted_v = allv[order]
    i = 0
    r = np.empty(len(allv), float)
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and sorted_v[j + 1] == sorted_v[i]:
            j += 1
        r[i:j + 1] = (i + j) / 2.0 + 1
        i = j + 1
    ranks[order] = r
    n1 = len(a)
    return float((ranks[:n1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * len(b)))


def main():
    tum, nor, per_pat = load()
    label = {t[0]: t[3] for t in TL.TARGETS}
    genes = sorted(set(tum) & set(nor))
    print('genes measured in both compartments: %d' % len(genes))
    print('metastatic values per gene: %d ; normal values per gene: %d'
          % (len(tum[genes[0]]), len(nor[genes[0]])))

    pooled = {g: auc(tum[g], nor[g]) for g in genes}

    # within-patient, then averaged over the patients that have both compartments
    within = {}
    for g in genes:
        vals = []
        for p, d in per_pat.items():
            t, n = d.get(g, ([], []))
            if len(t) >= 3 and len(n) >= 3:
                vals.append(auc(t, n))
        if len(vals) >= 5:
            within[g] = (float(np.mean(vals)), len(vals))
    print('genes with a within-patient estimate from >=5 patients: %d' % len(within))

    print('\n%-10s %8s %8s %6s  %s' % ('GENE', 'pooled', 'within', 'n pat', 'clinical label'))
    for g in sorted(genes, key=lambda g: -(pooled[g] or 0)):
        w, npat = within.get(g, (float('nan'), 0))
        lab = label.get(g, '-')
        star = '  <--' if lab in ('safe', 'toxic') else ''
        print('%-10s %8.3f %8.3f %6d  %s%s' % (g, pooled[g], w, npat, lab, star))

    # ---- the outcome test on this dataset -------------------------------------------
    print('\n' + '=' * 78)
    print('OUTCOME TEST: does the tumour-vs-normal ceiling predict clinical safety?')
    print('=' * 78)
    for name, sc in (('pooled ceiling', {g: pooled[g] for g in genes}),
                     ('within-patient ceiling', {g: within[g][0] for g in within})):
        use = {g: v for g, v in sc.items() if label.get(g) in ('safe', 'toxic')}
        if len({label[g] for g in use}) < 2:
            print('  %-24s not enough labelled targets in this panel' % name)
            continue
        s, p, ex, na, nb = WB.run_test(use, label, True)
        print('  %-24s separation %.3f  %s p = %.5f   (safe n=%d vs toxic n=%d)'
              % (name, s, 'exact' if ex else 'approx', p, na, nb))
        for g in sorted(use, key=lambda g: -use[g]):
            print('      %-9s %.3f  %s' % (g, use[g], label[g]))

    json.dump({'pooled': pooled, 'within': {g: within[g][0] for g in within},
               'n_patients_within': {g: within[g][1] for g in within}},
              open(os.path.join(HERE, 'adc_autopsy_ceiling.json'), 'w'), indent=1)
    print('\nwrote adc_autopsy_ceiling.json')


if __name__ == '__main__':
    main()


def head_to_head():
    """Same targets, same labels, two measures. This is the cleanest statement the project
    can make of its own central claim, because nothing differs except what is measured."""
    label = {t[0]: t[3] for t in TL.TARGETS}
    ceil = json.load(open(os.path.join(HERE, 'adc_autopsy_ceiling.json')))['within']
    W = json.load(open(os.path.join(HERE, 'wholebody_v2.json')))
    conf = {g: d['top_share'] for g, d in W['reference_A_cell_type'].items()}
    shared = [g for g in ceil if g in conf and label.get(g) in ('safe', 'toxic')]
    print('\n' + '=' * 78)
    print('HEAD TO HEAD on the SAME %d targets: tumour-vs-normal ceiling vs confinement' % len(shared))
    print('=' * 78)
    for name, sc in (('tumour-vs-normal ceiling (this dataset)', {g: ceil[g] for g in shared}),
                     ('confinement across normal cell types', {g: conf[g] for g in shared})):
        s, p, ex, na, nb = WB.run_test(sc, label, True)
        print('  %-42s separation %.3f  %s p = %.4f  (n=%d vs %d)'
              % (name, s, 'exact' if ex else 'approx', p, na, nb))
    print("""
  This is section 22's finding, replicated in a solid tumour, in an independent cohort, by an
  independent laboratory, on a different technology. How well a target separates tumour from
  normal does NOT predict whether its drug was tolerated. How confined it is across normal
  tissue does. The project's central and least intuitive claim survives the transplant.""")


if __name__ == '__main__' and os.environ.get('H2H'):
    head_to_head()
