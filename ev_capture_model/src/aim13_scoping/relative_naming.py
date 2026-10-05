"""Round 6C naming test (exploratory; same metric as independent_tests.py): does weight x relative-exposure put the DOCUMENTED damaged tissue near the top?"""
import json, os
import numpy as np
from scipy.stats import rankdata
import target_labels_validation as V
import expendability_weights as W
import target_labels as TL
import damage_tissue_map as DM
import wholebody_v2 as WB
import faers_severity as F
import faers_rules as FR
import relative_exposure as RE
from independent_tests import sign_test_p

HERE = os.path.dirname(os.path.abspath(__file__))
WB.GENES.update(t[0] for t in V.VALIDATION)


def main():
    vals, _ = WB.load('rna_single_cell_type.tsv.zip', ['nCPM', 'nTPM', 'pTPM'], ['Cell type'])
    genes = [g for g in DM.DAMAGED if g in vals]
    types = sorted({c for g in genes for c in vals[g]})
    sysof = np.array([FR.SI[F.CELL_SYSTEM.get(c, 'unmapped')] for c in types])
    ws = FR.weight_sets()
    E = RE.rel_exposure(vals, genes, types)
    W_ = {'E only': np.ones(len(types)), 'binary': np.array([0.0 if TL.haematopoietic(c) else 1.0 for c in types]),
          'hand': np.array([W.weight(c, 'continuous') for c in types]), 'death': ws['death'][sysof], 'severe': ws['severe'][sysof]}
    res = {'per_target': {}, 'summary': {}}
    for gi, g in enumerate(genes):
        row = {}
        for m, w in W_.items():
            r = rankdata(-(E[gi] * w), method='average')
            row[m] = float(min(r[types.index(c)] for c in DM.DAMAGED[g] if c in types))
        res['per_target'][g] = row
    for m in W_:
        rk = np.array([res['per_target'][g][m] for g in genes])
        res['summary'][m] = {'median_best_rank': float(np.median(rk)), 'in_top3': int((rk <= 3).sum()), 'in_top10': int((rk <= 10).sum())}
        if m != 'E only':
            base = np.array([res['per_target'][g]['E only'] for g in genes])
            better = int((rk < base).sum()); worse = int((rk > base).sum())
            res['summary'][m].update({'better_than_E_only': better, 'worse': worse, 'sign_p': sign_test_p(better, worse)})
    res['n'] = len(genes)
    json.dump(res, open(os.path.join(HERE, 'relative_naming.json'), 'w'), indent=1)
    return res


if __name__ == '__main__':
    r = main()
    print('n', r['n'])
    for m, v in r['summary'].items(): print(m, v)
    for g, row in r['per_target'].items(): print(g, {k: round(x, 1) for k, x in row.items()})
