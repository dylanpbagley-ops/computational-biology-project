"""
DO THE FOUR NORMAL REFERENCES EVEN AGREE WITH EACH OTHER?

The protein cross-check split: the finding replicates on immunohistochemistry (0.785,
p = 0.010) and does not replicate on mass spectrometry (0.650, p = 0.115). Coarseness was
tested as an explanation and REJECTED -- randomly pooling the 154-cell-type RNA reference down
to 20 groups still gives 0.837.

So the question becomes factual rather than rhetorical: how much do these references agree
about how confined a given protein is? This script measures that directly, with no
interpretation attached.
"""
import os, json, itertools
import numpy as np
from scipy.stats import spearmanr
import target_labels as TL

HERE = os.path.dirname(os.path.abspath(__file__))
LABEL = {t[0]: t[3] for t in TL.TARGETS}

W = json.load(open(os.path.join(HERE, 'wholebody_v2.json')))
P = json.load(open(os.path.join(HERE, 'protein_cross_check.json')))

refs = {'RNA single-cell (154)': {g: d['top_share'] for g, d in W['reference_A_cell_type'].items()},
        'RNA bulk tissue (51)': {g: d['top_share'] for g, d in W['reference_B_tissue'].items()}}
for name, r in P.items():
    short = ('PROTEIN IHC (76)' if 'immunohistochemistry, tissue' in name else
             'PROTEIN IHC strict' if 'uncertain' in name else
             'PROTEIN MS samples (20)' if '60 samples' in name else
             'PROTEIN MS tissues (20)')
    refs[short] = {g: d['top_share'] for g, d in r['per_gene'].items()}

names = list(refs)
print('pairwise Spearman agreement on confinement, over the targets each pair shares:\n')
print('  %-24s %-24s %6s %8s %6s' % ('reference A', 'reference B', 'rho', 'p', 'n'))
out = {}
for a, b in itertools.combinations(names, 2):
    shared = sorted(set(refs[a]) & set(refs[b]))
    if len(shared) < 6:
        continue
    x = [refs[a][g] for g in shared]
    y = [refs[b][g] for g in shared]
    rho, p = spearmanr(x, y)
    out['%s vs %s' % (a, b)] = {'rho': float(rho), 'p': float(p), 'n': len(shared)}
    print('  %-24s %-24s %6.3f %8.4f %6d' % (a, b, rho, p, len(shared)))

print("""
The row that matters is RNA single-cell versus PROTEIN MS. If those two references disagree
about which proteins are confined, then they are not two measurements of one quantity, and a
result holding in one while failing in the other is not a contradiction -- it is two different
quantities. Whether the framework's quantity is the useful one is then an open question, and
the write-up has to say so.""")

json.dump(out, open(os.path.join(HERE, 'reference_agreement.json'), 'w'), indent=1)
print('\nwrote reference_agreement.json')
