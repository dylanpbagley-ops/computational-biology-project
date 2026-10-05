"""POST-HOC addendum to validation_test.py (recorded as post-hoc): two validation genes (CTAG1B = NY-ESO-1, MAGEA12) have NO measured expression in any of the
reference cell types, so a share-based score is undefined and they were dropped from the pre-specified test. CTAG1B was a primary SAFE target.
Convention applied after seeing that it was undefined (not after seeing any ranking): a gene absent from all normal cell types is the safest possible, so it is scored above every other gene."""
import json, os
import numpy as np
import wholebody_v2 as WB
HERE = os.path.dirname(os.path.abspath(__file__))
r = json.load(open(os.path.join(HERE, 'validation_test.json')))
t = r['targets']
safe = [g for g, v in t.items() if v['tier'] == 'primary' and v['outcome'] == 'safe']
tox = [g for g, v in t.items() if v['tier'] == 'primary' and v['outcome'] == 'toxic']
top = max(v['rule'] for v in t.values()) + 1.0
a = [t[g]['rule'] for g in safe] + [top]          # CTAG1B, absent from normal cells
b = [t[g]['rule'] for g in tox]
p, exact = WB.exact_p(a, b)
out = {'n_safe': len(a), 'n_toxic': len(b), 'sep_with_absent_gene_as_safest': WB.sep(a, b), 'p_exact': float(p),
       'note': 'post-hoc convention; CTAG1B has no measured expression in the reference'}
json.dump(out, open(os.path.join(HERE, 'validation_addendum.json'), 'w'), indent=1)
print(out)
