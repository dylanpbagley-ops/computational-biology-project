"""Recompute the already published label-flip checks, without changing the panel."""
from pathlib import Path
import sys,json,itertools
import numpy as np
SRC=Path(__file__).resolve().parent.parent/'ev_capture_model/src/aim13_scoping'
sys.path.insert(0,str(SRC))
import label_fragility as LF
genes,score,lab=LF.load(); safe=lab==1
old=json.loads((SRC/'label_fragility.json').read_text())
result={}
for k in range(1,7):
    count=significant=0; amin=1.; pmax=0.
    for idx in itertools.combinations(range(len(genes)),k):
        l=safe.copy(); l[list(idx)]=~l[list(idx)]
        a,p=LF.rank_stats(score,l)
        count+=1; significant+=p<.05; amin=min(amin,a); pmax=max(pmax,p)
    row={'n_sets':count,'frac_p_lt_0.05':significant/count,'sep_min':amin,'worst_p':pmax}
    row['matches_stored']=all(abs(row[x]-old['F1'][str(k)][x])<1e-12 for x in list(row))
    result[str(k)]=row
    print(k,row,flush=True)
weak=[genes.index(g) for g in LF.WEAK_TO_SAFE+LF.WEAK_TO_TOX]
matched=0
for mask in range(256):
    l=safe.copy()
    for b,i in enumerate(weak):
        if mask>>b&1:l[i]=not l[i]
    a,p=LF.rank_stats(score,l)
    r=old['F2_all_subsets'][mask]
    matched+=abs(a-r['sep'])<1e-12 and abs(p-r['p'])<1e-12
result['weak_subsets_matching']=matched
Path(__file__).with_name('fragility_recheck.json').write_text(json.dumps(result,indent=2))
print('Weak subsets matched:',matched,flush=True)
