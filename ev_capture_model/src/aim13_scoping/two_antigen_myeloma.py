"""Two-antigen extension, replicated in multiple myeloma (GSE234261 CITE-seq, 63 proteins).

Same hierarchy, same held-out protocol as two_antigen_ceiling.py. Different disease,
different lab, different technology. Normal = T + myeloid + NK + normal B pooled.
"""
import os, sys, glob, json, itertools, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from two_antigen_ceiling import gates_and_lr, edges
ADT=os.path.abspath(os.path.join(HERE,'..','..','..','ev_capture_data','myeloma_adt'))
rng=np.random.default_rng(211); SUB=4000
def split(a):
    a=np.asarray(a,float)
    if len(a)>2*SUB: a=a[rng.choice(len(a),2*SUB,replace=False)]
    p=rng.permutation(len(a)); h=len(a)//2
    return a[p[:h]], a[p[h:]]
def clr(M):
    X=np.log(M+1.0); return X-X.mean(axis=1,keepdims=True)

files=sorted(glob.glob(os.path.join(ADT,'*.npz')))
ALL=[str(x) for x in np.load(files[0],allow_pickle=True)['names']]
ix={n:i for i,n in enumerate(ALL)}
# CIRCULARITY FIX: these 7 antigens DEFINE the gates -- CD3E/CD14/CD33/CD56/CD19 mark the
# normal compartments and CD138/CD269 mark the tumor -- so any pair containing them separates
# the two populations by construction. They are excluded from the pair set entirely.
GATING={'CD3E','CD14','CD33','CD56','CD19','CD138','CD269'}
NAMES=[n for n in ALL if n not in GATING]
print('excluding %d gating markers; %d antigens remain, %d pairs'
      % (len(GATING), len(NAMES), len(NAMES)*(len(NAMES)-1)//2))
acc={}
nused=0
for f in files:
    M=clr(np.load(f,allow_pickle=True)['adt'].astype(float))
    hi=lambda v,p: v>np.percentile(v,p)
    cd3,cd14,cd33,cd56,cd19,cd138,bcma=(M[:,ix[k]] for k in ['CD3E','CD14','CD33','CD56','CD19','CD138','CD269'])
    T=hi(cd3,70); MY=(~T)&(hi(cd14,80)|hi(cd33,80)); NK=(~T)&(~MY)&hi(cd56,85)
    NB=(~T)&(~MY)&(~NK)&hi(cd19,85)&(cd138<=np.percentile(cd138,70))
    TU=(~T)&(~MY)&(~NK)&(~NB)&(hi(cd138,70)|hi(bcma,80))
    if TU.sum()<300: continue
    NORM=M[T|MY|NK|NB]
    tr_t,te_t=split(M[TU]); tr_n,te_n=split(NORM)
    tr_t=tr_t[:,[ix[n] for n in NAMES]]; te_t=te_t[:,[ix[n] for n in NAMES]]
    tr_n=tr_n[:,[ix[n] for n in NAMES]]; te_n=te_n[:,[ix[n] for n in NAMES]]
    E=[edges(np.concatenate([tr_t[:,c],tr_n[:,c]])) for c in range(len(NAMES))]
    for i,j in itertools.combinations(range(len(NAMES)),2):
        r=gates_and_lr(tr_t,tr_n,te_t,te_n,i,j,E[i],E[j])
        acc.setdefault((NAMES[i],NAMES[j]),[]).append((r['AND'],r['NOT'],r['LR'],r['SINGLE_in_pair'],r['AND_train'],r['NOT_train'],r['LR_train'],r['SINGLE_train'],r['SINGLEHI'],r['SINGLEHI_train']))
    nused+=1
    print('  processed %s' % os.path.basename(f), flush=True)
rows=[]
for k,v in acc.items():
    v=np.array(v)
    rows.append({'pair':list(k),'AND':float(np.median(v[:,0])),'NOT':float(np.median(v[:,1])),
                 'LR':float(np.median(v[:,2])),'SINGLE':float(np.median(v[:,3])),
                 'AND_tr':float(np.median(v[:,4])),'NOT_tr':float(np.median(v[:,5])),
                 'LR_tr':float(np.median(v[:,6])),'SINGLE_tr':float(np.median(v[:,7])),
                 'SINGLEHI':float(np.median(v[:,8])),'SINGLEHI_tr':float(np.median(v[:,9]))})
json.dump({'n_samples':nused,'pairs':rows}, open(os.path.join(HERE,'two_antigen_myeloma.json'),'w'))
ghi=np.array([max(r['AND'],r['NOT'])-r['SINGLEHI'] for r in rows])
g=np.array([max(r['AND'],r['NOT'])-r['SINGLE'] for r in rows])
print('\nMYELOMA REPLICATION: %d samples, %d antigen pairs' % (nused,len(rows)))
print('  ENGINEERING gain (2 receptors vs 1): %d of %d (%.0f%%), median %+.3f, max %+.3f'
      % ((ghi>0).sum(),len(ghi),100*(ghi>0).mean(),np.median(ghi),ghi.max()))
print('  STATISTICAL gain (vs 1 receptor allowed to kill the low side): %d of %d (%.0f%%), median %+.3f, max %+.3f'
      % ((g>0).sum(),len(g),100*(g>0).mean(),np.median(g),g.max()))
bad=sum(1 for r in rows if r['LR_tr']<max(r['AND_tr'],r['NOT_tr'],r['SINGLE_tr'])-1e-9)
print('  IN-SAMPLE theorem violations (must be 0): %d of %d' % (bad,len(rows)))
print('  overfit gap, LR: %.3f ; best gate: %.3f' % (
      np.median([r['LR_tr']-r['LR'] for r in rows]), np.median([max(r['AND_tr'],r['NOT_tr'])-max(r['AND'],r['NOT']) for r in rows])))
nots=sum(1 for r in rows if r['NOT']>r['AND']+0.02)
print('  pairs where NOT beats AND: %d of %d (%.0f%%)' % (nots,len(rows),100*nots/len(rows)))
print('\nTOP 10 PAIRS BY CONTROLLED GAIN')
print('%-22s %9s %8s %8s %8s' % ('PAIR','1-antigen','AND','NOT','LR'))
for r in sorted(rows,key=lambda r:-(max(r['AND'],r['NOT'])-r['SINGLEHI']))[:10]:
    print('%-22s %9.3f %8.3f %8.3f %+8.3f' % ('+'.join(r['pair']),r['SINGLE'],r['AND'],r['NOT'],max(r['AND'],r['NOT'])-r['SINGLE']))
bc=[r for r in rows if 'CD319' in r['pair']]
if bc:
    b=max(bc,key=lambda r:max(r['AND'],r['NOT']))
    print('\nBest pair containing SLAMF7/CD319 (approved myeloma target, not a gating marker):')
    print('   %s  1-antigen %.3f -> best gate %.3f (gain %+.3f)' % ('+'.join(b['pair']),b['SINGLE'],max(b['AND'],b['NOT']),max(b['AND'],b['NOT'])-b['SINGLE']))
