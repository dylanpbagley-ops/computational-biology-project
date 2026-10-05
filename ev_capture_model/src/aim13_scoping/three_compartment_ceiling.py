"""
Therapeutic-index ceiling against THREE normal compartments:
  T  = bystander T cells, same file
  O  = other CD45+ leukocytes that are neither B nor T (myeloid, NK, etc.), same file
  NB = normal B cells from batch-matched reactive lymph node controls
A target must spare ALL of them, so the ceiling is the minimum across compartments.
Adding a compartment can only LOWER a ceiling, never raise it -- so every number here
is a valid upper bound and every exclusion is permanent.
"""
import os, glob, json, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
CACHE=os.path.abspath(os.path.join(HERE,'..','..','..','ev_capture_data','gated_cache'))
rng=np.random.default_rng(31); SUB=12000
def auc(x,y):
    n1,n2=len(x),len(y); z=np.concatenate([x,y]); o=np.argsort(z,kind='mergesort'); zs=z[o]
    r=np.empty(len(z)); i=0
    while i<len(z):
        j=i
        while j+1<len(z) and zs[j+1]==zs[i]: j+=1
        r[i:j+1]=(i+j)/2.0+1; i=j+1
    rr=np.empty(len(z)); rr[o]=r
    return (rr[:n1].sum()-n1*(n1+1)/2.0)/(n1*n2)
def sub(a): return a[rng.choice(len(a),SUB,replace=False)] if len(a)>SUB else a
files={os.path.basename(p).replace('.npz',''):p for p in glob.glob(os.path.join(CACHE,'*.npz'))}
ANT=list(np.load(next(iter(files.values())),allow_pickle=True)['antigens'])
bat=lambda k:k.split('_')[0]
tumors=sorted(k for k in files if '_DLB_' in k)
rlns=sorted(k for k in files if '_RLN_' in k and 'pooled' not in k)
normalB={b: np.concatenate([np.load(files[k],allow_pickle=True)['B'] for k in rlns if bat(k)==b],axis=0).astype(float)
         for b in sorted({bat(k) for k in rlns})}
res={a:{'T':[],'O':[],'NB':[]} for a in ANT}
for t in tumors:
    d=np.load(files[t],allow_pickle=True)
    B,T,O=d['B'].astype(float),d['T'].astype(float),d['O'].astype(float)
    if len(B)<500 or len(T)<500 or len(O)<500: print('skip',t); continue
    Bs,Ts,Os,Ns=sub(B),sub(T),sub(O),sub(normalB[bat(t)])
    for i,a in enumerate(ANT):
        res[a]['T'].append(auc(Bs[:,i],Ts[:,i]))
        res[a]['O'].append(auc(Bs[:,i],Os[:,i]))
        res[a]['NB'].append(auc(Bs[:,i],Ns[:,i]))
summ={a:{k:float(np.median(v)) for k,v in res[a].items()} for a in ANT}
for a in summ: summ[a]['joint3']=min(summ[a]['T'],summ[a]['O'],summ[a]['NB'])
VAL=['CD19','CD20','CD22','CD79B']
print('%-10s %7s %7s %7s %9s %s' % ('ANTIGEN','vs T','vs other','vs normB','JOINT-3','binding compartment'))
for a in sorted(ANT,key=lambda x:-summ[x]['joint3']):
    s=summ[a]; b=min(('T','O','NB'), key=lambda k:s[k])
    lbl={'T':'T cells','O':'myeloid/NK','NB':'normal B'}[b]
    tag='  <-- approved drug' if a in VAL else ''
    print('%-10s %7.3f %7.3f %7.3f %9.3f  %-12s%s' % (a,s['T'],s['O'],s['NB'],s['joint3'],lbl,tag))
print('\nBest 3-compartment ceiling in the panel: %.3f (%s)' % (
    max(summ[a]['joint3'] for a in ANT), max(ANT,key=lambda a:summ[a]['joint3'])))
print('Best among the four approved-drug targets: %.3f' % max(summ[a]['joint3'] for a in VAL))
# Monotonicity check, done EXACTLY: joint2 and joint3 are computed from the SAME
# AUC estimates, so min over a larger set is mathematically forced to be <=.
# (Comparing against a previously saved run instead would mix in resampling jitter.)
for a in ANT: summ[a]['joint2']=min(summ[a]['T'],summ[a]['NB'])
viol=[a for a in ANT if summ[a]['joint3'] > summ[a]['joint2'] + 1e-12]
drop=[(a, summ[a]['joint2'], summ[a]['joint3']) for a in ANT if summ[a]['joint3'] < summ[a]['joint2']-1e-12]
print('\nAdding the third compartment LOWERED %d of %d ceilings; raised %d (must be 0).'
      % (len(drop), len(ANT), len(viol)))
print('Largest drops caused by the new myeloid/NK compartment:')
for a,p,n in sorted(drop,key=lambda r:r[1]-r[2],reverse=True)[:8]:
    print('   %-10s %.3f -> %.3f  (%+.3f)' % (a,p,n,n-p))
# separately, compare to the earlier saved 2-compartment run to size the resampling jitter
prev=json.load(open(os.path.join(HERE,'joint_window_results.json')))['joint']
jit=np.array([summ[a]['joint2']-prev[a] for a in ANT])
print('\nRun-to-run jitter of the 2-compartment ceiling (different random subsample draws):')
print('   median %+.4f, max |jitter| %.4f  -- this is why the comparison above is done within a run'
      % (np.median(jit), np.abs(jit).max()))
json.dump(summ, open(os.path.join(HERE,'three_compartment_results.json'),'w'), indent=1)
