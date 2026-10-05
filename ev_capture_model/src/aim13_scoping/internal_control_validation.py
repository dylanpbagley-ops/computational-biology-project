"""
THE VALIDATION THAT MATTERS: does the cross-sample RLN normal-B comparator agree with a
SAME-TUBE internal normal control?

For the 18 light-chain-restricted patients that also retain a residual polyclonal
population, the minority light chain is normal B cells from the same patient, same tube,
same run, same staining -- the best possible matched control. If AUC(clone vs internal
normal) agrees with AUC(clone vs batch-matched RLN normal), the RLN comparator is not
producing a technical artifact and the section-3/4 results stand.
"""
import os, glob, json, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
CACHE=os.path.abspath(os.path.join(HERE,'..','..','..','ev_capture_data','gated_cache'))
rng=np.random.default_rng(17); SUB=12000
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
IK,IL=ANT.index('IgKap'),ANT.index('IgLam')
lc=json.load(open(os.path.join(HERE,'light_chain_calls.json')))
usable=lc['usable']; calls=lc['calls']
bat=lambda k:k.split('_')[0]
rlns=sorted(k for k in files if '_RLN_' in k and 'pooled' not in k)
normalB={}
for b in sorted({bat(k) for k in rlns}):
    m=np.concatenate([np.load(files[k],allow_pickle=True)['B'] for k in rlns if bat(k)==b],axis=0)
    normalB[b]=m[rng.choice(len(m),min(len(m),40000),replace=False)].astype(float)

pairs={a:{'internal':[],'rln':[]} for a in ANT}
for key in usable:
    d=np.load(files[key],allow_pickle=True); B=d['B'].astype(float); T=d['T'].astype(float)
    kc,lcut=np.percentile(T[:,IK],99.0),np.percentile(T[:,IL],99.0)
    k,l=B[:,IK],B[:,IL]
    kpos,lpos=(k>kc)&(l<=lcut),(l>lcut)&(k<=kc)
    clone,internal=(kpos,lpos) if calls[key]['call']=='kappa-restricted' else (lpos,kpos)
    C,I=sub(B[clone]),sub(B[internal]); N=sub(normalB[bat(key)])
    for i,a in enumerate(ANT):
        pairs[a]['internal'].append(auc(C[:,i],I[:,i]))
        pairs[a]['rln'].append(auc(C[:,i],N[:,i]))

print('%d light-chain-restricted patients with a same-tube internal normal control\n' % len(usable))
print('%-10s %14s %12s %8s' % ('ANTIGEN','internal(same tube)','RLN(cross)','diff'))
xs,ys=[],[]
for a in sorted(ANT, key=lambda x:-np.median(pairs[x]['internal'])):
    mi,mr=np.median(pairs[a]['internal']),np.median(pairs[a]['rln'])
    xs.append(mi); ys.append(mr)
    print('%-10s %14.3f %12.3f %8.3f' % (a,mi,mr,mr-mi))
xs,ys=np.array(xs),np.array(ys)
r=np.corrcoef(xs,ys)[0,1]
def rk(v): 
    o=np.argsort(np.argsort(v)).astype(float); return (o-o.mean())
sp=float(rk(xs)@rk(ys)/np.sqrt((rk(xs)@rk(xs))*(rk(ys)@rk(ys))))
print('\nAgreement between same-tube internal control and cross-sample RLN comparator:')
print('  Pearson r  = %.3f' % r)
print('  Spearman   = %.3f' % sp)
print('  median |difference| across 37 antigens = %.3f' % np.median(np.abs(ys-xs)))
print('  mean signed difference (RLN minus internal) = %+.3f' % np.mean(ys-xs))
VAL=['CD19','CD20','CD22','CD79B']
print('\nThe four approved-drug targets, same-tube internal control:')
for a in VAL:
    v=np.array(pairs[a]['internal'])
    print('  %-6s internal AUC median %.3f  IQR [%.3f, %.3f]   (RLN comparator gave %.3f)' % (
        a, np.median(v), np.percentile(v,25), np.percentile(v,75), np.median(pairs[a]['rln'])))
json.dump({a:{'internal_median':float(np.median(pairs[a]['internal'])),
              'rln_median':float(np.median(pairs[a]['rln']))} for a in ANT},
          open(os.path.join(HERE,'internal_control_validation.json'),'w'), indent=1)
