"""
Step 10: does the result depend on how cells were sorted?

The main analysis uses per-file PERCENTILE gates. This re-runs the whole
tumor-vs-normal ceiling with two structurally different, automatic gating methods
on arcsinh-transformed signal (the standard CyTOF transform):

  otsu   -- Otsu's method: pick the cutoff that best separates the marker's own
            histogram into two groups (maximises between-group variance). No
            percentile assumption at all.
  gmm    -- a 2-component 1-D Gaussian mixture fitted by EM; the cutoff is where
            the two fitted components are equally likely. Assumes a negative and a
            positive population, and finds them.

If the antigen ranking and the four approved-drug targets' scores are stable across
all three methods, the conclusions are not an artifact of the gating choice.
"""
import os, glob, json, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
CACHE=os.path.abspath(os.path.join(HERE,'..','..','..','ev_capture_data','gated_cache'))
RAW=os.path.abspath(os.path.join(HERE,'..','..','..','ev_capture_data'))
import sys; sys.path.insert(0, os.path.abspath(os.path.join(HERE,'..')))
import aim12_single_cell_cytof_reordering as a12
rng=np.random.default_rng(41); SUB=10000
asinh=lambda x: np.arcsinh(np.asarray(x,dtype=float)/5.0)   # standard CyTOF cofactor

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

def otsu(v, bins=256):
    v=asinh(v); h,e=np.histogram(v,bins=bins); h=h.astype(float)/h.sum()
    c=(e[:-1]+e[1:])/2
    w0=np.cumsum(h); w1=1-w0
    m0=np.cumsum(h*c)/np.maximum(w0,1e-12); m1=(np.sum(h*c)-np.cumsum(h*c))/np.maximum(w1,1e-12)
    var=w0*w1*(m0-m1)**2
    return float(c[int(np.nanargmax(var))])

def gmm2(v, iters=60):
    v=asinh(v); v=v[np.isfinite(v)]
    lo,hi=np.percentile(v,[20,80]); mu=np.array([lo,hi]); sd=np.array([v.std()/2+1e-6]*2); pi=np.array([.5,.5])
    for _ in range(iters):
        p=np.stack([pi[k]*np.exp(-0.5*((v-mu[k])/sd[k])**2)/sd[k] for k in range(2)])
        p=p/np.maximum(p.sum(0),1e-300)
        n=p.sum(1)+1e-12
        mu=(p*v).sum(1)/n; sd=np.sqrt(np.maximum((p*(v-mu[:,None])**2).sum(1)/n,1e-8)); pi=n/len(v)
    order=np.argsort(mu); mu,sd,pi=mu[order],sd[order],pi[order]
    xs=np.linspace(mu[0],mu[1],400)
    d=[pi[k]*np.exp(-0.5*((xs-mu[k])/sd[k])**2)/sd[k] for k in range(2)]
    return float(xs[int(np.argmin(np.abs(d[0]-d[1])))])

def gate_file(path, method):
    raw=open(path,'rb').read(); h=a12.parse_fcs3_header(raw); tk,st=a12.parse_fcs3_text_segment(raw,h)
    tot,npar=int(tk['TOT']),int(tk['PAR'])
    arr=np.frombuffer(raw,dtype='>f4',count=tot*npar,offset=h['data_start']).reshape(tot,npar)
    g=lambda m: arr[:,a12.find_parameter_index_by_stain_substring(st,m)].astype(np.float32)
    dna,via,bd,cd45=g('193Ir-DNA'),g('195Pt-Viability'),g('140Ce-Beads'),g('CD45')
    cd3,cd19,cd20=g('CD3'),g('CD19'),g('CD20')
    qc=((dna>np.percentile(dna,5))&(dna<np.percentile(dna,99.5))&(via<np.percentile(via,95))
        &(bd<np.percentile(bd,99))&(cd45>np.percentile(cd45,2)))
    if method=='percentile':
        # the main analysis's own rule, evaluated on the SAME files as the others
        c3,c19,c20 = (asinh(np.percentile(cd3[qc],50)), asinh(np.percentile(cd19[qc],50)),
                      asinh(np.percentile(cd20[qc],50)))
    else:
        f = otsu if method=='otsu' else gmm2
        c3,c19,c20 = f(cd3[qc]), f(cd19[qc]), f(cd20[qc])
    a3,a19,a20 = asinh(cd3), asinh(cd19), asinh(cd20)
    B=qc&((a19>c19)|(a20>c20))&(a3<=c3); T=qc&(a3>c3)&(a19<=c19)&(a20<=c20)
    return arr, st, B, T

ANT=list(np.load(glob.glob(os.path.join(CACHE,'*.npz'))[0],allow_pickle=True)['antigens'])
tum=sorted(glob.glob(os.path.join(RAW,'raw_cytof','*.fcs')))[:12]
rln=sorted(glob.glob(os.path.join(RAW,'raw_cytof_rln','*.fcs')))
rln=[p for p in rln if 'pooled' not in p][:6]
VAL=['CD19','CD20','CD22','CD79B']
out={}
for method in ['percentile','otsu','gmm']:
    nb=[]
    for p in rln:
        arr,st,B,T=gate_file(p,method)
        cols=np.stack([arr[:,a12.find_parameter_index_by_stain_substring(st,m)].astype(float) for m in ANT],1)
        nb.append(sub(cols[B]))
    NB=np.concatenate(nb,0); NB=sub(NB)
    acc={a:[] for a in ANT}
    for p in tum:
        arr,st,B,T=gate_file(p,method)
        if B.sum()<500 or T.sum()<500: continue
        cols=np.stack([arr[:,a12.find_parameter_index_by_stain_substring(st,m)].astype(float) for m in ANT],1)
        Bs,Ts=sub(cols[B]),sub(cols[T])
        for i,a in enumerate(ANT):
            acc[a].append(min(auc(Bs[:,i],Ts[:,i]), auc(Bs[:,i],NB[:,i])))
    out[method]={a: float(np.median(acc[a])) for a in ANT}
    print('%s gating done (%d tumors, %d RLN)' % (method, len(tum), len(rln)), flush=True)

pct=json.load(open(os.path.join(HERE,'three_compartment_results.json')))
out['percentile_fullcohort']={a: min(pct[a]['T'],pct[a]['NB']) for a in ANT}
def rk(d):
    o=sorted(ANT,key=lambda a:-d[a]); return {a:i for i,a in enumerate(o)}
def spear(d1,d2):
    r1,r2=rk(d1),rk(d2)
    x=np.array([r1[a] for a in ANT],float); y=np.array([r2[a] for a in ANT],float)
    x-=x.mean(); y-=y.mean(); return float(x@y/np.sqrt((x@x)*(y@y)))
print('\nRank agreement between gating methods (Spearman):')
print('  percentile vs Otsu : %.3f' % spear(out['percentile'],out['otsu']))
print('  percentile vs GMM  : %.3f' % spear(out['percentile'],out['gmm']))
print('  Otsu vs GMM        : %.3f' % spear(out['otsu'],out['gmm']))
print('\nThe four approved-drug targets, joint ceiling under each method:')
print('%-8s %12s %8s %8s' % ('','percentile','otsu','gmm'))
for a in VAL:
    print('%-8s %12.3f %8.3f %8.3f' % (a,out['percentile'][a],out['otsu'][a],out['gmm'][a]))
print('\n  percentile(same files) vs percentile(full cohort): %.3f'
      % spear(out['percentile'],out['percentile_fullcohort']))
print('\nBest target in the panel under each method:')
for m in ['percentile','otsu','gmm']:
    b=max(ANT,key=lambda a:out[m][a]); print('  %-11s %s (%.3f)' % (m,b,out[m][b]))
json.dump(out, open(os.path.join(HERE,'alt_gating_results.json'),'w'), indent=1)
