"""
Step 9: the therapeutic-index ceiling, replicated in MULTIPLE MYELOMA.

Different cancer, different lab, different technology (CITE-seq vs mass cytometry),
different targets. If the framework is measuring biology, it should behave consistently.

Normalisation: CLR (centred log ratio) per cell -- the standard for CITE-seq antibody
counts. Needed because cells differ in sequencing depth; without it a deeply-sequenced
cell looks high for every protein at once.

Populations (bone marrow):
  tumor   = malignant plasma cells: CD138+ and/or BCMA(CD269)+, CD3E-, CD14-, CD19-
  T       = CD3E+
  myeloid = CD14+ or CD33+ (monocytes/granulocytes), CD3E-
  NK      = CD56+ CD3E-
  normalB = CD19+ CD138-
A target must spare ALL normal compartments -> ceiling = min across them.
"""
import os, glob, json, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
ADT=os.path.abspath(os.path.join(HERE,'..','..','..','ev_capture_data','myeloma_adt'))
rng=np.random.default_rng(53); SUB=8000

def auc(x,y):
    n1,n2=len(x),len(y)
    if n1<50 or n2<50: return np.nan
    z=np.concatenate([x,y]); o=np.argsort(z,kind='mergesort'); zs=z[o]
    r=np.empty(len(z)); i=0
    while i<len(z):
        j=i
        while j+1<len(z) and zs[j+1]==zs[i]: j+=1
        r[i:j+1]=(i+j)/2.0+1; i=j+1
    rr=np.empty(len(z)); rr[o]=r
    return (rr[:n1].sum()-n1*(n1+1)/2.0)/(n1*n2)
def sub(a): return a[rng.choice(len(a),SUB,replace=False)] if len(a)>SUB else a
def clr(M):
    X=M+1.0
    return np.log(X) - np.log(X).mean(axis=1, keepdims=True)

files=sorted(glob.glob(os.path.join(ADT,'*.npz')))
NAMES=list(np.load(files[0],allow_pickle=True)['names'])
ix={n:i for i,n in enumerate(NAMES)}
res={n:{'T':[],'MY':[],'NK':[],'NB':[]} for n in NAMES}
pop_report=[]
for f in files:
    sid=os.path.basename(f)[:-4]
    M=clr(np.load(f,allow_pickle=True)['adt'].astype(float))
    q=lambda m,p: np.percentile(M[:,ix[m]],p)
    cd3,cd14,cd33,cd56,cd19,cd138,bcma = (M[:,ix[k]] for k in ['CD3E','CD14','CD33','CD56','CD19','CD138','CD269'])
    hi=lambda v,p=75: v > np.percentile(v,p)
    T  = hi(cd3,70)
    MY = (~T) & (hi(cd14,80) | hi(cd33,80))
    NK = (~T) & (~MY) & hi(cd56,85)
    NB = (~T) & (~MY) & (~NK) & hi(cd19,85) & (cd138<=np.percentile(cd138,70))
    TU = (~T) & (~MY) & (~NK) & (~NB) & (hi(cd138,70) | hi(bcma,80))
    pop_report.append((sid,len(M),int(TU.sum()),int(T.sum()),int(MY.sum()),int(NK.sum()),int(NB.sum())))
    if TU.sum()<100: print('  skip',sid,'too few tumor cells'); continue
    tu=sub(M[TU])
    for lab,mask in (('T',T),('MY',MY),('NK',NK),('NB',NB)):
        if mask.sum()<100: continue
        cm=sub(M[mask])
        for n in NAMES:
            a=auc(tu[:,ix[n]],cm[:,ix[n]])
            if not np.isnan(a): res[n][lab].append(a)
print('%-6s %7s %7s %7s %7s %6s %7s' % ('sample','cells','tumor','T','myeloid','NK','normalB'))
for r in pop_report: print('%-6s %7d %7d %7d %7d %6d %7d' % r)

summ={}
for n in NAMES:
    s={k:(float(np.median(v)) if v else np.nan) for k,v in res[n].items()}
    vals=[v for v in s.values() if not np.isnan(v)]
    s['ceiling']=float(min(vals)) if vals else np.nan
    summ[n]=s
APPROVED=['CD269','CD38','CD319']          # BCMA (ide-cel/cilta-cel/teclistamab), CD38 (daratumumab), SLAMF7 (elotuzumab)
TRIALS  =['CD307']                          # FcRH5 (cevostamab)
FAILED  =['CD138']                          # indatuximab ravtansine -- discontinued
print('\n%-9s %7s %7s %7s %7s %9s  %s' % ('PROTEIN','vs T','vs myel','vs NK','vs normB','CEILING','note'))
for n in sorted(NAMES,key=lambda x:-(summ[x]['ceiling'] if not np.isnan(summ[x]['ceiling']) else -9)):
    s=summ[n]
    note=''
    if n in APPROVED: note='  <-- APPROVED myeloma drug'
    elif n in TRIALS: note='  <-- in clinical trials'
    elif n in FAILED: note='  <-- tried and discontinued'
    f=lambda v: ('%7.3f'%v) if not np.isnan(v) else '      -'
    print('%-9s %s %s %s %s %9s%s' % (n,f(s['T']),f(s['MY']),f(s['NK']),f(s['NB']),
          ('%.3f'%s['ceiling']) if not np.isnan(s['ceiling']) else '-', note))
json.dump(summ, open(os.path.join(HERE,'myeloma_ceiling_results.json'),'w'), indent=1)
