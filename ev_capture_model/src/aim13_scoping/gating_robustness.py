import sys, os, glob, json, numpy as np
SRC=os.path.expanduser('~/mnt/outputs/ev_capture_model/src'); sys.path.insert(0,SRC)
import aim12_single_cell_cytof_reordering as a12
RAW=os.path.expanduser('~/mnt/outputs/ev_capture_data/raw_cytof')
rng=np.random.default_rng(1)
ANT=['CD19','CD20','CD22','CD79B','CD38','CD40','CD24','CD21','CD72','CD48','CD84','CD43','CD44',
     'CD27','CD25','CD23','CD80','CD86','CD83','CXCR4','CXCR5','PD-L1','FasR','TRAIL-R1','IL-4R',
     'CD49D','CD10','HLA-DR','HLA-ABC','HLA-E','CCR4','BTLA','CD45RA','IgM','IgD','IgKap','IgLam']
VAL=['CD19','CD20','CD22','CD79B']
BROAD=['CD44','HLA-ABC','CD48','CD27','CD43','CD38','CD84','HLA-E','CD25','FasR']
def auc(x,y):
    n1,n2=len(x),len(y); z=np.concatenate([x,y]); o=np.argsort(z,kind='mergesort'); zs=z[o]
    r=np.empty(len(z)); i=0
    while i<len(z):
        j=i
        while j+1<len(z) and zs[j+1]==zs[i]: j+=1
        r[i:j+1]=(i+j)/2.0+1; i=j+1
    rr=np.empty(len(z)); rr[o]=r
    return (rr[:n1].sum()-n1*(n1+1)/2.0)/(n1*n2)

files=sorted(glob.glob(os.path.join(RAW,'*.fcs')))[:14]
for PCT in [40,50,60,70]:
    acc={a:[] for a in ANT}; means={a:[] for a in ANT}
    for f in files:
        raw=open(f,'rb').read(); h=a12.parse_fcs3_header(raw); tk,st=a12.parse_fcs3_text_segment(raw,h)
        tot=int(tk['TOT']); npar=int(tk['PAR'])
        arr=np.frombuffer(raw,dtype='>f4',count=tot*npar,offset=h['data_start']).reshape(tot,npar)
        col=lambda m: arr[:,a12.find_parameter_index_by_stain_substring(st,m)].astype(np.float32)
        dna=col('193Ir-DNA'); via=col('195Pt-Viability'); bd=col('140Ce-Beads'); cd45=col('CD45')
        cd3=col('CD3'); cd19=col('CD19'); cd20=col('CD20')
        qc=(dna>np.percentile(dna,5))&(dna<np.percentile(dna,99.5))&(via<np.percentile(via,95))&(bd<np.percentile(bd,99))&(cd45>np.percentile(cd45,2))
        c3=np.percentile(cd3[qc],PCT); c19=np.percentile(cd19[qc],PCT); c20=np.percentile(cd20[qc],PCT)
        B=qc&((cd19>c19)|(cd20>c20))&(cd3<=c3); T=qc&(cd3>c3)&(cd19<=c19)&(cd20<=c20)
        if B.sum()<500 or T.sum()<500: continue
        bi=np.flatnonzero(B); ti=np.flatnonzero(T)
        if len(bi)>12000: bi=rng.choice(bi,12000,replace=False)
        if len(ti)>12000: ti=rng.choice(ti,12000,replace=False)
        for a in ANT:
            v=col(a); acc[a].append(auc(v[bi].astype(float),v[ti].astype(float))); means[a].append(float(v[B].mean()))
    ma={a:np.median(acc[a]) for a in ANT}; mm={a:np.median(means[a]) for a in ANT}
    ra={a:i+1 for i,a in enumerate(sorted(ANT,key=lambda x:-ma[x]))}
    rm={a:i+1 for i,a in enumerate(sorted(ANT,key=lambda x:-mm[x]))}
    sep=lambda r: sum(1 for v in VAL for b in BROAD if r[v]<r[b])/(len(VAL)*len(BROAD))
    print('gate pct=%d  n_used=%d | separation AUC-rank=%.3f  mean-rank=%.3f | VAL ranks(AUC)=%s' % (
        PCT, len(acc['CD19']), sep(ra), sep(rm), [ra[v] for v in VAL]))
