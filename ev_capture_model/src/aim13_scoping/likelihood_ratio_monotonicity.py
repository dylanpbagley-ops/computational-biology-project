"""
Check three things on REAL data:
 (1) any MONOTONE kill function's (normal-killed, tumor-killed) point lies on/under
     the ROC convex hull of the two per-cell distributions;
 (2) a NON-monotone (band-pass) kill function CAN beat every threshold, when the
     empirical likelihood ratio is non-monotone;
 (3) whether real antigens in this panel actually have non-monotone likelihood ratios.
"""
import sys, os, numpy as np
SRC=os.path.expanduser('~/mnt/outputs/ev_capture_model/src'); sys.path.insert(0,SRC)
import aim12_single_cell_cytof_reordering as a12
from aim10_car_t_kill_probability_model import kill_probability
RAW=os.path.expanduser('~/mnt/outputs/ev_capture_data/raw_cytof')
rng=np.random.default_rng(7)

f=os.path.join(RAW,'CytB2_DLB_P5904.fcs')
raw=open(f,'rb').read(); h=a12.parse_fcs3_header(raw); tk,st=a12.parse_fcs3_text_segment(raw,h)
tot=int(tk['TOT']); npar=int(tk['PAR'])
arr=np.frombuffer(raw,dtype='>f4',count=tot*npar,offset=h['data_start']).reshape(tot,npar)
col=lambda m: arr[:,a12.find_parameter_index_by_stain_substring(st,m)].astype(float)
cd3=col('CD3'); cd19=col('CD19'); cd20=col('CD20'); dna=col('193Ir-DNA'); via=col('195Pt-Viability')
qc=(dna>np.percentile(dna,5))&(dna<np.percentile(dna,99.5))&(via<np.percentile(via,95))
c3=np.percentile(cd3[qc],50); c19=np.percentile(cd19[qc],50); c20=np.percentile(cd20[qc],50)
B=qc&((cd19>c19)|(cd20>c20))&(cd3<=c3); T=qc&(cd3>c3)&(cd19<=c19)&(cd20<=c20)

def roc_hull(xt, xn, grid):
    """returns arrays (fpr, tpr) over thresholds, plus upper concave hull."""
    fpr=np.array([np.mean(xn> g) for g in grid]); tpr=np.array([np.mean(xt> g) for g in grid])
    pts=sorted(set(zip(np.round(fpr,10),np.round(tpr,10))))
    pts=[(0.0,0.0)]+pts+[(1.0,1.0)]
    pts=sorted(set(pts))
    # upper convex hull (monotone chain)
    hull=[]
    for p in pts:
        while len(hull)>=2:
            (x1,y1),(x2,y2)=hull[-2],hull[-1]
            if (x2-x1)*(p[1]-y1) - (y2-y1)*(p[0]-x1) >= 0: hull.pop()
            else: break
        hull.append(p)
    return fpr,tpr,np.array(hull)

def hull_y(hull,x):
    hx,hy=hull[:,0],hull[:,1]
    return np.interp(x,hx,hy)

print('=== (1) do monotone kill functions stay under the ROC convex hull? ===')
for marker in ['CD19','CD38','CD22','HLA-DR']:
    v=col(marker); xt=v[B]; xn=v[T]
    grid=np.unique(np.percentile(np.concatenate([xt,xn]), np.linspace(0,100,400)))
    fpr,tpr,hull=roc_hull(xt,xn,grid)
    worst=0.0
    # sweep the project's own kill_probability across a wide monotone family
    for k in np.logspace(-6,-1,12):
        for za in [0.2,1.0,5.0]:
            p_t=kill_probability(xt,1.0,k,za).mean(); p_n=kill_probability(xn,1.0,k,za).mean()
            worst=max(worst, p_t - hull_y(hull,p_n))
    # also a family of Hill (ultrasensitive) functions, incl. very steep
    for n_hill in [1,2,4,8,20]:
        for kd in np.percentile(np.concatenate([xt,xn]),[10,30,50,70,90]):
            kd=max(kd,1e-9)
            g=lambda x: x**n_hill/(x**n_hill+kd**n_hill)
            worst=max(worst, g(xt).mean()-hull_y(hull,g(xn).mean()))
    print('  %-8s AUC=%.4f   max excess above hull over 96 monotone kill fns = %+.2e' % (
        marker, np.mean(xt[:,None][rng.choice(len(xt),3000)][:,0][:,None] > xn[rng.choice(len(xn),3000)][None,:]), worst))

print()
print('=== (2)+(3) is the empirical likelihood ratio monotone for real antigens? ===')
print('    (non-monotone LR  =>  a band-pass recognizer can beat EVERY density threshold)')
for marker in ['CD19','CD20','CD22','CD79B','CD38','HLA-DR','CD44','CD27','CD48','CD25','IgM','CD21']:
    v=col(marker); xt=v[B]; xn=v[T]
    edges=np.unique(np.percentile(np.concatenate([xt,xn]), np.linspace(0,100,21)))
    if len(edges)<6: print('  %-8s (too few distinct values)'%marker); continue
    ht,_=np.histogram(xt,bins=edges); hn,_=np.histogram(xn,bins=edges)
    ht=ht/ht.sum(); hn=hn/hn.sum()
    lr=(ht+1e-9)/(hn+1e-9)
    incr=np.all(np.diff(lr)>=-1e-12)
    # is a threshold rule optimal? compare best threshold youden vs best LR-based rule
    grid=edges[1:-1]
    best_thr=max(np.mean(xt>g)-np.mean(xn>g) for g in grid)
    # LR rule: kill bins with lr > c, sweep c
    best_lr=0.0
    for c in np.unique(lr):
        sel=lr>=c
        best_lr=max(best_lr, ht[sel].sum()-hn[sel].sum())
    print('  %-8s LR monotone increasing? %-5s | best Youden threshold=%.3f  best LR-rule=%.3f  gain=%+.3f' % (
        marker, incr, best_thr, best_lr, best_lr-best_thr))
