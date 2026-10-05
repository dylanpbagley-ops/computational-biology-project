import sys, os, numpy as np, warnings
warnings.filterwarnings('ignore')
SRC=os.path.expanduser('~/mnt/outputs/ev_capture_model/src'); sys.path.insert(0,SRC)
import aim12_single_cell_cytof_reordering as a12
from aim10_car_t_kill_probability_model import kill_probability
RAW=os.path.expanduser('~/mnt/outputs/ev_capture_data/raw_cytof')
f=os.path.join(RAW,'CytB2_DLB_P5904.fcs')
raw=open(f,'rb').read(); h=a12.parse_fcs3_header(raw); tk,st=a12.parse_fcs3_text_segment(raw,h)
tot=int(tk['TOT']); npar=int(tk['PAR'])
arr=np.frombuffer(raw,dtype='>f4',count=tot*npar,offset=h['data_start']).reshape(tot,npar)
col=lambda m: arr[:,a12.find_parameter_index_by_stain_substring(st,m)].astype(float)
cd3=col('CD3'); cd19=col('CD19'); cd20=col('CD20'); dna=col('193Ir-DNA'); via=col('195Pt-Viability')
qc=(dna>np.percentile(dna,5))&(dna<np.percentile(dna,99.5))&(via<np.percentile(via,95))
c3=np.percentile(cd3[qc],50); c19=np.percentile(cd19[qc],50); c20=np.percentile(cd20[qc],50)
B=qc&((cd19>c19)|(cd20>c20))&(cd3<=c3); T=qc&(cd3>c3)&(cd19<=c19)&(cd20<=c20)

def hull_of(xt,xn):
    grid=np.unique(np.percentile(np.concatenate([xt,xn]),np.linspace(0,100,600)))
    pts=set([(0.0,0.0),(1.0,1.0)])
    for g in grid: pts.add((float(np.mean(xn>g)), float(np.mean(xt>g))))
    pts=sorted(pts); hull=[]
    for p in pts:
        while len(hull)>=2:
            (x1,y1),(x2,y2)=hull[-2],hull[-1]
            if (x2-x1)*(p[1]-y1)-(y2-y1)*(p[0]-x1)>=0: hull.pop()
            else: break
        hull.append(p)
    return np.array(hull)
def hy(hull,x): return np.interp(x,hull[:,0],hull[:,1])

tot_eval=0; tot_skip=0; worst_all=-9e9; worst_where=None
for marker in ['CD19','CD20','CD22','CD79B','CD38','HLA-DR','CD44','CD27']:
    v=col(marker); xt=v[B]; xn=v[T]; hull=hull_of(xt,xn)
    fam=[]
    for k in np.logspace(-6,-2,15):
        for za in [0.2,1.0,5.0]:
            fam.append(('project_kill k=%.1e z=%.1f'%(k,za), lambda x,k=k,za=za: kill_probability(x,1.0,k,za)))
    for n_hill in [1,2,4,8,20,50]:
        for q in [5,20,40,60,80,95]:
            kd=max(float(np.percentile(np.concatenate([xt,xn]),q)),1e-9)
            fam.append(('hill n=%d kd=q%d'%(n_hill,q), lambda x,n=n_hill,kd=kd: x**n/(x**n+kd**n)))
    for name,g in fam:
        pt=float(np.mean(g(xt))); pn=float(np.mean(g(xn)))
        if not (np.isfinite(pt) and np.isfinite(pn)): tot_skip+=1; continue
        tot_eval+=1
        exc=pt-hy(hull,pn)
        if exc>worst_all: worst_all=exc; worst_where=(marker,name,pt,pn)
print('monotone kill functions evaluated (finite): %d ; skipped (float overflow): %d' % (tot_eval, tot_skip))
print('LARGEST excess above the ROC convex hull across all of them: %+.3e' % worst_all)
print('  (at %s / %s ; tumor-killed=%.4f normal-killed=%.4f)' % worst_where)
print('  negative or zero  =>  no monotone kill function ever beat the bound')
