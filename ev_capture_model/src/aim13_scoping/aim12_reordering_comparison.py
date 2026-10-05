"""
Aim 12's reordering metric, recomputed three ways (ungated / B-gated / light-chain clone),
plus the exact continuum version enabled by escape_fraction == F(n*).
"""
import os, sys, glob, json, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE,'..')))
from aim10_car_t_kill_probability_model import fit_all_real_thresholds
DATA=os.path.abspath(os.path.join(HERE,'..','..','..','ev_capture_data','cd19_populations'))

fits=fit_all_real_thresholds()
KS=[fits['CLN-978_cytolysis']['k_intra'], fits['blinatumomab_cytolysis']['k_intra'], fits['blinatumomab_ifn_gamma']['k_intra']]
TAUS=np.arange(0.1,0.95,0.1)
nstar=lambda k,tau,z=1.0: np.log(1.0+tau/(z*(1.0-tau)))/(1.0*k)
GRID=[nstar(k,t) for k in KS for t in TAUS]           # the original 27-point grid, as n* values
pops={}
for p in sorted(glob.glob(os.path.join(DATA,'*.npz'))):
    d=np.load(p); pops[os.path.basename(p)[:-4]]={k:d[k] for k in d.files}
pids=sorted(pops)
print('patients: %d ; grid spans n* = %.1f to %.0f molecules/cell' % (len(pids), min(GRID), max(GRID)))

def reordering(kind, min_spread=0.02):
    ps=[p for p in pids if kind in pops[p] and len(pops[p][kind])>=500]
    arrs={p:pops[p][kind].astype(float) for p in ps}
    means=np.array([arrs[p].mean() for p in ps])
    naive=np.argsort(-means)                      # rank by mean CD19, descending
    meaningful=0
    max_spread=0.0
    for ns in GRID:
        esc=np.array([np.mean(arrs[p]<ns) for p in ps])   # escape fraction == CDF at n*
        model=np.argsort(esc)                             # lower escape = better, descending quality
        flipped=False
        for i in range(len(ps)):
            for j in range(i+1,len(ps)):
                a,b=ps[i],ps[j]
                ia,ib=list(naive).index(i), list(naive).index(j)
                ja,jb=list(model).index(i), list(model).index(j)
                if (ia<ib) != (ja<jb) and abs(esc[i]-esc[j])>=min_spread:
                    flipped=True; max_spread=max(max_spread,abs(esc[i]-esc[j]))
        if flipped: meaningful+=1
    return len(ps), meaningful, len(GRID), max_spread

print('\n%-10s %8s %28s %14s' % ('POPULATION','n','grid points w/ meaningful reorder','max spread'))
res={}
for kind in ['ungated','Bgated','clone']:
    n,m,g,sp=reordering(kind)
    res[kind]=dict(n=n,reordered=m,grid=g,frac=m/g,max_spread=sp)
    print('%-10s %8d %20d/%d (%5.1f%%) %13.1f pp' % (kind,n,m,g,100*m/g,100*sp))

print('\n--- the continuum version (no arbitrary grid) ---')
for kind in ['ungated','Bgated','clone']:
    ps=[p for p in pids if kind in pops[p] and len(pops[p][kind])>=500]
    arrs={p:pops[p][kind].astype(float) for p in ps}
    sweep=np.logspace(np.log10(min(GRID)), np.log10(max(GRID)), 400)
    cdf={p:np.array([np.mean(arrs[p]<s) for s in sweep]) for p in ps}
    crossing=0; total=0
    for i in range(len(ps)):
        for j in range(i+1,len(ps)):
            total+=1
            d=cdf[ps[i]]-cdf[ps[j]]
            d=d[np.abs(d)>=0.02]
            if len(d)>0 and (d.min()<0<d.max()): crossing+=1
    print('  %-8s %4d/%4d patient pairs (%.1f%%) have CDFs that CROSS within the real n* range'
          ' -> mean-based ranking is provably wrong for some achievable CAR sensitivity'
          % (kind, crossing, total, 100*crossing/total))
    res[kind]['crossing_pairs']=crossing; res[kind]['total_pairs']=total
json.dump(res, open(os.path.join(HERE,'aim12_reordering_comparison.json'),'w'), indent=1)
