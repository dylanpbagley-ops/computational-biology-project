"""Compute the Figure 5 worked example on B-GATED tumor cells (the version that survives sorting)."""
import os, sys, glob, json, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.abspath(os.path.join(HERE,'..')))
from aim10_car_t_kill_probability_model import fit_all_real_thresholds, escape_fraction
D=os.path.abspath(os.path.join(HERE,'..','..','..','ev_capture_data','cd19_populations'))
pops={os.path.basename(p)[:-4].split('_DLB_')[-1]: np.load(p)['Bgated'].astype(float)
      for p in sorted(glob.glob(os.path.join(D,'*.npz')))}
k=fit_all_real_thresholds()['CLN-978_cytolysis']['k_intra']
naive=sorted(pops,key=lambda p:-pops[p].mean())
esc={p: escape_fraction(pops[p],1.0,k,1.0,kill_threshold=0.5)['escape_fraction'] for p in pops}
model=sorted(pops,key=lambda p:esc[p])
print('B-GATED, CLN-978_cytolysis, kill_threshold=0.5 (n* = 81 molecules/cell)\n')
print('naive top-5 by mean :', [f'{p} (mean={pops[p].mean():.1f})' for p in naive[:5]])
print('model top-5 by escape:', [f'{p} (escape={100*esc[p]:.1f}%)' for p in model[:5]])
disp=set(naive[:5])-set(model[:5]); prom=set(model[:5])-set(naive[:5])
print('\nDisplaced (drops out of top 5):', disp)
print('Promoted  (new entry to top 5):', prom)
for p in sorted(disp|prom):
    print('   %-7s mean=%8.2f  escape=%5.1f%%  naive rank=%2d  model rank=%2d'
          % (p, pops[p].mean(), 100*esc[p], naive.index(p)+1, model.index(p)+1))
if disp and prom:
    d,q=sorted(disp)[0],sorted(prom)[0]
    print('\nDirection check: promoted %s has mean %.1f, displaced %s has mean %.1f -> %s'
          % (q,pops[q].mean(),d,pops[d].mean(),
             'promoted has the LOWER mean (the point holds)' if pops[q].mean()<pops[d].mean()
             else 'promoted has the HIGHER mean (point does NOT hold for this pair)'))
json.dump({'naive_top5':[[p,float(pops[p].mean())] for p in naive[:5]],
           'model_top5':[[p,float(esc[p])] for p in model[:5]],
           'displaced':sorted(disp),'promoted':sorted(prom),
           'all_means':{p:float(pops[p].mean()) for p in pops},
           'all_escape':{p:float(esc[p]) for p in pops}},
          open(os.path.join(HERE,'gated_reordering_example.json'),'w'), indent=1)
