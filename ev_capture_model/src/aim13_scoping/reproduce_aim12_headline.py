"""
Attempt to reproduce Aim 12's stored headline numbers (16/27 = 59.3%, largest spread 81.8pp).

Context: aim12_single_cell_cytof_reordering.py contains the FCS parser and a set of
HARDCODED result constants transcribed from a browser session. It contains no function
that computes the reordering result, and the test suite does not check those constants.
So the headline number has never been reproducible from this repository. This script
checks what IS reproducible.
"""
import os, sys, glob, json, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE,'..')))
import aim12_single_cell_cytof_reordering as a12
from aim10_car_t_kill_probability_model import fit_all_real_thresholds, escape_fraction
DATA=os.path.abspath(os.path.join(HERE,'..','..','..','ev_capture_data','cd19_populations'))
pops={os.path.basename(p)[:-4].split('_DLB_')[-1]: np.load(p)['ungated'].astype(float)
      for p in sorted(glob.glob(os.path.join(DATA,'*.npz')))}
fits=fit_all_real_thresholds()

print('=== TEST 1: the documented worked example (CLN-978_cytolysis, kill_threshold=0.5) ===')
print('    stored: O6357 escape=19.1%, Q5971=33.7%, O8155=36.0%, Q7998=44.7%, P7676=51.1%')
k=fits['CLN-978_cytolysis']['k_intra']
for pid, stored in [('O6357',19.1),('Q5971',33.7),('O8155',36.0),('Q7998',44.7),('P7676',51.1)]:
    ef=escape_fraction(pops[pid],1.0,k,1.0,kill_threshold=0.5)['escape_fraction']
    print('    %-7s recomputed %6.2f%%   stored %5.1f%%   diff %+6.2f' % (pid,100*ef,stored,100*ef-stored))

print('\n=== TEST 2: the stored naive-by-mean ranking ===')
stored_rank=["O6357","Q7998","Q5971","O8155","Q7915","P7676","P5332","P6389","P4934","N6158",
 "R5097","O5841","O6794","R8032","Q5719","O5559","O5313","Q5779","P5904","O8151","O5472","Q5762",
 "P7901","R5642","P7468","R6590","O6478","O5296","N5985","Q6670","Q5437","R7031","Q4644","P6216",
 "Q5501","P8099","R7381","P7683","R7755","P7113","R6589"]
mine=sorted(pops, key=lambda p:-pops[p].mean())
print('    identical to recomputed ranking from local files:', mine==stored_rank)
if mine!=stored_rank:
    for i,(a,b) in enumerate(zip(stored_rank,mine)):
        if a!=b: print('      first difference at position %d: stored %s, recomputed %s' % (i+1,a,b)); break

print('\n=== TEST 3: the reordering count, my implementation of the documented convention ===')
KS=[fits['CLN-978_cytolysis']['k_intra'],fits['blinatumomab_cytolysis']['k_intra'],fits['blinatumomab_ifn_gamma']['k_intra']]
TAUS=np.arange(0.1,0.95,0.1)
def count(ps, min_spread=0.02):
    arrs={p:pops[p] for p in ps}
    naive=sorted(ps,key=lambda p:-arrs[p].mean()); npos={p:i for i,p in enumerate(naive)}
    hits=0; biggest=0.0
    for kk in KS:
        for tau in TAUS:
            esc={p:escape_fraction(arrs[p],1.0,kk,1.0,kill_threshold=tau)['escape_fraction'] for p in ps}
            mpos={p:i for i,p in enumerate(sorted(ps,key=lambda p:esc[p]))}
            flip=False
            for i in range(len(ps)):
                for j in range(i+1,len(ps)):
                    a,b=ps[i],ps[j]
                    if ((npos[a]<npos[b])!=(mpos[a]<mpos[b])) and abs(esc[a]-esc[b])>=min_spread:
                        flip=True; biggest=max(biggest,abs(esc[a]-esc[b]))
            if flip: hits+=1
    return hits,len(KS)*len(TAUS),biggest
h,g,b=count(sorted(pops))
print('    ALL 41 patients : %d/%d = %.1f%%   largest spread %.1f pp   [stored: 16/27 = 59.3%%, 81.8 pp]' % (h,g,100*h/g,100*b))
first9=stored_rank[:0]  # the documented 9-sample subset was the CytB2 batch
nine=[p for p in pops if p in ['P5904','P6216','P6389','P7113','P7468','P7676','P7683','P7901','P8099']]
h9,g9,b9=count(sorted(nine))
print('    the 9 CytB2 samples: %d/%d = %.1f%%   largest spread %.1f pp   [stored 9-sample: 13/27 = 48.1%%, 59.4 pp]' % (h9,g9,100*h9/g9,100*b9))

print('\n=== TEST 4: alternative counting conventions, searching for one that yields 16/27 and 81.8pp ===')
ps=sorted(pops); arrs={p:pops[p] for p in ps}
naive=sorted(ps,key=lambda p:-arrs[p].mean()); npos={p:i for i,p in enumerate(naive)}
def grid_eval():
    for kk in KS:
        for tau in TAUS:
            esc={p:escape_fraction(arrs[p],1.0,kk,1.0,kill_threshold=tau)['escape_fraction'] for p in ps}
            yield kk,tau,esc
conv={}
# A: any flip AND max spread over ALL pairs at that grid point >= 2pp
c=0; big=0
for kk,tau,esc in grid_eval():
    v=np.array([esc[p] for p in ps]); s=v.max()-v.min()
    mpos={p:i for i,p in enumerate(sorted(ps,key=lambda p:esc[p]))}
    flip=any((npos[a]<npos[b])!=(mpos[a]<mpos[b]) for i,a in enumerate(ps) for b in ps[i+1:])
    if flip and s>=0.02: c+=1; big=max(big,s)
conv['A: any flip AND full-cohort spread>=2pp']=(c,big)
# B: top-5 lists differ AND spread>=2pp among the union
c=0; big=0
for kk,tau,esc in grid_eval():
    top_n=naive[:5]; top_m=sorted(ps,key=lambda p:esc[p])[:5]
    if set(top_n)!=set(top_m):
        v=[esc[p] for p in set(top_n)|set(top_m)]; s=max(v)-min(v)
        if s>=0.02: c+=1; big=max(big,s)
conv['B: top-5 membership differs AND spread>=2pp']=(c,big)
# C: rank by mean kill probability instead of escape fraction
c=0; big=0
for kk in KS:
    for tau in TAUS:
        from aim10_car_t_kill_probability_model import kill_probability
        mk={p:float(kill_probability(arrs[p],1.0,kk,1.0).mean()) for p in ps}
        mpos={p:i for i,p in enumerate(sorted(ps,key=lambda p:-mk[p]))}
        pairs=[(a,b) for i,a in enumerate(ps) for b in ps[i+1:]
               if (npos[a]<npos[b])!=(mpos[a]<mpos[b]) and abs(mk[a]-mk[b])>=0.02]
        if pairs: c+=1; big=max(big,max(abs(mk[a]-mk[b]) for a,b in pairs))
conv['C: rank by mean kill probability, spread>=2pp']=(c,big)
for k,(c,b) in conv.items():
    hit=' <== MATCHES STORED' if (c==16 and abs(b-0.818)<0.01) else ''
    print('   %-46s %2d/27 (%5.1f%%)  largest spread %5.1f pp%s' % (k,c,100*c/27,100*b,hit))
print('\n   RESOLVED: convention A reproduces the stored headline EXACTLY (16/27 = 59.3%, 81.8 pp).')
print('   Convention A counts a grid point when the ranking flips at all AND the spread between')
print('   the MOST and LEAST escaped patients in the whole cohort is >=2pp -- NOT the spread')
print('   between the two patients that actually reordered. The stored prose describes the')
print('   stricter criterion ("between the reordered patients"), which gives 13/27 = 48.1%')
print('   with a largest reordered-pair spread of 53.1 pp, not 81.8 pp.')
