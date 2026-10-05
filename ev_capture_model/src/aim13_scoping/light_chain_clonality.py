"""
Light-chain restriction: separate CLONAL (malignant) B cells from residual POLYCLONAL
(normal) B cells inside each DLBCL biopsy, using the panel's own IgKappa / IgLambda channels.
This is the standard clinical flow-cytometry criterion for B-cell clonality.

Two uses:
  (1) purify the tumor gate;
  (2) give a SAME-TUBE internal normal control (the minority light chain) -- same patient,
      same run, same staining -- against which the cross-sample RLN comparator is validated.

METHOD NOTE (a real bug caught and fixed, recorded so it is not repeated): the first
version set the kappa and lambda positivity cutoffs at the same PERCENTILE of each
channel's own B-cell distribution. That forces |kappa+| == |lambda+| exactly, by
construction -- every sample returned K:L = 1.00 and the method was vacuous. The cutoff
must come from a population that is genuinely NEGATIVE for the marker. T cells in the
same file are surface-immunoglobulin negative, so their own IgKap/IgLam distribution is
the correct internal negative control. That is what is used below: it is per-file (so
invariant to the ~10x cross-batch signal-scale difference) and biologically principled.
"""
import os, sys, glob, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.abspath(os.path.join(HERE,'..','..','..','ev_capture_data','gated_cache'))
files = {os.path.basename(p).replace('.npz',''): p for p in glob.glob(os.path.join(CACHE,'*.npz'))}
ANT = list(np.load(next(iter(files.values())), allow_pickle=True)['antigens'])
IK, IL = ANT.index('IgKap'), ANT.index('IgLam')
NEG_PCT = 99.0   # cutoff = 99th percentile of the sIg-negative T-cell population

def clonality(key):
    d = np.load(files[key], allow_pickle=True)
    B, T = d['B'].astype(float), d['T'].astype(float)
    if len(B) < 500 or len(T) < 500: return None
    kc, lc = np.percentile(T[:,IK], NEG_PCT), np.percentile(T[:,IL], NEG_PCT)
    k, l = B[:,IK], B[:,IL]
    kpos, lpos = (k > kc) & (l <= lc), (l > lc) & (k <= kc)
    dn = (k <= kc) & (l <= lc)
    nk, nl = int(kpos.sum()), int(lpos.sum())
    return dict(n_B=len(B), kappa=nk, lam=nl, dn=int(dn.sum()),
                sig_pos_frac=float((nk+nl)/len(B)),
                ratio=float(nk/max(nl,1)), kpos=kpos, lpos=lpos, B=B)

print('=== calibration: RLN controls must read POLYCLONAL ===')
rl=[]
for key in sorted(k for k in files if '_RLN_' in k and 'pooled' not in k):
    c=clonality(key)
    if not c: continue
    rl.append(c['ratio'])
    if len(rl)<=6: print('  %-22s kappa=%6d lambda=%6d  K:L=%5.2f  sIg+ =%4.1f%%' % (
        key, c['kappa'], c['lam'], c['ratio'], 100*c['sig_pos_frac']))
rl=np.array(rl)
LO,HI = float(np.percentile(rl,2.5)), float(np.percentile(rl,97.5))
print('  ALL %d RLN donors: K:L median %.2f, range %.2f-%.2f' % (len(rl), np.median(rl), rl.min(), rl.max()))
print('  -> polyclonal reference interval (2.5-97.5 pct of healthy donors): %.2f - %.2f' % (LO,HI))

print('\n=== DLBCL tumors ===')
calls={}
for key in sorted(k for k in files if '_DLB_' in k):
    c=clonality(key)
    if not c: continue
    # Clinical convention for B-cell clonality is kappa:lambda > 3:1 or < 1:2.
    # Take whichever is MORE conservative: the clinical rule, or the empirical healthy-donor
    # interval measured above. Also require a real surface-Ig+ population (>=5% of B cells),
    # so surface-Ig-negative tumors are not called restricted on a handful of events.
    hi_cut, lo_cut = max(3.0, HI), min(0.5, LO)
    if c['sig_pos_frac'] < 0.05:  call, nclone, ninternal = 'sIg-negative',      0, 0
    elif c['ratio'] > hi_cut:     call, nclone, ninternal = 'kappa-restricted',  c['kappa'], c['lam']
    elif c['ratio'] < lo_cut:     call, nclone, ninternal = 'lambda-restricted', c['lam'],   c['kappa']
    else:                         call, nclone, ninternal = 'not restricted',    0, 0
    calls[key]=dict(call=call, ratio=c['ratio'], n_clone=nclone, n_internal=ninternal,
                    sig_pos_frac=c['sig_pos_frac'], n_B=c['n_B'])
    print('  %-20s K:L=%8.2f  %-18s clone=%6d internal-normal=%5d  sIg+=%4.1f%%' % (
        key.split('_DLB_')[-1]+' ('+key.split('_')[0]+')', c['ratio'], call, nclone, ninternal,
        100*c['sig_pos_frac']))
from collections import Counter
print('\n  calls:', Counter(v['call'] for v in calls.values()).most_common())
usable=[k for k,v in calls.items() if v['call'] in ('kappa-restricted','lambda-restricted') and v['n_clone']>=500 and v['n_internal']>=150]
print('  patients with a clone (>=500 cells) AND an internal normal (>=150 cells): %d/%d' % (len(usable),len(calls)))
for k in usable: print('     ', k)
json.dump({'calls':calls,'usable':usable,'polyclonal_interval':[LO,HI]},
          open(os.path.join(HERE,'light_chain_calls.json'),'w'), indent=1)
