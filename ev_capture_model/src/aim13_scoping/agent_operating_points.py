"""
Turn the ceiling into a statement about REAL agents.

n* = ln(1 + tau/(z(1-tau))) / (sigma*k_intra) is the antigen-density threshold above which
a cell is reliably killed. At tau=0.5 it equals exactly the measured activation thresholds
from Aim 10's real calibrations: 81 (CLN-978), 843 (blinatumomab cytolysis), 8,690
(blinatumomab IFN-gamma) molecules/cell. For each antigen, placing n* on the tumor and
normal CD19 distributions gives the operating point an agent of that sensitivity reaches.
"""
import os, sys, glob, json, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE,'..')))
CACHE=os.path.abspath(os.path.join(HERE,'..','..','..','ev_capture_data','gated_cache'))
rng=np.random.default_rng(23)
files={os.path.basename(p).replace('.npz',''):p for p in glob.glob(os.path.join(CACHE,'*.npz'))}
ANT=list(np.load(next(iter(files.values())),allow_pickle=True)['antigens'])
bat=lambda k:k.split('_')[0]
tumors=sorted(k for k in files if '_DLB_' in k)
rlns=sorted(k for k in files if '_RLN_' in k and 'pooled' not in k)
normalB={b: np.concatenate([np.load(files[k],allow_pickle=True)['B'] for k in rlns if bat(k)==b],axis=0).astype(float)
         for b in sorted({bat(k) for k in rlns})}
AGENTS=[('CLN-978 (T-cell engager)',81.0),('blinatumomab, cytolysis',843.0),('blinatumomab, IFN-gamma',8690.0)]

print('Operating points at the three REAL measured activation thresholds.')
print('"tumor killed" = fraction of tumor cells above n*; "normal B killed" = same for normal B cells.')
print('NOTE: CyTOF signal is not molecule-calibrated, so n* is applied on a same-order-of-magnitude')
print('but not unit-confirmed scale. Read as a shape/placement statement, not an absolute prediction.\n')
for aname, ns in AGENTS:
    print('--- %s : n* = %.0f ---' % (aname, ns))
    for anti in ['CD19','CD20','CD22','CD79B']:
        i=ANT.index(anti)
        tk=[]; nk=[]
        for t in tumors:
            B=np.load(files[t],allow_pickle=True)['B'].astype(float)
            if len(B)<500: continue
            tk.append(np.mean(B[:,i]>ns)); nk.append(np.mean(normalB[bat(t)][:,i]>ns))
        print('    %-6s tumor killed %5.1f%%   normal B killed %5.1f%%   window %+5.1f pp' % (
            anti, 100*np.median(tk), 100*np.median(nk), 100*(np.median(tk)-np.median(nk))))
    print()
