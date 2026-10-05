"""
(A) NEGATIVE CONTROL: normal-B vs normal-B across independent RLN samples in the same
    batch, and pooled-batch-control vs individual RLN. Should sit at AUC~0.5 for every
    antigen if the gating + rank statistic are not manufacturing differences.
(B) CROSS-BATCH CONTROL: same, but across different batches -- quantifies how much of an
    AUC could be batch artifact.
(C) JOINT WINDOW: a target must spare BOTH normal compartments, so the achievable
    therapeutic index is bounded by the WORSE of the two comparators.
"""
import os, sys, glob, json, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
CACHE=os.path.abspath(os.path.join(HERE,'..','..','..','ev_capture_data','gated_cache'))
rng=np.random.default_rng(5); SUB=12000
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
files={os.path.basename(p).replace('.npz',''):p for p in glob.glob(os.path.join(CACHE,'*.npz'))}
ANT=list(np.load(next(iter(files.values())),allow_pickle=True)['antigens'])
rlns=sorted(k for k in files if '_RLN_' in k and 'pooled' not in k)
pool=sorted(k for k in files if 'pooled' in k)
bat=lambda k:k.split('_')[0]
getB=lambda k: np.load(files[k],allow_pickle=True)['B']

print('=== (A) within-batch normal-B vs normal-B (independent donors) ===')
within=[]
for b in sorted({bat(k) for k in rlns}):
    ks=[k for k in rlns if bat(k)==b]
    for i in range(len(ks)):
        for j in range(i+1,len(ks)):
            x,y=sub(getB(ks[i])),sub(getB(ks[j]))
            within.append([auc(x[:,m].astype(float),y[:,m].astype(float)) for m in range(len(ANT))])
within=np.array(within)
dev=np.abs(within-0.5)
print('  %d donor pairs. |AUC-0.5| across all antigens: median %.3f, 90th pct %.3f, max %.3f' % (
    len(within), np.median(dev), np.percentile(dev,90), dev.max()))
print('  antigens with the largest normal-vs-normal spread (biological donor variation):')
mm=np.median(dev,axis=0); order=np.argsort(-mm)[:6]
for i in order: print('     %-10s median |AUC-0.5| = %.3f' % (ANT[i], mm[i]))

print('\n=== (B) cross-batch normal-B vs normal-B (batch artifact ceiling) ===')
cross=[]
bs=sorted({bat(k) for k in rlns})
for i in range(len(bs)):
    for j in range(i+1,len(bs)):
        a=[k for k in rlns if bat(k)==bs[i]][:2]; c=[k for k in rlns if bat(k)==bs[j]][:2]
        for k1 in a:
            for k2 in c:
                x,y=sub(getB(k1)),sub(getB(k2))
                cross.append([auc(x[:,m].astype(float),y[:,m].astype(float)) for m in range(len(ANT))])
cross=np.array(cross); dev2=np.abs(cross-0.5)
print('  %d cross-batch pairs. |AUC-0.5|: median %.3f, 90th pct %.3f, max %.3f' % (
    len(cross), np.median(dev2), np.percentile(dev2,90), dev2.max()))
print('  -> cross-batch inflation over within-batch: median %+.3f' % (np.median(dev2)-np.median(dev)))

print('\n=== (C) joint window: achievable index is bounded by the WORSE comparator ===')
r=json.load(open(os.path.join(HERE,'two_comparator_results.json')))['summary']
VAL=['CD19','CD20','CD22','CD79B']
rows=[(a, r[a]['medT'], r[a]['medNB'], min(r[a]['medT'], r[a]['medNB'])) for a in ANT]
rows.sort(key=lambda x:-x[3])
print('%-10s %9s %11s %11s' % ('ANTIGEN','vs T','vs normalB','JOINT(min)'))
for a,t,n,j in rows:
    tag='  <-- approved drug exists' if a in VAL else ''
    print('%-10s %9.3f %11.3f %11.3f%s' % (a,t,n,j,tag))
print('\nBest joint window of any antigen in the panel: %.3f (%s)' % (rows[0][3], rows[0][0]))
print('Best joint window among the 4 approved-drug targets: %.3f' % max(min(r[a]['medT'],r[a]['medNB']) for a in VAL))
json.dump({'within_batch_dev_median':float(np.median(dev)),'cross_batch_dev_median':float(np.median(dev2)),
           'joint':{a:min(r[a]['medT'],r[a]['medNB']) for a in ANT}}, open(os.path.join(HERE,'joint_window_results.json'),'w'))
