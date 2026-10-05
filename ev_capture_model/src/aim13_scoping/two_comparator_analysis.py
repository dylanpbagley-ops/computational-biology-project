"""
Therapeutic-index ceiling per antigen, under TWO comparators:
  (1) bystander non-B leukocytes (CD3+ T cells) from the SAME file   -> batch-cancelling
  (2) normal B cells from batch-matched reactive lymph node controls -> the real
      on-target/off-tumor compartment for B-lineage targets
AUC = P(random tumor B cell carries more antigen than a random comparator cell),
computed by exact rank statistic. AUC is invariant to any monotone per-channel
rescaling, so it is immune to the ~10x cross-batch signal-scale difference.
"""
import os, sys, glob, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.abspath(os.path.join(HERE,'..','..','..','ev_capture_data','gated_cache'))
rng = np.random.default_rng(3)
SUB = 12000

def auc(x, y):
    n1, n2 = len(x), len(y)
    z = np.concatenate([x, y]); o = np.argsort(z, kind='mergesort'); zs = z[o]
    r = np.empty(len(z)); i = 0
    while i < len(z):
        j = i
        while j+1 < len(z) and zs[j+1] == zs[i]: j += 1
        r[i:j+1] = (i+j)/2.0 + 1; i = j+1
    rr = np.empty(len(z)); rr[o] = r
    return (rr[:n1].sum() - n1*(n1+1)/2.0) / (n1*n2)

def sub(a):
    if len(a) > SUB: return a[rng.choice(len(a), SUB, replace=False)]
    return a

files = {os.path.basename(p).replace('.npz',''): p for p in glob.glob(os.path.join(CACHE,'*.npz'))}
ANT = list(np.load(next(iter(files.values())), allow_pickle=True)['antigens'])
tumors = sorted(k for k in files if '_DLB_' in k)
rlns   = sorted(k for k in files if '_RLN_' in k and 'pooled' not in k)
pooled = sorted(k for k in files if 'pooled' in k)
batch = lambda k: k.split('_')[0]

# batch-matched normal-B reference: pool all RLN B cells from the same batch
normalB = {}
for b in sorted({batch(k) for k in rlns}):
    mats = [np.load(files[k], allow_pickle=True)['B'] for k in rlns if batch(k)==b]
    m = np.concatenate(mats, axis=0)
    normalB[b] = m[rng.choice(len(m), min(len(m), 40000), replace=False)]
    print('normal-B reference  batch %s: %d RLN samples, %d cells' % (b, sum(1 for k in rlns if batch(k)==b), len(m)))

res = {a: {'vsT': [], 'vsNormalB': []} for a in ANT}
per_patient = {}
for t in tumors:
    d = np.load(files[t], allow_pickle=True)
    B, T = d['B'], d['T']
    nb = normalB[batch(t)]
    if len(B) < 500 or len(T) < 500: print('skip', t); continue
    Bs, Ts, Ns = sub(B), sub(T), sub(nb)
    rec = {}
    for i, a in enumerate(ANT):
        v1 = auc(Bs[:,i].astype(float), Ts[:,i].astype(float))
        v2 = auc(Bs[:,i].astype(float), Ns[:,i].astype(float))
        res[a]['vsT'].append(v1); res[a]['vsNormalB'].append(v2)
        rec[a] = {'vsT': float(v1), 'vsNormalB': float(v2), 'mean_tumor': float(B[:,i].mean())}
    per_patient[t] = rec
    print('  %-22s CD19 vsT=%.3f vsNormalB=%.3f | CD20 vsNormalB=%.3f | CD10 vsNormalB=%.3f' % (
        t, rec['CD19']['vsT'], rec['CD19']['vsNormalB'], rec['CD20']['vsNormalB'], rec['CD10']['vsNormalB']), flush=True)

summary = {a: {'medT': float(np.median(res[a]['vsT'])), 'medNB': float(np.median(res[a]['vsNormalB'])),
               'iqrNB': [float(np.percentile(res[a]['vsNormalB'],25)), float(np.percentile(res[a]['vsNormalB'],75))]}
           for a in ANT}
json.dump({'summary': summary, 'per_patient': per_patient}, open(os.path.join(HERE,'two_comparator_results.json'),'w'))
print('\n%-10s %10s %12s %20s' % ('ANTIGEN','AUC vs T','AUC vs normalB','IQR vs normalB'))
for a in sorted(ANT, key=lambda x: -summary[x]['medNB']):
    s = summary[a]
    print('%-10s %10.3f %12.3f      [%.3f, %.3f]' % (a, s['medT'], s['medNB'], s['iqrNB'][0], s['iqrNB'][1]))
