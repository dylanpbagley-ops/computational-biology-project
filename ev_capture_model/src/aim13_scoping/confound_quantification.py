import json, os, numpy as np
d = json.load(open(os.path.expanduser('~/scratch/cohort_scan.json')))
ids = sorted(d)
def rank(x):
    x=np.asarray(x,float); o=np.argsort(np.argsort(x)); return o.astype(float)
def spearman(a,b):
    ra,rb=rank(a),rank(b); ra=ra-ra.mean(); rb=rb-rb.mean()
    return float((ra@rb)/np.sqrt((ra@ra)*(rb@rb)))
allm=[d[i]['cd19_mean_all'] for i in ids]
gm  =[d[i].get('cd19_mean_B',np.nan) for i in ids]
bf  =[d[i]['bfrac'] for i in ids]
tf  =[d[i]['tfrac'] for i in ids]
dna =[d[i]['dna_med'] for i in ids]
print('n samples', len(ids))
print('Spearman(ungated mean CD19, B-gated mean CD19) = %.3f' % spearman(allm,gm))
print('Spearman(ungated mean CD19, B-cell fraction)   = %.3f' % spearman(allm,bf))
print('Spearman(ungated mean CD19, T-cell fraction)   = %.3f' % spearman(allm,tf))
print('Spearman(ungated mean CD19, DNA median[batch]) = %.3f' % spearman(allm,dna))
print('Spearman(B-gated mean CD19, DNA median[batch]) = %.3f' % spearman(gm,dna))
print('B-cell fraction range: %.3f - %.3f' % (min(bf),max(bf)))
print('DNA median range: %.1f - %.1f  (ratio %.1fx)' % (min(dna),max(dna),max(dna)/min(dna)))
# batch groups
import collections
g=collections.defaultdict(list)
for i in ids: g[d[i]['batch']].append(d[i]['dna_med'])
for k in sorted(g): print('  batch %s: n=%d  DNA median = %.1f' % (k,len(g[k]),np.median(g[k])))
# top-5 comparison
o1=[ids[i] for i in np.argsort(allm)[::-1][:5]]
o2=[ids[i] for i in np.argsort(gm)[::-1][:5]]
print('top5 ungated :', [x.split("_")[-1] for x in o1])
print('top5 B-gated :', [x.split("_")[-1] for x in o2])
# how many patients move >5 rank positions
r1,r2=rank(allm),rank(gm)
mv=np.abs(r1-r2)
print('patients moving >5 rank positions: %d/%d ; max move %d' % ((mv>5).sum(),len(ids),mv.max()))
# variance explained by batch on ungated log CD19
lg=np.log(np.array(allm))
bats=np.array([d[i]['batch'] for i in ids])
tot=lg.var()
within=np.mean([ (lg[bats==b]-lg[bats==b].mean())**2 @ np.ones((bats==b).sum())/(bats==b).sum() for b in set(bats)])
ss_w=sum(((lg[bats==b]-lg[bats==b].mean())**2).sum() for b in set(bats))
print('R^2 of batch on log(ungated CD19 mean) = %.3f' % (1-ss_w/((lg-lg.mean())**2).sum()))
lgB=np.log(np.array(gm))
ss_wB=sum(((lgB[bats==b]-lgB[bats==b].mean())**2).sum() for b in set(bats))
print('R^2 of batch on log(B-gated CD19 mean) = %.3f' % (1-ss_wB/((lgB-lgB.mean())**2).sum()))
