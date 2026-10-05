import json, os, numpy as np
d = json.load(open(os.path.expanduser('~/scratch/pilot_roc.json')))
pids = sorted(d)
ants = [k for k in d[pids[0]] if isinstance(d[pids[0]][k], dict)]
rows=[]
for a in ants:
    auc=[d[p][a]['auc'] for p in pids if a in d[p]]
    mb =[d[p][a]['mean_B'] for p in pids if a in d[p]]
    rows.append((a, float(np.median(auc)), float(np.percentile(auc,25)), float(np.percentile(auc,75)),
                 float(np.median(mb)), len(auc)))
print('%-10s %7s %15s %12s' % ('ANTIGEN','medAUC','IQR','med mean_B'))
for r in sorted(rows, key=lambda r:-r[1]):
    print('%-10s %7.3f  [%.3f,%.3f]  %12.1f  n=%d' % (r[0],r[1],r[2],r[3],r[4],r[5]))
print()
# ranking comparison
byauc=[r[0] for r in sorted(rows,key=lambda r:-r[1])]
bymean=[r[0] for r in sorted(rows,key=lambda r:-r[4])]
print('TOP 8 by tumor-vs-bystander AUC  :', byauc[:8])
print('TOP 8 by mean tumor expression   :', bymean[:8])
def rk(l): return {a:i for i,a in enumerate(l)}
ra,rm=rk(byauc),rk(bymean)
xs=np.array([ra[a] for a in ants]); ys=np.array([rm[a] for a in ants])
xs=xs-xs.mean(); ys=ys-ys.mean()
print('Spearman(AUC rank, mean rank) = %.3f' % float(xs@ys/np.sqrt((xs@xs)*(ys@ys))))
print()
print('Biggest disagreements (rank_by_mean minus rank_by_AUC):')
diff=sorted(ants, key=lambda a: ra[a]-rm[a])
for a in diff[:5]: print('   %-10s AUC rank %2d  mean rank %2d  (AUC says MUCH better than mean says)'%(a,ra[a]+1,rm[a]+1))
for a in diff[-5:]: print('   %-10s AUC rank %2d  mean rank %2d  (mean says better than AUC says)'%(a,ra[a]+1,rm[a]+1))
