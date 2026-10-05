"""Post-review reproduction, not a new preregistered confirmation."""
import ast, csv, hashlib, json, math
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parent
SRC=ROOT.parent/'ev_capture_model/src/aim13_scoping'
tree=ast.parse((SRC/'target_labels.py').read_text())
targets=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='TARGETS' for t in n.targets))
labels={t[0]:t[3] for t in targets if t[3] in ('safe','toxic')}
shares=json.loads((SRC/'wholebody_v2.json').read_text())['reference_A_cell_type']
flags=json.loads((SRC/'flag_cleanup.json').read_text())['targets']
scores={g:shares[g]['top_share']+int(flags[g]['M2_flag']) for g in labels}

def compute(lab):
    safe=[g for g in lab if lab[g]=='safe']; toxic=[g for g in lab if lab[g]=='toxic']
    wins=sum(scores[a]>scores[b] for a in safe for b in toxic)
    ties=sum(scores[a]==scores[b] for a in safe for b in toxic)
    assert ties==0 and len(set(scores[g] for g in lab))==len(lab)
    order=sorted(lab,key=scores.get)
    observed=sum(i+1 for i,g in enumerate(order) if lab[g]=='safe')
    n=len(lab); m=len(safe); dp=[{} for _ in range(m+1)];dp[0][0]=1
    for rank in range(1,n+1):
        for k in range(min(rank,m),0,-1):
            for total,count in list(dp[k-1].items()):
                dp[k][total+rank]=dp[k].get(total+rank,0)+count
    assert sum(dp[m].values())==math.comb(n,m)
    p=F(sum(c for total,c in dp[m].items() if total>=observed),math.comb(n,m))
    pairs=m*len(toxic)
    return dict(safe=m,toxic=len(toxic),concordant=wins,discordant=pairs-wins,ties=ties,auc=wins/pairs,exact_one_sided_p=float(p),exact_p_fraction=str(p))

scenarios=[('original',[],{}),('exclude_FAP',['FAP'],{}),('exclude_FOLR1',['FOLR1'],{}),('exclude_both',['FAP','FOLR1'],{}),('exclude_FAP_relabel_FOLR1',['FAP'],{'FOLR1':'safe'}),('plus_NECTIN4_toxic',['FAP'],{'FOLR1':'safe','NECTIN4':'toxic'}),('plus_F3_excluded',['FAP','F3'],{'FOLR1':'safe'})]
rows=[]
for name,drop,changes in scenarios:
    lab={g:v for g,v in labels.items() if g not in drop};lab.update(changes)
    rows.append(dict(scenario=name,**compute(lab)))

def lr_auc(p,q):
    return sum(p[i]*q[j]*(F(1) if p[i]/q[i]>p[j]/q[j] else F(1,2) if p[i]/q[i]==p[j]/q[j] else F(0)) for i in range(len(p)) for j in range(len(q)))
q=[F(1,3)]*3;p1=[F(6,10),F(3,10),F(1,10)];p2=[F(55,100),F(4,10),F(5,100)]
assert lr_auc(p1,q)==lr_auc(p2,q)==F(2,3)
out={'status':'Post hoc reproduction of supplied sensitivity scenarios; not clinical label adjudication','source_sha256':{f:hashlib.sha256((SRC/f).read_bytes()).hexdigest() for f in ['target_labels.py','wholebody_v2.json','flag_cleanup.json']},'scenarios':rows,'equal_LR_AUC_example':{'normal':list(map(str,q)),'case1':list(map(str,p1)),'case2':list(map(str,p2)),'LR_AUC_both':'2/3','FPR':'1/3','optimal_TPR_case1':'3/5','optimal_TPR_case2':'11/20'}}
(ROOT/'recheck_results.json').write_text(json.dumps(out,indent=2))
with (ROOT/'label_sensitivity.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(json.dumps(out,indent=2))
