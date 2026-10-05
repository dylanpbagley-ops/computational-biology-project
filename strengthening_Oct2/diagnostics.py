"""Fixed-rule, post-discovery diagnostics; no fitting or outcome relabeling."""
import itertools,json,math,pathlib
import numpy as np
from prepare import literal
from frontiers import tau_from_share
ROOT=pathlib.Path(__file__).resolve().parent;SRC=ROOT.parent/'ev_capture_model/src/aim13_scoping'
def auc(safe,toxic):
    a=np.asarray(safe)[:,None];b=np.asarray(toxic)[None,:]
    if not a.size or not b.size:raise ValueError('Both outcomes required')
    return float(np.mean((a>b)+.5*(a==b)))
def exact_p(safe,toxic):
    a=list(safe);b=list(toxic);x=a+b;m=len(a)
    # Twice the midrank retains exact integer arithmetic even with ties.
    ranks=[2*sum(v<t for v in x)+sum(v==t for v in x)+1 for t in x]
    observed=sum(ranks[:m]);dp=[{} for _ in range(m+1)];dp[0][0]=1
    for i,rank in enumerate(ranks):
        for k in range(min(i+1,m),0,-1):
            for total,count in list(dp[k-1].items()):dp[k][total+rank]=dp[k].get(total+rank,0)+count
    assert sum(dp[m].values())==math.comb(len(x),m)
    return sum(c for r,c in dp[m].items() if r>=observed)/math.comb(len(x),m)
def paired_bootstrap(a,b,y,n=20000,seed=20261002):
    a=np.asarray(a);b=np.asarray(b);y=np.asarray(y,bool);s=np.flatnonzero(y);t=np.flatnonzero(~y);rng=np.random.default_rng(seed)
    si=rng.choice(s,size=(n,len(s)));ti=rng.choice(t,size=(n,len(t)))
    def samples(x):
        aa=x[si][:,:,None];bb=x[ti][:,None,:]
        return ((aa>bb)+.5*(aa==bb)).mean(axis=(1,2))
    delta=samples(a)-samples(b)
    return {'difference':auc(a[y],a[~y])-auc(b[y],b[~y]),'ci95':np.quantile(delta,[.025,.975]).tolist(),'draws':n,'seed':seed}
def main():
    targets=literal(SRC/'target_labels.py','TARGETS');labels={t[0]:t[3] for t in targets if t[3] in ('safe','toxic')};genes=list(labels)
    raw=json.loads((SRC/'wholebody_v2.json').read_text())['reference_A_cell_type'];flags=json.loads((SRC/'flag_cleanup.json').read_text())['targets']
    share=np.array([raw[g]['top_share'] for g in genes]);blood=np.array([float(flags[g]['M2_flag']) for g in genes]);score=share+blood;y=np.array([labels[g]=='safe' for g in genes])
    tau=np.array([tau_from_share(x,154) for x in share]);models={'frozen_rule':score,'confinement':share,'blood_flag':blood,'tau':tau,'tau_plus_blood':tau+blood}
    baseline={k:{'auc':auc(v[y],v[~y]),'exact_one_sided_p':exact_p(v[y],v[~y])} for k,v in models.items()}
    comparisons={k:paired_bootstrap(score,v,y) for k,v in models.items() if k!='frozen_rule'}
    groups={'CAR-T':'CD19 TNFRSF17 ERBB2 CA9 FAP GPC3','antibody':'MS4A1 SLAMF7 CD38 CLDN18 CD52 CD47','ADC':'CD22 CD79B TNFRSF8 CD33 TACSTD2 NECTIN4 MET NCAM1 FOLR1 F3 CD44','bispecific':'GPRC5D DLL3','radioligand':'FOLH1','TCR-T':'CEACAM5','mixed':'MSLN EGFR EPCAM'}
    modality={g:k for k,v in groups.items() for g in v.split()};assert set(modality)==set(genes)
    strata={};pairs_total=0;wins_total=0
    for name in groups:
        ids=np.array([i for i,g in enumerate(genes) if modality[g]==name]);s=ids[y[ids]];t=ids[~y[ids]];pairs=len(s)*len(t)
        strata[name]={'safe':len(s),'toxic':len(t),'pairs':pairs,'auc':auc(score[s],score[t]) if pairs else None}
        if name!='mixed' and pairs:pairs_total+=pairs;wins_total+=auc(score[s],score[t])*pairs
    diagnostics={}
    disputed=['FAP','FOLR1','EGFR','EPCAM','NCAM1','CD47','F3']
    for name,drop in [('exclude_FAP_FOLR1',disputed[:2]),('exclude_seven_source_uncertainties',disputed)]:
        keep=np.array([g not in drop for g in genes]);s=score[keep&y];t=score[keep&~y]
        diagnostics[name]={'safe':len(s),'toxic':len(t),'auc':auc(s,t),'exact_one_sided_p':exact_p(s,t),'excluded':drop}
    rows=[]
    for bits in itertools.product([False,True],repeat=len(disputed)):
        yy=y.copy()
        for g,v in zip(disputed,bits):yy[genes.index(g)]=v
        rows.append({'assigned_safe':[g for g,v in zip(disputed,bits) if v],'auc':auc(score[yy],score[~yy])})
    out={'baseline':baseline,'paired_comparisons':comparisons,'modality_strata':strata,'within_modality_pooled_concordance':wins_total/pairs_total,'within_modality_pairs':pairs_total,'deletion_diagnostics':diagnostics,'disputed_assignments':{'n':len(rows),'min':min(rows,key=lambda r:r['auc']),'max':max(rows,key=lambda r:r['auc']),'interpretation':'Fixed-membership uncertainty envelope; no probability on assignments and no label adjudication.'},'status':'Post-discovery diagnostics conditional on the original selected labels; no independently validated outcome prediction.'}
    (ROOT/'diagnostic_results.json').write_text(json.dumps(out,indent=2));(ROOT/'disputed_assignments.json').write_text(json.dumps(rows,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':main()
