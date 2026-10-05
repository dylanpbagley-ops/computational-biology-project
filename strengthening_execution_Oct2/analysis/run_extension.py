"""Selected narrow B extension, including preserved negative findings."""
import ast,csv,hashlib,json,math,pathlib,sys,platform
import numpy as np
from run_b_pilot import load,verify_inputs,GENES,POLICIES,HERE,ROOT,OUT
from thresholds import select_threshold,fraction_above,auc
from startup import table,read

def freeze():
    p=HERE/'analysis_protocol.md';h=hashlib.sha256(p.read_bytes()).hexdigest();lock=HERE/'frozen_specification.sha256'
    if lock.exists() and lock.read_text().strip()!=h:raise ValueError('Locked protocol changed')
    lock.write_text(h+'\n');return h

def bootstrap(groups,g,rng):
    ts=[v[g] for v in groups['MM'].values()];ns=[v[g] for v in groups['NBM'].values()]
    thresholds=np.r_[0,np.unique(np.concatenate(ts+ns)),np.inf]
    tails=lambda arrays:np.array([(len(x)-np.searchsorted(np.sort(x),thresholds,'left'))/len(x) for x in arrays])
    tt=tails(ts);nt=tails(ns)
    values=[]
    for _ in range(2000):
        ti=rng.integers(0,len(ts),len(ts));ni=rng.integers(0,len(ns),len(ns))
        i=int(np.flatnonzero(nt[ni].mean(axis=0)<=.1+1e-14)[0])
        values.append((float(tt[ti,i].mean()),float(nt[ni,i].max()),float(tt[ti,i].min())))
    a=np.array(values)
    return dict(gene=g,draws=2000,seed=20261002,tumor_mean_ci_low=float(np.quantile(a[:,0],.025)),tumor_mean_ci_high=float(np.quantile(a[:,0],.975)),normal_worst_percentile_low=float(np.quantile(a[:,1],.025)),normal_worst_percentile_high=float(np.quantile(a[:,1],.975)),tumor_worst_percentile_low=float(np.quantile(a[:,2],.025)),tumor_worst_percentile_high=float(np.quantile(a[:,2],.975)),scope='donor bootstrap retunes equal-donor threshold; percentile sampling diagnostics, not guaranteed bounds on donor extremes')

def historical():
    src=ROOT/'ev_capture_model/src/aim13_scoping'
    tree=ast.parse((src/'target_labels.py').read_text(encoding='utf-8'))
    targets=None
    for node in tree.body:
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='TARGETS' for t in node.targets):targets=ast.literal_eval(node.value)
    if targets is None:raise ValueError('Literal TARGETS not found')
    scores=json.loads((src/'wholebody_v2.json').read_text())['reference_A_cell_type']
    flags=json.loads((src/'flag_cleanup.json').read_text())['targets']
    rows=[]
    for t in targets:
        if t[3] not in ['safe','toxic']:continue
        share=scores[t[0]]['top_share'];flag=float(flags[t[0]]['M2_flag'])
        rows.append(dict(gene=t[0],historical_outcome=t[3],top_share=share,blood_flag=flag,tau=(154-1/share)/153,frozen_rule=share+flag))
    result={'status':'Historical selected panel, unchanged labels; no clinical adjudication or independent validation','favorable':sum(r['historical_outcome']=='safe' for r in rows),'unfavorable':sum(r['historical_outcome']=='toxic' for r in rows),'metrics':{k:auc([r[k] for r in rows if r['historical_outcome']=='safe'],[r[k] for r in rows if r['historical_outcome']=='toxic']) for k in ['frozen_rule','top_share','blood_flag','tau']},'score_direction':'larger is historically favorable','tau_equivalence':'known algebraic ranking identity, not independent success','amended_evaluation':'not performed; outcome/negative adequacy and independent review gate remains closed','locked_predictions':'preserved; registration timestamp/independent post-lock clinical outcomes not verified in this execution'}
    if (result['favorable'],result['unfavorable'])!=(19,11):raise ValueError('Historical class count changed')
    table(HERE/'results/historical_inputs.tsv',rows)
    (HERE/'results/historical_reproduction.json').write_text(json.dumps(result,indent=2))
    return result

def main():
    protocol_hash=freeze();verify_inputs();groups=load();rng=np.random.default_rng(20261002)
    intervals=[bootstrap(groups,g,rng) for g in GENES]
    table(HERE/'results/b_bootstrap.tsv',intervals)
    raw=load('raw');secondary=[]
    for g in GENES:
        for budget in [.05,.10,.20,.50]:
            for policy in POLICIES:secondary.append(dict(gene=g,scale='raw_UMI',budget=budget,policy=policy,**select_threshold([v[g] for v in raw['MM'].values()],[v[g] for v in raw['NBM'].values()],budget,policy)))
    table(HERE/'results/b_raw_sensitivity.tsv',secondary)
    held=read(HERE/'results/b_held_out.tsv');donor_rows=[]
    for g in GENES:
        for policy in POLICIES:
            for stage,field,metric in [('NBM','held_normal','normal_fraction'),('MM','held_tumor','tumor_fraction')]:
                for donor in sorted(groups[stage]):
                    rs=[r for r in held if r['gene']==g and r['policy']==policy and r[field]==donor]
                    # These average over partner-held-out folds, not independent replicates.
                    values=[float(r[metric]) for r in rs]
                    donor_rows.append(dict(gene=g,policy=policy,stage=stage,donor=donor,partner_folds=len(rs),mean_held_out_fraction=float(np.mean(values)),min_fraction=min(values),max_fraction=max(values)))
    table(HERE/'results/b_held_out_donors.tsv',donor_rows)
    adt=ROOT/'ev_capture_data/myeloma_adt/01.npz'
    with np.load(adt,allow_pickle=True) as d:adt_audit={'file':str(adt.relative_to(ROOT)), 'keys':list(d.files),'shape':list(d['adt'].shape),'names':list(map(str,d['names'])),'replication_eligibility':'different protein assay and population comparator; target-marker-based tumor gates in legacy source; cannot replicate healthy versus malignant plasma-cell RNA comparison'}
    (HERE/'results/independent_assay_audit.json').write_text(json.dumps(adt_audit,indent=2))
    old=historical()
    report={'protocol_sha256':protocol_hash,'python':platform.python_version(),'numpy':np.__version__,'seed':20261002,'protected_files_unchanged':verify_inputs(),'bootstrap':intervals,'historical_metrics':old['metrics'],'raw_primary':[r for r in secondary if r['budget']==.1],'external_replication':'Not feasible from inspected local ADT source; mismatch explicitly retained'}
    (HERE/'results/extension_results.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
