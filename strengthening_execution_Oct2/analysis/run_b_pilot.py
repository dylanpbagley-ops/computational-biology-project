"""Run the frozen B pilot on matched real donor measurements, preserving inputs."""
import csv,gzip,hashlib,json,pathlib,sys,collections,platform
import numpy as np
from thresholds import select_threshold, fraction_above, auc
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]))
from startup import table

HERE=pathlib.Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OUT=HERE.parent
GENES=['TNFRSF17','GPRC5D','SLAMF7']
POLICIES=['pooled','equal_donor','donor_envelope']

def verify_inputs():
    recorded=json.loads((OUT/'input_hashes.json').read_text())
    changed=[p for p,h in recorded.items() if hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=h]
    if changed: raise RuntimeError('Protected input changed: '+str(changed))
    return len(recorded)

def load(scale='normalized'):
    folder=ROOT/'ev_capture_data/healthy_bone_marrow_boiarsky'
    d=json.loads((folder/'extracted_markers.json').read_text())
    with gzip.open(folder/'GSE193531_cell-level-metadata.csv.gz','rt') as f: meta=list(csv.DictReader(f))
    if [r['index'] for r in meta]!=d['cell_ids']:raise ValueError('Cell order mismatch')
    counts=np.array([float(r['n_counts']) for r in meta])
    if not np.isfinite(counts).all() or (counts<=0).any():raise ValueError('Invalid library totals')
    groups={'MM':{},'NBM':{}}
    for stage,kind in [('MM','neoplastic'),('NBM','normal')]:
        for sample in sorted({r['sample_ID'] for r in meta if r['disease_stage']==stage}):
            ix=np.array([i for i,r in enumerate(meta) if r['sample_ID']==sample and r['normal_or_neoplastic']==kind])
            if not len(ix):raise ValueError('Empty eligible sample')
            groups[stage][sample]={g:(np.array(d['expr'][g],float)[ix]/counts[ix]*10000 if scale=='normalized' else np.array(d['expr'][g],float)[ix]) for g in GENES}
    return groups

def main():
    protected=verify_inputs()
    spec=OUT/'feasibility/candidate_b/pilot_specification.md'
    spec_hash=hashlib.sha256(spec.read_bytes()).hexdigest()
    lock=OUT/'feasibility/candidate_b/pilot_specification.sha256'
    if lock.exists() and lock.read_text().strip()!=spec_hash:raise RuntimeError('Specification changed')
    lock.write_text(spec_hash+'\n')
    groups=load(); mm=list(groups['MM']); nbm=list(groups['NBM'])
    summaries=[]; donor=[]; folds=[]; pairwise=[]
    for g in GENES:
        ts=[groups['MM'][s][g] for s in mm];ns=[groups['NBM'][s][g] for s in nbm]
        for stage,samples in groups.items():
            for s,v in samples.items():
                a=v[g]; donor.append(dict(gene=g,sample_id=s,stage=stage,cells=len(a),mean_UMI_per_10000=float(a.mean()),detected_fraction=float((a>0).mean())))
        for t in mm:
            for n in nbm:pairwise.append(dict(gene=g,tumor_sample=t,normal_sample=n,auc=auc(groups['MM'][t][g],groups['NBM'][n][g])))
        ratio=float(np.mean([a.mean() for a in ts])/np.mean([a.mean() for a in ns]))
        for budget in [.05,.10,.20,.50]:
            for policy in POLICIES:
                r=select_threshold(ts,ns,budget,policy)
                summaries.append(dict(gene=g,budget=budget,policy=policy,**r,equal_donor_mean_ratio=ratio))
        # Candidate values are learned from training data only, inside each fold.
        for held_t in mm:
            for held_n in nbm:
                train_t=[s for s in mm if s!=held_t];train_n=[s for s in nbm if s!=held_n]
                assert held_t not in train_t and held_n not in train_n and not set(train_t)&set(train_n)
                for policy in POLICIES:
                    r=select_threshold([groups['MM'][s][g] for s in train_t],[groups['NBM'][s][g] for s in train_n],.10,policy)
                    normal=fraction_above(groups['NBM'][held_n][g],r['threshold'])
                    folds.append(dict(gene=g,held_tumor=held_t,held_normal=held_n,policy=policy,threshold=r['threshold'],normal_fraction=normal,tumor_fraction=fraction_above(groups['MM'][held_t][g],r['threshold']),normal_budget_exceeded=normal>.10+1e-14,train_tumor_ids=','.join(train_t),train_normal_ids=','.join(train_n)))
    table(HERE/'results/b_summary.tsv',summaries)
    table(HERE/'results/b_donors.tsv',donor)
    table(HERE/'results/b_pairwise_auc.tsv',pairwise)
    table(HERE/'results/b_held_out.tsv',folds)
    metrics=[]
    for g in GENES:
        for policy in POLICIES:
            rs=[r for r in folds if r['gene']==g and r['policy']==policy]
            metrics.append(dict(gene=g,policy=policy,folds=len(rs),distinct_MM=8,distinct_NBM=9,exceeded_folds=sum(r['normal_budget_exceeded'] for r in rs),normal_mean=float(np.mean([r['normal_fraction'] for r in rs])),normal_worst=max(r['normal_fraction'] for r in rs),tumor_mean=float(np.mean([r['tumor_fraction'] for r in rs])),tumor_worst=min(r['tumor_fraction'] for r in rs)))
    table(HERE/'results/b_held_out_summary.tsv',metrics)
    result={'protected_files_unchanged':verify_inputs(),'specification_sha256':spec_hash,'numpy':np.__version__,'python':platform.python_version(),'NBM_samples':nbm,'MM_samples':mm,'primary_summary':[r for r in summaries if r['budget']==.1],'held_out_summary':metrics,'interpretation':'Exact empirical deterministic thresholds on RNA; 72 dependent folds per gene/policy, 17 biological samples; no external or clinical validation'}
    (HERE/'results/b_pilot.json').write_text(json.dumps(result,indent=2,allow_nan=False))
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
