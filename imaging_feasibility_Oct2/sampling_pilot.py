"""Prespecified retrospective sampling benchmark; no efficacy inference."""
import csv,gzip,hashlib,json,pathlib
import numpy as np
ROOT=pathlib.Path(__file__).resolve().parent

def profile(x,method):
    return x.mean(axis=0) if method=='mean' else np.quantile(x,[.1,.5,.9],axis=0).ravel()

def ranks(x):
    _,i,c=np.unique(x,return_inverse=True,return_counts=True)
    return (np.cumsum(c)-(c-1)/2)[i]

def main():
    chosen=json.loads((ROOT/'selection.json').read_text()); manifest=json.loads((ROOT/'download_manifest.json').read_text())
    plates=['BR00117011','BR00117012'];data={};finite=0;total=0
    for p in plates:
        for r in chosen:
            well=r['well_position'];data[p,well]=[]
            for s in range(1,10):
                path=ROOT/f'{p}_{well}_{s}.csv.gz'
                source=next(v for v in manifest if (v['plate'],v['well'],v['site'])==(p,well,s))
                assert hashlib.sha256(path.read_bytes()).hexdigest()==source['sha256_extract']
                with gzip.open(path,'rt') as f:
                    a=np.loadtxt(f,delimiter=',',skiprows=1,ndmin=2)[:,2:]
                assert len(a)==source['objects'];total+=len(a)
                mask=np.isfinite(a).all(axis=1);finite+=int(mask.sum());a=a[mask]
                assert len(a)>0 and a.shape[1]==8
                data[p,well].append(a)
    assert len(manifest)==216
    ctrl=[r['well_position'] for r in chosen if not r['broad_sample']]
    ref={p:np.concatenate([a for w in ctrl for a in data[p,w]]) for p in plates}
    scale=np.maximum(np.quantile(ref[plates[0]],.75,axis=0)-np.quantile(ref[plates[0]],.25,axis=0),1e-6)
    rows=[];rng=np.random.default_rng(20261002)
    for p in plates:
        for r in chosen:
            well=r['well_position'];sites=data[p,well];full=np.concatenate(sites)
            for method in ['mean','quantiles']:
                sc=scale if method=='mean' else np.tile(scale,3)
                rf=profile(ref[p],method);fp=profile(full,method)
                full_distance=float(np.abs((fp-rf)/sc).mean())
                for repeat in range(30):
                    chosen_sites=rng.choice(9,3,replace=False);sub=np.concatenate([sites[i] for i in chosen_sites])
                    random_cells=full[rng.choice(len(full),len(sub),replace=False)]
                    for sampling,a in [('three_sites',sub),('random_cells',random_cells)]:
                        dist=float(np.abs((profile(a,method)-rf)/sc).mean())
                        rows.append({'plate':p,'well':well,'condition':r['broad_sample'] or 'control',
                            'method':method,'sampling':sampling,'repeat':repeat,'cells':len(a),
                            'full_distance':full_distance,'sample_distance':dist,'absolute_error':abs(dist-full_distance)})
    with (ROOT/'sampling_results.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summaries=[]
    for method in ['mean','quantiles']:
        subset=[r for r in rows if r['method']==method]
        full=[r['full_distance'] for r in subset if r['repeat']==0 and r['sampling']=='three_sites']
        rms=float(np.sqrt(np.mean(np.square(full))))
        errors={}
        for sampling in ['three_sites','random_cells']:
            rr=[r for r in subset if r['sampling']==sampling]
            errors[sampling]={'median_absolute_error':float(np.median([r['absolute_error'] for r in rr])),
                'median_error_over_rms':float(np.median([r['absolute_error'] for r in rr])/rms)}
        distances={p:[next(r['full_distance'] for r in subset if r['plate']==p and r['well']==w['well_position']) for w in chosen] for p in plates}
        # Correlation is descriptive across 12 wells, only six non-control conditions.
        correlation=float(np.corrcoef(ranks(distances[plates[0]]),ranks(distances[plates[1]]))[0,1])
        ratio=errors['three_sites']['median_absolute_error']/errors['random_cells']['median_absolute_error']
        summaries.append({'method':method,'errors':errors,'site_to_cell_error_ratio':ratio,
            'field_structure_flag':bool(ratio>1.2),'plate_rank_correlation':correlation,'full_distances':distances})
    audit={'sites':len(manifest),'wells_per_plate':len(chosen),'plates':2,'treated_conditions':6,'control_wells_per_plate':6,
        'objects_total':total,'objects_finite':finite,'objects_removed':total-finite,
        'min_objects_per_site':min(len(a) for arrays in data.values() for a in arrays),
        'max_objects_per_site':max(len(a) for arrays in data.values() for a in arrays),
        'source_bytes':sum(v['source_bytes'] for v in manifest),'summaries':summaries,
        'protocol_sha256':hashlib.sha256((ROOT/'precommitment.md').read_bytes()).hexdigest(),
        'script_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()}
    (ROOT/'results.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit,indent=2))

if __name__=='__main__':main()
