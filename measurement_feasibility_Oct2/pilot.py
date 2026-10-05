"""Exploratory healthy-cell annotation benchmark; numpy-only, no therapy modeling."""
import csv, gzip, hashlib, json, pathlib, time
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent
DATA = ROOT.parent / 'ev_capture_data/healthy_pbmc_citeseq'

def auc(y, x):
    _, inv, cnt = np.unique(x, return_inverse=True, return_counts=True)
    rank = np.cumsum(cnt) - (cnt-1)/2
    p = y.sum(); q = len(y)-p
    return float((rank[inv][y].sum()-p*(p+1)/2)/(p*q))

def fit(x, y, w):
    mean = np.average(x, axis=0, weights=w)
    sd = np.sqrt(np.average((x-mean)**2, axis=0, weights=w))
    sd = np.maximum(sd, 1e-5)
    a = np.column_stack([np.ones(len(x)), (x-mean)/sd])
    w = w/w.sum(); b = np.zeros(a.shape[1]); reg = np.eye(len(b))*.01; reg[0,0]=0
    for _ in range(60):
        p = 1/(1+np.exp(-np.clip(a@b,-35,35)))
        g = a.T@(w*(p-y))+reg@b
        h = a.T@((w*p*(1-p))[:,None]*a)+reg+np.eye(len(b))*1e-9
        step = np.linalg.solve(h,g); b -= step
        if np.max(np.abs(step))<1e-7: break
    assert np.isfinite(b).all()
    return mean, sd, b

def predict(model,x):
    mean,sd,b=model
    return b[0]+((x-mean)/sd)@b[1:]

def run():
    assert auc(np.array([False,True]),np.array([0,1])) == 1
    assert auc(np.array([False,True]),np.array([1,1])) == .5
    xx=np.arange(20.)[:,None]; yy=xx[:,0]>9
    mm=fit(xx,yy,np.ones(20)); assert auc(yy,predict(mm,xx))==1
    z=np.load(DATA/'healthy_adt.npz',allow_pickle=False)
    raw=z['adt']; names=z['names']
    with gzip.open(ROOT.parent/'direction_landscape_Oct1/pbmc_features.tsv.gz','rt') as f:
        official=np.array([l.strip().split('\t')[0] for l in f])
    assert np.array_equal(names,official), 'Feature order mismatch'
    with gzip.open(DATA/'GSE164378_sc.meta.data_3P.csv.gz','rt') as f: meta=list(csv.DictReader(f))
    donor=np.array([r['donor'] for r in meta]);label=np.array([r['celltype.l1'] for r in meta]);l2=np.array([r['celltype.l2'] for r in meta])
    with gzip.open(DATA/'GSM5008738_ADT_3P-barcodes.tsv.gz','rt') as f: barcodes=[l.strip() for l in f]
    assert barcodes==[r[''] for r in meta]
    assert np.array_equal(donor,[r['donor'] for r in meta])
    assert np.array_equal(label,[r['celltype.l1'] for r in meta])
    assert np.isfinite(raw).all() and raw.min()>=0
    checked=0; repairs=[]
    with gzip.open(DATA/'GSM5008738_ADT_3P-matrix.mtx.gz','rt') as f:
        for line in f:
            if not line.startswith('%'): dims=list(map(int,line.split()));break
        assert dims[:2]==[raw.shape[1],raw.shape[0]]
        for i,line in enumerate(f):
            if i%10007==0:
                r,c,v=line.split();rr=int(r)-1;cc=int(c)-1
                if raw[cc,rr]!=float(v):
                    assert i==0 and rr==0 and cc==0 and raw[cc,rr]==0 and float(v)==66
                    repairs.append({'row':rr,'column':cc,'old':0,'new':66,'reason':'legacy loader skiprows=3 skipped first entry'})
                    raw[cc,rr]=float(v)
                checked+=1
    audit={'shape':list(raw.shape),'feature_order_verified':True,'barcode_order_verified':True,
           'sampled_raw_nonzero_coordinates_verified':checked,'raw_coordinates_seen':i+1,
           'label_independence':False,'numpy':np.__version__,'working_copy_repairs':repairs}
    (ROOT/'data_audit.json').write_text(json.dumps(audit,indent=2))
    good=(label!='other')&(l2!='Doublet')
    raw=np.log1p(raw[good]); donor=donor[good]; label=label[good]
    candidates=np.array([i for i,n in enumerate(names) if 'IgG' not in n and 'IgG' not in n.replace('Rag','Rat')])
    donors=np.unique(donor); rows=[]; counts={}
    for task in ['DC','B']:
        y=label==task
        counts[task]={d:{'positive':int(y[donor==d].sum()),'total':int((donor==d).sum())} for d in donors}
        for fold,test in enumerate(donors):
            cal=donors[(fold+1)%len(donors)]; train=donors[(donors!=test)&(donors!=cal)]
            assert test not in train and cal not in train
            for budget in [250,1000]:
                indices=[];weights=[]
                for di,d in enumerate(train):
                    ix=np.flatnonzero(donor==d); rng=np.random.default_rng(20261002+fold*100+di)
                    ix=rng.choice(ix,min(budget,len(ix)),replace=False)
                    pos=np.flatnonzero((donor==d)&y)
                    missing=max(0,min(50,len(pos))-int(y[ix].sum()))
                    if missing: ix=np.union1d(ix,rng.choice(np.setdiff1d(pos,ix),missing,replace=False))
                    p_all=y[donor==d].mean(); p_sample=y[ix].mean()
                    w=np.where(y[ix],p_all/p_sample,(1-p_all)/(1-p_sample))/len(ix)
                    indices.extend(ix);weights.extend(w)
                ix=np.array(indices); w=np.array(weights); x=raw[ix][:,candidates]; yt=y[ix]
                mp=np.average(x[yt],axis=0,weights=w[yt]);mn=np.average(x[~yt],axis=0,weights=w[~yt])
                vp=np.average((x[yt]-mp)**2,axis=0,weights=w[yt]);vn=np.average((x[~yt]-mn)**2,axis=0,weights=w[~yt])
                direction=np.where(mp>=mn,1,-1)
                aa=np.array([[auc(y[ix[donor[ix]==d]],raw[ix[donor[ix]==d],j]*direction[k]) for k,j in enumerate(candidates)] for d in train])
                scores={'mean_effect':np.abs(mp-mn)/np.sqrt((vp+vn)/2+1e-8),
                        'mean_auc':aa.mean(axis=0),'worst_auc':aa.min(axis=0)}
                for method,score in scores.items():
                    for k in [1,2,4]:
                        selected=candidates[np.argsort(-score,kind='stable')[:k]]
                        model=fit(raw[ix][:,selected],yt,w)
                        ci=np.flatnonzero(donor==cal);ti=np.flatnonzero(donor==test)
                        cp=predict(model,raw[ci][:,selected]); threshold=np.quantile(cp[~y[ci]],.99,method='higher')
                        pred=predict(model,raw[ti][:,selected])>threshold; truth=y[ti]
                        tp=int((pred&truth).sum());fp=int((pred&~truth).sum()); p=int(truth.sum());n=len(truth)-p
                        rows.append({'task':task,'test_donor':test,'calibration_donor':cal,'budget':budget,
                            'method':method,'features':k,'selected':','.join(names[selected]),'train_n':len(ix),
                            'tp':tp,'fp':fp,'positive':p,'negative':n,'recall':tp/p,'fpr':fp/n,
                            'precision':tp/(tp+fp) if tp+fp else None})
            print('Completed',task,test,flush=True)
    with (ROOT/'results.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    summary=[]
    for task in ['DC','B']:
        for budget in [250,1000]:
            for k in [1,2,4]:
                for method in ['mean_effect','mean_auc','worst_auc']:
                    r=[v for v in rows if (v['task'],v['budget'],v['features'],v['method'])==(task,budget,k,method)]
                    summary.append({'task':task,'budget':budget,'features':k,'method':method,
                        'median_recall':float(np.median([v['recall'] for v in r])),
                        'min_recall':min(v['recall'] for v in r),'median_fpr':float(np.median([v['fpr'] for v in r])),
                        'fpr_pass_donors':sum(v['fpr']<=.015 for v in r),
                        'median_precision':float(np.median([v['precision'] or 0 for v in r]))})
    primary={m:[v for v in rows if (v['task'],v['budget'],v['features'],v['method'])==('DC',1000,4,m)] for m in ['mean_effect','mean_auc','worst_auc']}
    contrast={}
    for baseline in ['mean_effect','mean_auc']:
        delta=[a['recall']-b['recall'] for a,b in zip(primary['worst_auc'],primary[baseline])]
        pa=sum(a['fpr']<=.015 for a in primary['worst_auc']);pb=sum(a['fpr']<=.015 for a in primary[baseline])
        contrast[baseline]={'paired_recall_differences':delta,'median_difference':float(np.median(delta)),
            'positive_donors':sum(v>0 for v in delta),'robust_fpr_pass':pa,'baseline_fpr_pass':pb,
            'screen_pass':bool(np.median(delta)>=.05 and sum(v>0 for v in delta)>=6 and pa>=6 and pa>=pb)}
    out={'counts':counts,'summary':summary,'primary_contrasts':contrast,
         'protocol_sha256':hashlib.sha256((ROOT/'precommitment.md').read_bytes()).hexdigest(),
         'script_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()}
    (ROOT/'results.json').write_text(json.dumps(out,indent=2)); print(json.dumps(contrast,indent=2),flush=True)

if __name__=='__main__':run()
