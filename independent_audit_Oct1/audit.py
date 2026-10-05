"""Independent fixed-rule audit. Does not modify original project artifacts."""
from pathlib import Path
import sys, json, csv, hashlib, itertools
import numpy as np
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent
SRC = ROOT / 'ev_capture_model/src/aim13_scoping'
sys.path.insert(0, str(SRC))
import wholebody_v2 as WB
import target_labels as TL
import lineage_map as LM
import weighted_ceiling_theorem as WC
import adc_autopsy_ceiling as AC

def j(name):
    return json.loads((SRC/name).read_text(encoding='utf-8'))

def auc(a,b):
    d=np.asarray(a)[:,None]-np.asarray(b)[None,:]
    return float(((d>0)+.5*(d==0)).mean())

def boot_pair(sa,ta,sb=None,tb=None):
    rng=np.random.default_rng(20261001)
    sa,ta=np.asarray(sa),np.asarray(ta)
    ia=rng.integers(len(sa),size=(20000,len(sa)))
    ib=rng.integers(len(ta),size=(20000,len(ta)))
    delta=sa[ia][:,:,None]-ta[ib][:,None,:]
    v=((delta>0)+.5*(delta==0)).mean(axis=(1,2))
    res={'auc':auc(sa,ta),'ci95':np.quantile(v,[.025,.975]).tolist()}
    if sb is not None:
        sb,tb=np.asarray(sb),np.asarray(tb)
        d=sb[ia][:,:,None]-tb[ib][:,None,:]
        v2=((d>0)+.5*(d==0)).mean(axis=(1,2))
        res['paired_difference']=auc(sa,ta)-auc(sb,tb)
        res['paired_difference_ci95']=np.quantile(v-v2,[.025,.975]).tolist()
    return res

def upper_hull(points):
    h=[]
    for p in sorted(set(points)):
        while len(h)>1:
            a,b=h[-2:]
            cross=(b[0]-a[0])*(p[1]-a[1])-(b[1]-a[1])*(p[0]-a[0])
            if cross>=-1e-12: h.pop()
            else: break
        h.append(p)
    return h

def threshold_points(t,n):
    return [(float(sum(n[k:])),float(sum(t[k:]))) for k in range(len(t)+1)]

def math_checks():
    t=np.array([.45,.1,.45]); n=np.array([.1,.8,.1])
    points=threshold_points(t,n); h=upper_hull(points)
    allrules=[np.array(v) for v in itertools.product([0,1],repeat=3)]
    arb=max(float(f@(t-n)) for f in allrules)
    mono=max(y-x for x,y in points)
    res={'arbitrary_vs_monotone':{'pT':t.tolist(),'pN':n.tolist(),'antigen_threshold_hull':h,
        'arbitrary_rule':[1,0,1],'arbitrary_point':[.2,.9],
        'monotone_hull_ceiling_at_normal_kill_0.2':float(np.interp(.2,*np.array(h).T)),
        'best_arbitrary_T_minus_H':arb,'best_monotone_T_minus_H':mono}}
    n=np.full(4,.25); a=n.copy(); b=np.array([.45,.05,.05,.45])
    def ranked_auc(t,n): return float(sum(t[i]*(sum(n[:i])+.5*n[i]) for i in range(len(t))))
    ha,hb=upper_hull(threshold_points(a,n)),upper_hull(threshold_points(b,n))
    res['same_auc_different_ceiling']={'auc_A':ranked_auc(a,n),'auc_B':ranked_auc(b,n),
        'hull_A':ha,'hull_B':hb,'normal_kill_0.25_tumor_kill_A':float(np.interp(.25,*np.array(ha).T)),
        'normal_kill_0.25_tumor_kill_B':float(np.interp(.25,*np.array(hb).T))}
    res['existing_weighted_verifier']=dict(zip(['worst_excess','worst_np_gap','checks'],WC.verify()))
    return res

def main():
    stored=j('wholebody_v2.json'); flags=j('flag_cleanup.json')['targets']
    labels={r[0]:r[3] for r in TL.TARGETS}
    vals,_=WB.load('rna_single_cell_type.tsv.zip',['nCPM','nTPM','pTPM'],['Cell type'])
    raw=WB.confinement(vals)
    genes=sorted(flags); safe=[g for g in genes if labels[g]=='safe']; toxic=[g for g in genes if labels[g]=='toxic']
    score={}; discrepancies=[]
    for g in genes:
        d=vals[g]; fam={}
        for c,v in d.items(): fam[LM.lineage(c)]=fam.get(LM.lineage(c),0)+v
        flag=max(fam,key=fam.get)=='blood'
        score[g]=raw[g]['top_share']+float(flag)
        if flag!=flags[g]['M2_flag']: discrepancies.append(g)
    result={'raw_hpa':{'resolved_genes':len(raw),'cell_type_counts':sorted(set(map(len,vals.values()))),
        'max_top_share_difference':max(abs(raw[g]['top_share']-stored['reference_A_cell_type'][g]['top_share']) for g in raw),
        'flag_discrepancies':discrepancies},'operational_counts':{lab:sum(x==lab for x in labels.values()) for lab in sorted(set(labels.values()))}}
    s=[score[g] for g in safe]; t=[score[g] for g in toxic]
    c_s=[raw[g]['top_share'] for g in safe]; c_t=[raw[g]['top_share'] for g in toxic]
    result['main']=boot_pair(s,t,c_s,c_t)
    result['main']['exact_p']=float(mannwhitneyu(s,t,alternative='greater',method='exact').pvalue)
    result['confinement']=boot_pair(c_s,c_t)
    result['confinement']['exact_p']=float(mannwhitneyu(c_s,c_t,alternative='greater',method='exact').pvalue)
    result['pairs']={'concordant':sum(a>b for a in s for b in t),'discordant':sum(a<b for a in s for b in t),'tied':sum(a==b for a in s for b in t)}
    rows=list(csv.DictReader((SRC/'target_labels_frozen.tsv').open(encoding='utf-8'),delimiter='\t'))
    old={r['gene']:r['outcome'] for r in rows}
    result['frozen_tsv']={'counts':{v:sum(x==v for x in old.values()) for v in sorted(set(old.values()))},
        'missing_operational_genes':sorted(set(labels)-set(old)),
        'disagreeing_labels':[g for g in old if g in labels and old[g]!=labels[g]]}
    ceil=j('adc_autopsy_ceiling.json')['within']
    _,_,patients=AC.load()
    recomputed={}
    for g in ceil:
        per=[]
        for profile in patients.values():
            a,b=profile.get(g,([],[]))
            if len(a)>=3 and len(b)>=3:
                per.append(float(mannwhitneyu(a,b).statistic/(len(a)*len(b))))
        if len(per)>=5: recomputed[g]=float(np.mean(per))
    result['raw_autopsy_comparator']={'genes_recomputed':len(recomputed),
        'max_difference_from_stored':max(abs(recomputed[g]-ceil[g]) for g in ceil)}
    shared=sorted(g for g in ceil if g in raw and labels.get(g) in ('safe','toxic'))
    ss=[g for g in shared if labels[g]=='safe']; tt=[g for g in shared if labels[g]=='toxic']
    result['head_to_head']=boot_pair([raw[g]['top_share'] for g in ss],[raw[g]['top_share'] for g in tt],
        [ceil[g] for g in ss],[ceil[g] for g in tt])
    result['head_to_head'].update(n_safe=len(ss),n_toxic=len(tt),genes=shared,
        comparator_auc=auc([ceil[g] for g in ss],[ceil[g] for g in tt]))
    hashes=[]
    for hf in ['label_fragility_precommitment.md.hashes','locked_predictions_check_precommitment.md.hashes']:
        for line in (SRC/hf).read_text().splitlines():
            h,name=line.split(); hashes.append({'file':name,'matches':hashlib.sha256((SRC/name).read_bytes()).hexdigest()==h})
    result['hashes']=hashes
    result['mathematics']=math_checks()
    result['notes']=['Bootstrap uncertainty conditional on fixed labels and panel; it does not include selection or label uncertainty.',
        'Head-to-head compares existing predictors on existing labels; no independent safety replication.',
        'No original labels, predictions, rules or outputs changed.']
    (OUT/'audit_results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
