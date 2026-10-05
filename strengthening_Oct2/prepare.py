"""Preserve historical artifacts and prepare score-free review inputs."""
import ast,csv,hashlib,json,pathlib,random,re,shutil
ROOT=pathlib.Path(__file__).resolve().parent; BASE=ROOT.parent
SRC=BASE/'ev_capture_model/src/aim13_scoping'; DOC=BASE/'ev_capture_model/docs'
def literal(path,name):
    return next(ast.literal_eval(n.value) for n in ast.parse(path.read_text(encoding='utf-8-sig')).body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id==name for t in n.targets))
def tsv(path,rows,fields):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
def main():
    protocol=ROOT/'analysis_precommitment.md'
    (ROOT/'analysis_precommitment.md.hashes').write_text(hashlib.sha256(protocol.read_bytes()).hexdigest()+'  analysis_precommitment.md\n')
    files=[SRC/n for n in ['target_labels.py','target_labels_frozen.tsv','target_labels_validation.py','prospective_predictions.json','wholebody_v2.json','flag_cleanup.json','weighted_ceiling_theorem.py','weighted_ceiling_theorem.json']]
    files+=list(SRC.glob('*.hashes'))
    files += [BASE/'HANDOFF_START_HERE_Oct1.md',BASE/'Research_Direction_Reset_Oct2.md',DOC/'Aim13_Results_Two_Comparator_Ceiling.md',DOC/'ISEF_Paper_Skeleton.md',DOC/'Label_Review_Sheet.md',BASE/'Poster, Paper, and Write Up/ISEF_Paper_Skeleton.md']
    manifest={}
    for f in files:
        rel=f.relative_to(BASE);dest=ROOT/'before'/rel
        if dest.exists():raise RuntimeError('Snapshot exists: refusing to replace '+str(dest))
        dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(f,dest)
        manifest[str(rel)]=hashlib.sha256(f.read_bytes()).hexdigest()
    (ROOT/'before_hashes.json').write_text(json.dumps(manifest,indent=2))
    targets=literal(SRC/'target_labels.py','TARGETS'); labels={t[0]:t[3] for t in targets}
    fields=['gene','common','agent','outcome','evidence','expected_normal_site']
    tsv(ROOT/'operational_labels_44_Oct2.tsv',[dict(zip(fields,t)) for t in targets],fields)
    old=list(csv.DictReader((SRC/'target_labels_frozen.tsv').open(encoding='utf-8-sig'),delimiter='\t'))
    added=sorted(set(labels)-{r['gene'] for r in old})
    changed=[r['gene'] for r in old if r['outcome']!=labels[r['gene']]]
    pool=[]
    for name in ['trial_agents.tsv','trial_agents_2.tsv','approved_agents.tsv','approved_agents_2.tsv']:
        rows=list(csv.DictReader((BASE/'ev_capture_data/round7_agent_data'/name).open(encoding='utf-8-sig'),delimiter='\t'))
        for row in rows: row['source_file']=name
        pool+=rows
    sources={}
    for r in pool:
        links=[]
        for pmid in re.findall(r'\b\d{7,8}\b',r.get('pmid','')):links.append('https://pubmed.ncbi.nlm.nih.gov/'+pmid+'/')
        if r.get('source_url','').startswith('http'):links.append(r['source_url'])
        sources.setdefault(r['gene'],set()).update(links)
    selected=[t for t in targets if t[3] in ('safe','toxic')]; random.Random(20261002).shuffle(selected)
    review=[]; key=[]
    for i,t in enumerate(selected,1):
        rid=f'R{i:02d}';key.append({'review_id':rid,'gene':t[0],'historical_outcome':t[3]})
        review.append({'review_id':rid,'gene':t[0],'historical_agent_text_to_reconcile':t[2],'candidate_source_links':' | '.join(sorted(sources.get(t[0],set()))),'exact_agent_and_study':'','human_evidence':'','normal_tissue_event':'','severity_and_dose_limitation':'','attribution_and_uncertainty':'','adequacy_of_negative_evidence':'','eligible_endpoint':'','proposed_category':'','reviewer_reason_and_source_location':'','reviewer_name_date':''})
    tsv(ROOT/'reviewer_form_score_free.tsv',review,list(review[0]))
    tsv(ROOT/'review_key_do_not_send.tsv',key,list(key[0]))
    v=literal(SRC/'target_labels_validation.py','VALIDATION');vg={t[0] for t in v};genes={r['gene'] for r in pool}
    inventory=[{'gene':g,'evidence_rows':sum(r['gene']==g for r in pool),'previous_use':'discovery_or_original_exclusion' if g in labels else 'validation_or_hedged_or_excluded' if g in vg else 'previously_extracted_not_in_those_tables','source_links':' | '.join(sorted(sources.get(g,set())))} for g in sorted(genes)]
    tsv(ROOT/'existing_evidence_pool.tsv',inventory,list(inventory[0]))
    out={'original_frozen_rows':len(old),'operational_rows':len(targets),'omitted_genes':added,'shared_outcomes_changed':changed,'review_rows':len(review),'pool_rows':len(pool),'pool_unique_genes':len(genes),'pool_unique_gene_agent_strings':len({(r['gene'],r['agent']) for r in pool}),'pool_genes_outside_original_and_validation':sorted(genes-set(labels)-vg),'status':'Source inventory only; no new outcome eligibility adjudicated; pool has previous AI exposure'}
    (ROOT/'preparation_results.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':main()
