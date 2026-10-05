"""Version evidence corrections and inventory prior source rows without scoring them."""
import csv,json,pathlib
from prepare import literal,tsv
ROOT=pathlib.Path(__file__).resolve().parent;BASE=ROOT.parent;SRC=BASE/'ev_capture_model/src/aim13_scoping'
targets=literal(SRC/'target_labels.py','TARGETS');discovery={t[0] for t in targets};validation={t[0] for t in literal(SRC/'target_labels_validation.py','VALIDATION')}
locked={r['gene'] for r in json.loads((SRC/'prospective_predictions.json').read_text())['predictions']}
corrections={
 'CEACAM5':('CEA-directed TCR-T (Parkhurst 2011)','Agent modality corrected; outcome still subject to consistent endpoint review','https://pubmed.ncbi.nlm.nih.gov/21157437/'),
 'FAP':('FAP-directed intervention; historical marrow/cachexia source is preclinical','Cited evidence cannot establish a human clinical-toxic outcome','https://pubmed.ncbi.nlm.nih.gov/23712432/'),
 'FOLR1':('mirvetuximab soravtansine','Recorded corneal on-target attribution contradicted; amended outcome unresolved','https://pubmed.ncbi.nlm.nih.gov/30413525/'),
 'PTPRC':('Anti-CD45 radioconjugate 131I-apamistamab exists in human trials','Historical no-human-agent exclusion reason contradicted; transplant-conditioning eligibility unresolved','https://pubmed.ncbi.nlm.nih.gov/39298738/')}
records=[]
for t in targets:
    correction=corrections.get(t[0],(t[2],'No final source adjudication in this pass',''))
    records.append({'gene':t[0],'historical_agent':t[2],'historical_outcome':t[3],'amended_agent_or_evidence_identity':correction[0],'source_correction':correction[1],'verified_source_link':correction[2],'adjudicated_outcome':'PENDING_REVIEW','status':'Evidence amendment only; not a replacement scoring-label table'})
tsv(ROOT/'evidence_amendments_Oct2.tsv',records,list(records[0]))
index=[]
for file in ['trial_agents.tsv','trial_agents_2.tsv','approved_agents.tsv','approved_agents_2.tsv']:
    rows=list(csv.DictReader((BASE/'ev_capture_data/round7_agent_data'/file).open(encoding='utf-8-sig'),delimiter='\t'))
    for line,r in enumerate(rows,2):
        gene=r['gene'];attribution=r.get('on_target_attribution',r.get('attribution_stated',''))
        index.append({'evidence_id':f'E{len(index)+1:03d}','gene_or_target_string':gene,'agent_as_extracted':r['agent'],'modality_as_extracted':r.get('modality',''),'source_file':file,'source_line':line,'pmid_or_url':r.get('pmid',r.get('source_url','')),'prior_source_access_flag':r.get('source_type',''),'attribution_as_extracted':attribution,'already_in_original_table':gene in discovery,'already_in_validation_table':gene in validation,'already_in_locked_predictions':gene in locked,'identity_issue':'dual target; do not treat as one gene' if '/' in gene else 'biosynthesis gene is not the GD2 surface antigen; review mapping' if gene=='B4GALNT1' else 'check isoform and antigen identity','eligible_endpoint':'PENDING_REVIEW','independent_validation_status':'Not established; previously AI-extracted source pool'})
tsv(ROOT/'agent_evidence_index.tsv',index,list(index[0]));print(json.dumps({'evidence_amendment_rows':len(records),'prior_source_index_rows':len(index),'new_outcomes_scored':0}))
