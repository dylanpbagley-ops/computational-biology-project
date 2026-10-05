"""Serialize reviewed source observations and clearly labeled inherited screening."""
import csv,hashlib,json,pathlib,sys,collections
OUT=pathlib.Path(__file__).resolve().parents[1]; ROOT=OUT.parent
sys.path.insert(0,str(OUT));from startup import table,read
sys.path.insert(0,str(OUT/'analysis'));from thresholds import validate_count

SOURCE='source_id primary_url DOI_or_PMID publication_date access_date access_level inspected_locations prior_project_exposure search_id notes'.split()
STUDY='study_id cohort_id source_id gene antigen_isoform agent_id agent_name modality disease_context route combination_context treated_n cohort_definition observation_window overlap_group'.split()
OBS='observation_id study_id cohort_id source_id event_name normal_tissue event_count denominator denominator_unit grade injury_cessation death protocol_DLT attribution_code attribution_location depletion consequences supportive_care recovery ascertainment exposure_activity missingness_censoring observation_state adequacy_limits extractor review_status'.split()

def main():
    batch=json.loads((OUT/'evidence/fresh_batch.json').read_text())
    for name,fields in [('sources',SOURCE),('studies',STUDY),('observations',OBS)]:
        rs=[{k:r.get(k,'NR') for k in fields} for r in batch[name]]
        if name=='observations':
            for r in rs:
                if isinstance(r['event_count'],int):validate_count(r['event_count'],r['denominator'],'patients')
                if r['event_count']=='NR' and r['observation_state']=='explicit_nonobservation':raise ValueError('Missing is not zero')
        table(OUT/f'evidence/{name}.tsv',rs,fields)
    endpoint=OUT/'evidence/endpoint_v03.md'
    (OUT/'evidence/endpoint_v03.sha256').write_text(hashlib.sha256(endpoint.read_bytes()).hexdigest()+'\n')
    coverage={r['gene']:r for r in read(ROOT/'strengthening_Oct2/clinical_review_pass3/original_panel_coverage.tsv')}
    triage={r['gene']:r for r in read(ROOT/'strengthening_Oct2/negative_evidence_pass5/favorable_case_triage.tsv')}
    originals=read(ROOT/'strengthening_Oct2/operational_labels_44_Oct2.tsv')
    screened=[];review=[]
    for i,r in enumerate(originals):
        g=r['gene']; c=coverage.get(g)
        reason=triage[g]['unresolved_requirement'] if g in triage else (c['review_priority'] if c else 'Efficacy failure/exclusion does not establish injury absence; human event eligibility unreviewed')
        if g=='FAP':reason='Cited outcome is preclinical; ineligible evidence for human endpoint; no general gene verdict'
        if g=='PTPRC':reason='Old exclusion rationale contradicted by SIERRA human trial; conditioning scope needs review'
        screened.append(dict(candidate_id=f'OP{i+1:02}',source_id='inherited_source_inventory',study_id='not_extracted',gene=g,eligibility='ineligible_cited_evidence' if g=='FAP' else 'pending',reason=reason,prior_use_class='historical original outcome panel' if c else 'historical efficacy failure or excluded operational row',decision_date='2026-10-02',evidence_state='ineligible' if g=='FAP' else 'indeterminate_pending_event_review',source_urls=c['source_urls'] if c else 'not_extracted'))
        review.append(dict(record_id=f'OP{i+1:02}',gene=g,agent_identity=r['agent'],source_urls=c['source_urls'] if c else '',initial_eligibility='',initial_event='',initial_grade='',initial_attribution='',initial_negative_adequacy='',initial_judgment='',reviewer='',review_date='',source_location=''))
    table(OUT/'evidence/screening.tsv',screened)
    # This reviewer index has no outcomes, scores, expected organs or assistant judgments.
    table(OUT/'evidence/reviewer_packet/operational_index.tsv',review)
    event_review=[]
    sources={s['source_id']:s['primary_url'] for s in batch['sources']}
    for s in batch['studies']:
        event_review.append(dict(study_id=s['study_id'],cohort_id=s['cohort_id'],gene=s['gene'],agent_name=s['agent_name'],combination_context=s['combination_context'],primary_url=sources[s['source_id']],initial_event='',initial_grade='',initial_attribution='',initial_negative_adequacy='',initial_judgment='',reviewer='',review_date='',source_location=''))
    table(OUT/'evidence/reviewer_packet/batch_studies.tsv',event_review)
    overlaps=[dict(record_a='S-CD19-POOL',record_b='ELIANA original 29385370',relationship='pooled extension includes original trial',evidence='Methods; Table 1; references',overlap_group='ELIANA_ENSIGN',uncertainty='patient mapping not extracted; do not add denominators'),dict(record_a='S-CD33-FINAL',record_b='ALFA-0701 original 22482940',relationship='final update of same trial',evidence='Discussion cites initial cutoff',overlap_group='ALFA0701',uncertainty='as-treated sets differ; no silent substitution'),dict(record_a='E003',record_b='E004',relationship='final trial includes first three patients',evidence='CAIX Introduction explicitly describes earlier first three patients',overlap_group='CAIX12',uncertainty='not a second agent or independent cohort'),dict(record_a='CO-CA9-12',record_b='CO-CA9-3',relationship='different cohorts/regimens within same target/trial',evidence='CAIX abstract; regimen paragraph',overlap_group='CAIX12',uncertainty='separate patients but related development and same target')]
    table(OUT/'evidence/overlaps.tsv',overlaps)
    frame=read(OUT/'feasibility/candidate_a/pilot_frame.tsv')
    comparisons=[]
    notes={
      'B4GALNT1':('mapping_unresolved','GD2 is a glycolipid; biosynthesis gene is not directly measured antigen','different agents/diseases/routes; no matched intended-target endpoint'),
      'CA9':('known_context_case','same CAR with G250 pretreatment; earlier/final duplicate','4/8 cessation vs 0/4 liver nonobservation; window inadequate; already published intervention'),
      'CD19':('primary_trials_needed','7 regulatory product descriptions, not 7 trial cohorts','one pooled primary follow-up opened; other agents not endpoint-adjudicated'),
      'CD22':('same_agent_reports_and_other_agent','moxetumomab multiple studies and label; inotuzumab distinct ADC','different disease/transplant/payload contexts; mechanism and negatives unresolved'),
      'CD276':('context_and_attribution_unresolved','CAR/antibody/ADCs with different disease and route','severe immune/payload effects not equivalent to target-recognition injury; prior locked-list overlap'),
      'CD33':('clinical_attribution_unresolved','GO combination ADC; vadastuximab; single-patient CAR report','marrow injury cessation verified but target mechanism unresolved'),
      'CD44':('isoform_and_missing_report','CD44v6 ADC schedules; CAR source NOT_FOUND','gene-level RNA not isoform-matched; nonreported CAR outcomes cannot be negative'),
      'CD52':('same_agent_different_indication','alemtuzumab brands in CLL/MS; not distinct agents','cancer versus autoimmune indication; primary matched injury monitoring not extracted')}
    for g in sorted({r['gene_or_target_string'] for r in frame}):
        rs=[r for r in frame if r['gene_or_target_string']==g]
        state,identity,gap=notes.get(g,('control_unadjudicated','single-agent workload control','primary cohort/event extraction needed; no outcome metric'))
        comparisons.append(dict(gene=g,source_rows=len(rs),regulatory_rows=sum(r['pmid_or_url'].startswith('https://api.fda') for r in rs),distinct_descriptions=len({r['agent_as_extracted'] for r in rs}),independent_targets=1,usable_comparison='no independent adjudicated comparable binary pair',status=state,identity=identity,endpoint_tissue_grade_attribution_window_care_gap=gap))
    table(OUT/'feasibility/candidate_a/comparability.tsv',comparisons)
    inherited=read(ROOT/'strengthening_Oct2/clinical_review_pass2/clinical_evidence.tsv')+read(ROOT/'strengthening_Oct2/clinical_review_pass3/extraction.tsv')
    (OUT/'evidence/inherited_extractions.json').write_text(json.dumps({'status':'Inherited extraction; mixed event families not converted to incidence observations or fresh full-text review','records':inherited},indent=2))
    ledger=[]
    with (ROOT/'strengthening_Oct2/weakness_ledger.csv').open() as f:
        for r in csv.DictReader(f):
            ledger.append(dict(issue_id=r['id'],issue=r['issue'],inherited_status=r['status'],current_requirement=r['remaining_requirement'],execution_state='open_dependency' if r['id'] in ['W06','W07','W08','W09','W10','W12','W13','W14','W15','W16','W21','W22'] else 'preserve_existing_repair; no new closure',evidence='passes 2-5/current status; forms do not close adjudication'))
    table(OUT/'evidence/weakness_status.tsv',ledger)
    print(json.dumps({'screened_operational':len(screened),'original_panel':len(coverage),'fresh_event_rows':len(batch['observations']),'fresh_study_units':len(batch['studies']),'pilot_frame_rows':len(frame),'regulatory_rows':sum(r['pmid_or_url'].startswith('https://api.fda') for r in frame),'pilot_targets':len(comparisons),'amended_labels_created':False},indent=2))

if __name__=='__main__':main()
