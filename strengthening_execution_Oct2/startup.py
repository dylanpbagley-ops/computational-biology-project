"""Non-destructive startup inventory. Run once; never replace the initial hashes."""
import csv, gzip, hashlib, json, pathlib, collections, datetime

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = pathlib.Path(__file__).resolve().parent

def table(path, rows, fields=None):
    rows = list(rows)
    if fields is None: fields = list(rows[0])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader(); w.writerows(rows)

def read(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def main():
    src = ROOT/'ev_capture_model/src/aim13_scoping'
    protected = [src/n for n in ['target_labels.py','target_labels_frozen.tsv','wholebody_v2.json','flag_cleanup.json','prospective_predictions.json']]
    protected += [ROOT/'strengthening_Oct2/operational_labels_44_Oct2.tsv']
    for folder in ['clinical_review_pass2','clinical_review_pass3','uncertainty_pass4','negative_evidence_pass5']:
        protected.extend(p for p in (ROOT/'strengthening_Oct2'/folder).iterdir() if p.is_file())
    protected += [ROOT/'strengthening_Oct2'/n for n in ['agent_evidence_index.tsv','evidence_pool_with_prediction_overlap.tsv','diagnostic_results.json','weakness_ledger.csv']]
    protected += list((ROOT/'ev_capture_data/healthy_bone_marrow_boiarsky').glob('*'))
    manifest = OUT/'input_hashes.json'
    if manifest.exists(): raise RuntimeError('Initial manifest exists: resume rather than overwrite startup')
    manifest.write_text(json.dumps({str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in protected if p.is_file()},indent=2),encoding='utf-8')
    rows = read(ROOT/'strengthening_Oct2/agent_evidence_index.tsv')
    groups = collections.defaultdict(list)
    for r in rows: groups[r['gene_or_target_string']].append(r)
    # Candidate labels are leads, not established distinct agents or cohorts.
    multi = sorted(g for g, rs in groups.items() if len({r['agent_as_extracted'] for r in rs}) >= 2 and '/' not in g and '+' not in g)[:8]
    controls=[]; modalities=set()
    for g in sorted(groups):
        rs=groups[g]
        if len({r['agent_as_extracted'] for r in rs})!=1 or '/' in g or '+' in g: continue
        m=rs[0]['modality_as_extracted']
        if m not in modalities: controls.append(g); modalities.add(m)
        if len(controls)==4: break
    selected=[dict(r, pilot_group='multi_agent_lead' if g in multi else 'single_agent_control') for g in multi+controls for r in groups[g]]
    table(OUT/'feasibility/candidate_a/pilot_frame.tsv',selected)
    (OUT/'feasibility/pilot_selection.md').write_text(f'''# Pilot selection locked before new performance calculations
October 2, 2026; post-discovery development data, not preregistered validation.

Candidate A: first eight alphabetically ordered single-target strings with at least two distinct extracted agent/study descriptions. All source rows for those strings are retained. These descriptions may be duplicate reports, not distinct agents. Dual strings containing / or + are not assigned a single-gene outcome. B4GALNT1/GD2 remains a flagged biosynthesis proxy until resolved.

Multi-agent leads: {', '.join(multi)}.
Controls: first alphabetical single-agent-description target from each of four distinct extracted modality groups: {', '.join(controls)}. Modality labels are inherited and need normalization. No outcome or score used for selection. This is a workload rule, not a sample-size adequacy test.

Candidate B: inspect existing Boiarsky GSE193531 metadata and 11-gene count extraction. Candidate primary comparison is malignant MM plasma cells versus healthy NBM plasma cells for TNFRSF17, GPRC5D and SLAMF7, established myeloma antigens available in the extraction. Dominant proposed uncertainty is finite biological donor sampling. No RNA-to-protein conversion, inferred harm weight, reconstructed joint distribution or absolute calibration. Establish sample identities, normalization and closest prior-art overlap before calculations. One donor is one biological unit; cells estimate its distribution. A leave-donor-out check is internal evaluation, not an independent dataset.

No candidate performance metric is authorized by this frame alone. Clinical amended labels require their gate; B requires a specified estimand and assumptions. Deadline: late January school fair, mid-March state fair; weekly effort unknown/unfixed (user response).
''',encoding='utf-8')
    p=ROOT/'ev_capture_data/healthy_bone_marrow_boiarsky'
    with gzip.open(p/'GSE193531_cell-level-metadata.csv.gz','rt') as f: meta=list(csv.DictReader(f))
    d=json.loads((p/'extracted_markers.json').read_text())
    assert [r['index'] for r in meta]==d['cell_ids']
    audit=[]
    for (sample,stage,kind), n in sorted(collections.Counter((r['sample_ID'],r['disease_stage'],r['normal_or_neoplastic']) for r in meta).items()):
        audit.append(dict(sample_id=sample,disease_stage=stage,normal_or_neoplastic=kind,cells=n,unit='sample; donor identity to confirm in primary methods',measurement='raw RNA UMI counts',normalization='not yet applied; n_counts available',detection_limit='zero means not detected, not biological absence',scope='CD138-selected plasma cells; no whole-body safety inference'))
    table(OUT/'feasibility/candidate_b/sample_audit.tsv',audit)
    table(OUT/'decision_log.tsv',[dict(date='2026-10-02',decision='Execute in existing non-Git workspace; preserve durable ledger',evidence='No .git; user specified root and no agents',alternatives='Generic worktree/commit/scratch deletion workflow',consequences='No Git recovery; protect original artifacts by hashes'),dict(date='2026-10-02',decision='Keep clinical outcome gate closed; start both feasibility pilots',evidence='Pass 5: negative adequacy and independent review pending',alternatives='Force binary labels or postpone extension until review',consequences='No amended performance claim; extension may proceed independently'),dict(date='2026-10-02',decision='Use standard library and available NumPy for initial work',evidence='Bundled Python NumPy 2.3.5; SciPy import failed',alternatives='Install dependencies before data inspection',consequences='Do not rely on unavailable optimizer; runtime limits recorded')])
    print(json.dumps({'protected_files':len(json.loads(manifest.read_text())), 'pilot_rows':len(selected),'multi_agent_leads':multi,'controls':controls,'metadata_cells':len(meta),'samples':len({r['sample_ID'] for r in meta}),'NBM_samples':len({r['sample_ID'] for r in meta if r['disease_stage']=='NBM'}),'MM_samples':len({r['sample_ID'] for r in meta if r['disease_stage']=='MM'})},indent=2))

if __name__=='__main__': main()
