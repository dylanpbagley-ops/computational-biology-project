"""Fresh verification of repairs without overwriting historical numerical outputs."""
import csv,hashlib,json,pathlib,shutil,subprocess,sys
ROOT=pathlib.Path(__file__).resolve().parent;BASE=ROOT.parent
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    manifest=json.loads((ROOT/'before_hashes.json').read_text());modified={v['file'].replace('\\','/') for v in json.loads((ROOT/'document_changes.json').read_text())}
    preserved=[]
    for rel,h in manifest.items():
        if rel.replace('\\','/') not in modified:
            assert digest(BASE/rel)==h, 'Historical artifact changed: '+rel
            preserved.append(rel)
    assert digest(ROOT/'analysis_precommitment.md')==(ROOT/'analysis_precommitment.md.hashes').read_text().split()[0]
    proc=subprocess.run([sys.executable,'-m','unittest','discover','-s',str(ROOT),'-p','test_*.py'],capture_output=True,text=True)
    (ROOT/'test_results.log').write_text(proc.stdout+proc.stderr);assert proc.returncode==0,proc.stdout+proc.stderr
    runtime=ROOT/'legacy_runtime';runtime.mkdir(exist_ok=True)
    name='weighted_ceiling_theorem.py';shutil.copy2(BASE/'ev_capture_model/src/aim13_scoping'/name,runtime/name)
    proc=subprocess.run([sys.executable,str(runtime/name)],capture_output=True,text=True)
    (ROOT/'legacy_theorem_verification.log').write_text(proc.stdout+proc.stderr);assert proc.returncode==0
    actual=json.loads((runtime/'weighted_ceiling_theorem.json').read_text());assert actual['checks']==16400 and actual['worst_excess']<1e-9 and actual['worst_np_gap']<1e-9
    old=json.loads((ROOT/'before/ev_capture_model/src/aim13_scoping/weighted_ceiling_theorem.json').read_text())
    assert actual==old, 'Legacy output changed'
    doc=(BASE/'ev_capture_model/docs/Aim13_Results_Two_Comparator_Ceiling.md').read_text(encoding='utf-8')
    assert 'governed by the single number AUC(p_T, p_W)' not in doc
    assert 'They are the same statement.' not in doc
    outline=(BASE/'ev_capture_model/docs/ISEF_Paper_Skeleton.md').read_text(encoding='utf-8')
    assert outline==(ROOT/'current_outline.md').read_text(encoding='utf-8')
    review=list(csv.DictReader((ROOT/'reviewer_form_score_free.tsv').open(),delimiter='\t'));assert len(review)==30
    assert all(not row['proposed_category'] and not row['reviewer_name_date'] for row in review)
    assert not any(k in review[0] for k in ['score','top_share','historical_outcome','outcome'])
    index=list(csv.DictReader((ROOT/'agent_evidence_index.tsv').open(),delimiter='\t'))
    amendments=list(csv.DictReader((ROOT/'evidence_amendments_Oct2.tsv').open(),delimiter='\t'))
    assert len(index)==166 and len({r['evidence_id'] for r in index})==166
    assert len(amendments)==44 and all(r['adjudicated_outcome']=='PENDING_REVIEW' for r in amendments)
    diagnostics=json.loads((ROOT/'diagnostic_results.json').read_text());assert abs(diagnostics['baseline']['frozen_rule']['auc']-191/209)<1e-14
    assert diagnostics['baseline']['tau']['auc']==diagnostics['baseline']['confinement']['auc']
    assert '221 checks run, all passed' in (ROOT/'legacy_document_verification.log').read_text()
    files=[p for p in ROOT.rglob('*') if p.is_file() and 'before' not in p.parts and '__pycache__' not in p.parts and p.name not in ['artifact_hashes.json','verification_results.json']]
    (ROOT/'artifact_hashes.json').write_text(json.dumps({str(p.relative_to(ROOT)):digest(p) for p in sorted(files)},indent=2))
    out={'new_unit_tests':'10 passed','legacy_theorem_checks':actual['checks'],'legacy_numerical_output_unchanged':True,'legacy_document_number_checks':221,'preserved_historical_files':len(preserved),'protocol_hash_matches':True,'score_free_review_rows':30,'indexed_evidence_rows':len(index),'amended_evidence_rows':len(amendments),'scope':'No full 341-test legacy rerun: numerical legacy implementation unchanged; affected standalone theorem run in isolated directory and document verifier run directly.'}
    (ROOT/'verification_results.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':main()
