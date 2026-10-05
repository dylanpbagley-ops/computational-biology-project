import csv, hashlib, json, pathlib, zipfile
root=pathlib.Path(__file__).resolve().parent
p=root/'majzner_supplement.docx'
if p.exists() and not zipfile.is_zipfile(p):
    p.rename(root/'majzner_supplement_download_response.html')
    (root/'majzner_supplement_access_note.json').write_text(json.dumps({'status':'Supplement not obtained','reason':'HTTP 200 returned 1817 bytes of HTML, not the advertised 27.9MB DOCX. Do not treat response as supplementary data.','provenance':'majzner_supplement.docx.provenance.json'},indent=2))
tables={}
for p in (root/'barisa_source').glob('*.txt'):
    rows=list(csv.reader(p.open(encoding='utf-8-sig',errors='replace'),delimiter='\t'))
    tables[p.name]={'rows':len(rows),'max_columns':max(map(len,rows),default=0)}
(root/'barisa_table_inventory.json').write_text(json.dumps(tables,indent=2),encoding='utf-8')
hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in root.rglob('*') if p.is_file() and p.name!='evidence_hashes.json'}
(root/'evidence_hashes.json').write_text(json.dumps(hashes,indent=2),encoding='utf-8')
latest=json.loads((root/'aml_latest_metadata.json').read_text())
print(json.dumps({'B7H3_source_tables':len(tables),'FCS_files_in_downloaded_archives':sum(p.suffix.lower()=='.fcs' for d in ['barisa_source','aml_code'] for p in (root/d).rglob('*') if p.is_file()),'AML_latest_record':latest['id'],'AML_latest_modified':latest['modified'],'LT15_training_file_found':any(p.name.lower()=='lt15.data' for p in (root/'aml_code').rglob('*')),'files_hashed':len(hashes)},indent=2))
