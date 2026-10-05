from pathlib import Path, PurePosixPath
import csv, hashlib, json
STAGE=Path(__file__).resolve().parent
ROOT=STAGE.parents[1]
MAX_FILE=25*1024*1024
MAX_TOTAL=75*1024*1024
REQUIRED=['decision.md','STATUS.md','candidate_gates.tsv','data_cards.tsv','data_manifest.json','output_hashes.json','preservation_before.json','search_log.tsv','sources.tsv','validate_artifacts.py','prior_art/STATUS.md','prior_art/search_log.tsv','prior_art/sources.tsv','prior_art/review.md','flag_review/STATUS.md','flag_review/review.md','flag_review/search_log.tsv','flag_review/sources.tsv','pilot/protocol.json','pilot/hash_definition_note.md','pilot/preservation_before.json','pilot/inspect_transport.py','pilot/transport_helpers.py','pilot/test_transport.py','pilot/gene_effects.tsv','pilot/guide_eligibility.tsv','pilot/transport_flags.tsv','pilot/summary.json','pilot/decision.md']
def fail(message): raise SystemExit('FAIL: '+message)
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
 return h.hexdigest()
def reject(value): raise ValueError('non-finite JSON constant '+value)
def parse_json(path):
 try: return json.loads(path.read_text(encoding='utf-8'),parse_constant=reject)
 except Exception as e: fail(f'invalid JSON {path.relative_to(STAGE)}: {e}')
def safe_stage_path(rel):
 if not isinstance(rel,str) or not rel or '\\' in rel: fail(f'unsafe relative path: {rel!r}')
 pp=PurePosixPath(rel)
 if pp.is_absolute() or any(x in ('.','..') for x in pp.parts): fail(f'unsafe relative path: {rel!r}')
 p=(STAGE/Path(*pp.parts)).resolve()
 try: p.relative_to(STAGE)
 except ValueError: fail(f'path escapes stage: {rel!r}')
 if not p.is_file(): fail(f'missing stage file: {rel}')
 return p
def verify_tsv(rel,ncols):
 p=safe_stage_path(rel)
 try:
  with p.open('r',encoding='utf-8',newline='') as f:
   rows=csv.reader(f,delimiter=chr(9)); header=next(rows)
   if len(header)!=ncols: fail(f'{rel}: header width {len(header)} != {ncols}')
   for i,row in enumerate(rows,2):
    if len(row)!=ncols: fail(f'{rel}:{i}: row width {len(row)} != {ncols}')
 except StopIteration: fail(f'{rel}: empty TSV')
def main():
 for rel in REQUIRED: safe_stage_path(rel)
 for p in STAGE.rglob('*.json'): parse_json(p)
 preservation=parse_json(safe_stage_path('preservation_before.json'))
 for item in preservation.get('protected_inputs',[]):
  p=(ROOT/item['path']).resolve()
  if not p.is_file(): fail('missing protected input '+item['path'])
  if p.stat().st_size!=item['bytes'] or sha(p)!=item['sha256']: fail('protected input changed: '+item['path'])
 if not preservation.get('authorized_mutable'): fail('authorized mutable baseline was not recorded')
 for rel,width in [('search_log.tsv',7),('sources.tsv',9),('data_cards.tsv',11),('candidate_gates.tsv',13),('prior_art/search_log.tsv',6),('prior_art/sources.tsv',10),('flag_review/search_log.tsv',6),('flag_review/sources.tsv',10),('pilot/gene_effects.tsv',26),('pilot/guide_eligibility.tsv',12),('pilot/transport_flags.tsv',25)]: verify_tsv(rel,width)
 manifest=parse_json(safe_stage_path('data_manifest.json'))
 files=manifest.get('files'); resources=manifest.get('resources')
 if not isinstance(files,list) or not files or not isinstance(resources,list): fail('manifest resources/files missing')
 if manifest.get('resource_count')!=len(resources) or len(resources)>4: fail('resource count mismatch/budget exceeded')
 listed,total,download_count=set(),0,0
 for item in files:
  rel=item.get('path')
  if rel in listed: fail('duplicate manifest path '+str(rel))
  listed.add(rel); p=safe_stage_path(rel); size=p.stat().st_size
  if size!=item.get('bytes') or sha(p)!=item.get('sha256'): fail('wrong size/hash for '+rel)
  if item.get('downloaded'):
   download_count+=1; total+=size
   if size>MAX_FILE: fail('download exceeds per-file cap: '+rel)
   if not item.get('official_url'): fail('download lacks official URL: '+rel)
 actual_data={p.relative_to(STAGE).as_posix() for p in (STAGE/'data').rglob('*') if p.is_file()}
 if actual_data!=listed: fail(f'unlisted/missing data files: actual={sorted(actual_data)} manifest={sorted(listed)}')
 if total>MAX_TOTAL or total!=manifest.get('downloaded_bytes_total'): fail('aggregate download cap/bytes mismatch')
 if download_count!=manifest.get('download_count'): fail('download count mismatch')
 pp_path=safe_stage_path('pilot/protocol.json')
 pp=parse_json(pp_path)
 protocol_hash=pp.get('protocol_sha256')
 protocol_file_hash=sha(pp_path)
 if protocol_file_hash!='9cf1278621e7b520bfd2bbfbac0a35065bbe60fb4f6b7ddd1f467e75add6b268': fail('frozen pilot protocol full-file SHA mismatch')
 if protocol_hash!='45103a84c6d99f9bea919f8ae5b3511f0fbf515b1e52b5d4f976cdea6b2167a5': fail('frozen pilot protocol stored hash changed')
 single=dict(pp); single.pop('protocol_sha256',None)
 single_hash=hashlib.sha256(json.dumps(single,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode('utf-8')).hexdigest()
 canonical=dict(pp); canonical.pop('protocol_sha256',None); canonical.pop('protocol_hash_definition',None)
 calculated=hashlib.sha256(json.dumps(canonical,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode('utf-8')).hexdigest()
 if single_hash!='c6a95fe8cbd20698872932361ce33eea112922805a69b2f6e956cefce124ed02': fail('protocol single-exclusion hash changed')
 if protocol_hash!=calculated: fail('stored protocol hash does not match both-metadata-field exclusion')
 note=safe_stage_path('pilot/hash_definition_note.md').read_text(encoding='utf-8')
 for token in (protocol_file_hash,protocol_hash,single_hash,'wording is inconsistent'):
  if token not in note: fail('hash-definition erratum does not acknowledge frozen inconsistency')
 ps=parse_json(safe_stage_path('pilot/preservation_before.json'))
 if ps.get('protocol',{}).get('sha256')!=protocol_hash: fail('pilot preservation protocol hash mismatch')
 summary=parse_json(safe_stage_path('pilot/summary.json'))
 if summary.get('protocol_sha256')!=protocol_hash: fail('summary protocol hash mismatch')
 if summary.get('reviewed_plan_sha256')!=ps.get('reviewed_plan',{}).get('sha256'): fail('reviewed plan hash mismatch')
 expected={x['path']:x['sha256'] for x in files if x.get('downloaded')}
 if summary.get('scientific_reproduction_checks',{}).get('raw_input_hashes')!=expected: fail('summary raw-input hash set mismatch')
 frozen={x['path']:x['sha256'] for x in ps.get('frozen_prior_art',[])}
 if summary.get('scientific_reproduction_checks',{}).get('prior_art_hashes')!=frozen: fail('summary prior-art hash set mismatch')
 hashes=parse_json(safe_stage_path('output_hashes.json')); entries=hashes.get('files')
 if not isinstance(entries,list): fail('output hash inventory has no files array')
 inventory=set()
 for item in entries:
  rel=item.get('path')
  if rel in inventory: fail('duplicate output hash path '+str(rel))
  inventory.add(rel); p=safe_stage_path(rel)
  if p.stat().st_size!=item.get('bytes') or sha(p)!=item.get('sha256'): fail('output hash mismatch: '+rel)
 actual_stage={p.relative_to(STAGE).as_posix() for p in STAGE.rglob('*') if p.is_file() and p.name!='output_hashes.json'}
 if inventory!=actual_stage: fail(f'output inventory mismatch: actual={sorted(actual_stage)} listed={sorted(inventory)}')
 print(f"PASS: required files, finite JSON, 11 TSVs, protocol hash-definition erratum/input/prior-art hashes, {len(resources)}-resource/{download_count}-download manifest ({total} bytes), and {len(inventory)} hashed stage outputs")
if __name__=='__main__': main()
