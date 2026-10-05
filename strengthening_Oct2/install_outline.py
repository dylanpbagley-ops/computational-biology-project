import pathlib,json,hashlib
ROOT=pathlib.Path(__file__).resolve().parent;BASE=ROOT.parent
changes=json.loads((ROOT/'document_changes.json').read_text());manifest=json.loads((ROOT/'before_hashes.json').read_text())
pairs=[('ev_capture_model/docs/ISEF_Paper_Skeleton.md',(ROOT/'current_outline.md').read_text(encoding='utf-8')),('Poster, Paper, and Write Up/ISEF_Paper_Skeleton.md','# Superseded local copy\n\nThe current advisory outline is in ev_capture_model/docs/ISEF_Paper_Skeleton.md in the outputs root. The previous version of this file is preserved under strengthening_Oct2/before/. Start with PROJECT_STATUS_Oct2.md.\n')]
for name,content in pairs:
    p=BASE/name;before=hashlib.sha256(p.read_bytes()).hexdigest();assert before==manifest[str(pathlib.Path(name))]
    p.write_text(content,encoding='utf-8');changes.append(dict(file=name,before=before,after=hashlib.sha256(p.read_bytes()).hexdigest()))
(ROOT/'document_changes.json').write_text(json.dumps(changes,indent=2));print('Installed canonical advisory outline and pointer.')
