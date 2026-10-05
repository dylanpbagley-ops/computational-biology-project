import hashlib, json, pathlib, zipfile
from html.parser import HTMLParser

root = pathlib.Path(__file__).resolve().parent
inventories = {}
for basename, subfolder in [('barisa_source_data.zip','barisa_source'),('aml_code.zip','aml_code')]:
    destination = (root / subfolder).resolve()
    entries = []
    with zipfile.ZipFile(root / basename) as archive:
        for member in archive.infolist():
            if member.is_dir() or member.filename.startswith('__MACOSX/'):
                continue
            body = archive.read(member)
            entries.append({'name':member.filename, 'bytes':len(body), 'sha256':hashlib.sha256(body).hexdigest()})
            path = (destination / member.filename).resolve()
            if not path.is_relative_to(destination):
                raise ValueError('Unsafe archive path')
            path.parent.mkdir(parents=True,exist_ok=True)
            path.write_bytes(body)
    inventories[basename] = entries
(root/'file_inventory.json').write_text(json.dumps(inventories,indent=2),encoding='utf-8')
print(json.dumps({name:len(entries) for name,entries in inventories.items()}))
for name in ['3a.txt','3f.txt','3g.txt','4c.txt','6a.txt','6d.txt','6e.txt','6f.txt','6g_B.txt','6i.txt','6n.txt']:
    path = root/'barisa_source'/name
    if path.exists():
        print('\nFILE:',name, '\n', path.read_text(encoding='utf-8-sig',errors='replace')[:2200])

class Text(HTMLParser):
    def __init__(self):
        super().__init__(); self.parts=[]; self.skip=0
    def handle_starttag(self,tag,attrs):
        if tag in ('script','style'): self.skip+=1
        if tag in ('p','h1','h2','h3','h4','li','section'): self.parts.append('\n')
    def handle_endtag(self,tag):
        if tag in ('script','style'): self.skip=max(0,self.skip-1)
    def handle_data(self,data):
        if not self.skip: self.parts.append(data)
t=Text();t.feed((root/'barisa_article.html').read_text(encoding='utf-8',errors='replace'))
(root/'barisa_article.txt').write_text(''.join(t.parts),encoding='utf-8')
