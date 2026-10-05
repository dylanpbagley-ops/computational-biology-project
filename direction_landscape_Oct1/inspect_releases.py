import json,pathlib,zipfile,csv
from html.parser import HTMLParser
root=pathlib.Path(__file__).resolve().parent
class Parser(HTMLParser):
    def __init__(self):super().__init__();self.text=[];self.links=[];self.skip=0
    def handle_starttag(self,tag,attrs):
        if tag in ['script','style']:self.skip+=1
        if tag in ['p','h1','h2','h3','h4','li','section']:self.text.append('\n')
        if tag=='a':
            href=dict(attrs).get('href','')
            if any(t in href.lower() for t in ['mediaobjects','file_stream','github','zenodo','source-data']):self.links.append(href)
    def handle_endtag(self,tag):
        if tag in ['script','style']:self.skip=max(0,self.skip-1)
    def handle_data(self,value):
        if not self.skip:self.text.append(value)
for p in root.glob('*.html'):
    parser=Parser();parser.feed(p.read_text(encoding='utf-8',errors='replace'));p.with_suffix('.txt').write_text(''.join(parser.text),encoding='utf-8')
    print(p.name,json.dumps(sorted(set(parser.links))))
inventory={}
for p in root.glob('*.zip'):
    destination=(root/p.stem).resolve();entries=[]
    with zipfile.ZipFile(p) as z:
        for i in z.infolist():
            if i.is_dir() or i.filename.startswith('__MACOSX/'):continue
            path=(destination/i.filename).resolve()
            if not path.is_relative_to(destination):raise ValueError('Unsafe ZIP path')
            if i.file_size>80_000_000:raise ValueError('Large member')
            entries.append({'name':i.filename,'bytes':i.file_size});path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(z.read(i))
    inventory[p.name]=entries;print(p.name,len(entries),json.dumps(entries[:14]))
(root/'archive_inventory.json').write_text(json.dumps(inventory,indent=2))
