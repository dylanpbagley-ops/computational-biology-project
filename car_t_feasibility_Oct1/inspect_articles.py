import csv, hashlib, json, pathlib, re
from html.parser import HTMLParser
root=pathlib.Path(__file__).resolve().parent
class Text(HTMLParser):
    def __init__(self): super().__init__();self.parts=[];self.skip=0;self.links=[]
    def handle_starttag(self,tag,attrs):
        if tag in ('script','style'): self.skip+=1
        if tag in ('p','h1','h2','h3','h4','li','section'): self.parts.append('\n')
        if tag=='a':
            href=dict(attrs).get('href','')
            if any(t in href.lower() for t in ('supp','mmc','pone','zip','xlsx','csv','nihms','github')): self.links.append(href)
    def handle_endtag(self,tag):
        if tag in ('script','style'): self.skip=max(0,self.skip-1)
    def handle_data(self,data):
        if not self.skip:self.parts.append(data)
for name in ['majzner','aml']:
    t=Text();t.feed((root/f'{name}_article.html').read_text(encoding='utf-8',errors='replace'))
    (root/f'{name}_article.txt').write_text(''.join(t.parts),encoding='utf-8')
    print(name,'supplement links',json.dumps(sorted(set(t.links))))
data=list(csv.DictReader((root/'LT6.csv').open()))
summary={'rows':len(data),'columns':list(data[0]),'unique_values':{col:sorted(set(r[col] for r in data)) for col in ['effector','target','ratio','time']},'independent_replicate_column':False,'antigen_density_column':False,'sha256':hashlib.sha256((root/'LT6.csv').read_bytes()).hexdigest(),'origin':'Member LT6.csv in Zenodo 12823159 data.zip; CRC checked by Python zipfile; extracted from HTTP byte range documented in aml_data_zip_tail.bin.provenance.json. Full archive checksum not verified because full 2.4GB archive was not downloaded.'}
(root/'LT6_inventory.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary))
