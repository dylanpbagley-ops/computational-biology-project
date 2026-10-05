import pathlib,zipfile,xml.etree.ElementTree as ET,json,csv
root=pathlib.Path(__file__).resolve().parent
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
with zipfile.ZipFile(root/'hiv_combinations_supplement.docx') as z:
    doc=ET.fromstring(z.read('word/document.xml'))
tables=[]
for i,table in enumerate(doc.findall('.//w:tbl',ns)):
    rows=[]
    for row in table.findall('w:tr',ns):
        rows.append([''.join(c.itertext()) if False else ''.join(c.findall('.//w:t',ns)[j].text or '' for j in range(len(c.findall('.//w:t',ns)))) for c in row.findall('w:tc',ns)])
    tables.append({'table':i,'rows':len(rows),'max_columns':max(map(len,rows),default=0),'first_rows':rows[:5]})
    with (root/f'hiv_supp_table_{i}.csv').open('w',newline='',encoding='utf-8') as f:csv.writer(f).writerows(rows)
(root/'hiv_supp_table_inventory.json').write_text(json.dumps(tables,indent=2),encoding='utf-8')
text='\n'.join(''.join(x.text or '' for x in p.findall('.//w:t',ns)) for p in doc.findall('.//w:p',ns))
(root/'hiv_supplement_text.txt').write_text(text,encoding='utf-8')
print(json.dumps(tables,indent=2)[:18000])
