"""Descriptive source-table reproduction, not predictive validation."""
import csv, hashlib, json, pathlib, statistics
root=pathlib.Path(__file__).resolve().parent
records=[]
for filename,measurement in [('6d.txt','specific_cytotoxicity_percent'),('6e.txt','T_cell_count'),('6i.txt','residual_tumor_antigen_gMFI')]:
    path=root/'barisa_source'/filename
    rows=list(csv.reader(path.open(encoding='utf-8-sig'),delimiter='\t'))
    groups={}
    for col,(construct,rep) in enumerate(zip(rows[0],rows[1])):
        if col and construct and rep:
            groups.setdefault((construct,rep.split('_')[0]),[]).append(col)
    for row in rows[2:]:
        if not row or not row[0]: continue
        for (construct,donor),cols in groups.items():
            vals=[float(row[c]) for c in cols if c<len(row) and row[c].strip()]
            records.append({'source':filename,'measurement':measurement,'day':float(row[0]),'construct':construct,'source_donor_label':donor,'wells':len(vals),'mean':statistics.mean(vals)})
with (root/'descriptive_donor_means.csv').open('w',newline='',encoding='utf-8') as f:
    writer=csv.DictWriter(f,fieldnames=records[0].keys());writer.writeheader();writer.writerows(records)
summary=[]
for measurement in sorted(set(r['measurement'] for r in records)):
    for construct in sorted(set(r['construct'] for r in records if r['measurement']==measurement)):
        for day in sorted(set(r['day'] for r in records if r['measurement']==measurement and r['construct']==construct)):
            vals=[r['mean'] for r in records if r['measurement']==measurement and r['construct']==construct and r['day']==day]
            summary.append({'measurement':measurement,'construct':construct,'day':day,'source_donor_labels':len(vals),'equal_label_mean':statistics.mean(vals),'label_means':vals})
result={'status':'Descriptive reproduction only. Three source donor labels; identity and replication discrepancy with paper require clarification. No hypothesis test or forecast score. No clipping normalized cytotoxicity.','summary':summary,'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [root/'barisa_source'/n for n in ['6d.txt','6e.txt','6i.txt']]}}
(root/'reproduction_results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
