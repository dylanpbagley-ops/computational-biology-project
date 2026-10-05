"""Download only public morphology tables and retain selected numeric columns."""
import concurrent.futures,csv,datetime,gzip,hashlib,io,json,pathlib,random,time,urllib.request
ROOT=pathlib.Path(__file__).resolve().parent
FEATURES=['AreaShape_Area','AreaShape_Eccentricity','AreaShape_Solidity']+['Intensity_MeanIntensity_'+c for c in ['DNA','AGP','ER','Mito','RNA']]
with (ROOT.parent/'direction_landscape_Oct1/cp_platemap.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
rows=sorted(rows,key=lambda x:x['well_position']);seen=set();treated=[];controls=[]
for r in rows:
    if not r['broad_sample']:controls.append(r)
    elif r['broad_sample'] not in seen:treated.append(r);seen.add(r['broad_sample'])
rng=random.Random(20261002); chosen=rng.sample(treated,6)+rng.sample(controls,6)
(ROOT/'selection.json').write_text(json.dumps(chosen,indent=2))
BASE='https://cellpainting-gallery.s3.amazonaws.com/cpg0000-jump-pilot/source_4/workspace/analysis/2020_11_04_CPJUMP1/'

def get(item):
    plate,well,site=item;name=f'{plate}_{well}_{site}';out=ROOT/(name+'.csv.gz');prov=ROOT/(name+'.json')
    if out.exists() and prov.exists():return json.loads(prov.read_text())
    url=f'{BASE}{plate}/analysis/{plate}-{well}-{site}/Nuclei.csv'
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url,timeout=45) as response:
                assert response.status==200
                body=response.read(25_000_001)
            assert len(body)<25_000_001
            break
        except Exception:
            if attempt==2:raise
            time.sleep(1)
    reader=csv.reader(io.StringIO(body.decode('utf-8-sig')));header=next(reader)
    idx=[header.index(f) for f in ['ImageNumber','ObjectNumber']+FEATURES]
    slim=[];ids=set()
    for row in reader:
        v=[row[i] for i in idx];identity=tuple(v[:2]);assert identity not in ids;ids.add(identity);slim.append(v)
    assert slim
    with gzip.open(out,'wt',newline='') as f:
        w=csv.writer(f);w.writerow(['ImageNumber','ObjectNumber']+FEATURES);w.writerows(slim)
    meta={'url':url,'sha256_source':hashlib.sha256(body).hexdigest(),'source_bytes':len(body),
        'plate':plate,'well':well,'site':site,'objects':len(slim),'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'sha256_extract':hashlib.sha256(out.read_bytes()).hexdigest()}
    prov.write_text(json.dumps(meta,indent=2));return meta

if __name__=='__main__':
    tasks=[(p,r['well_position'],s) for p in ['BR00117011','BR00117012'] for r in chosen for s in range(1,10)]
    results=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
        for future in concurrent.futures.as_completed([pool.submit(get,t) for t in tasks]):
            results.append(future.result())
            if len(results)%18==0:print('Retrieved',len(results),'of',len(tasks),'sites',flush=True)
    (ROOT/'download_manifest.json').write_text(json.dumps(results,indent=2))
    print('Source bytes',sum(x['source_bytes'] for x in results),'objects',sum(x['objects'] for x in results),flush=True)
