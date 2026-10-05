"""Read public primary-source material with bounded downloads and provenance."""
import argparse,hashlib,json,pathlib,urllib.request,datetime
p=argparse.ArgumentParser();p.add_argument('url');p.add_argument('name');args=p.parse_args()
root=pathlib.Path(__file__).resolve().parent
if pathlib.Path(args.name).name!=args.name:raise ValueError('Use basename')
request=urllib.request.Request(args.url,headers={'User-Agent':'Mozilla/5.0 (research data inspection)'})
with urllib.request.urlopen(request,timeout=45) as response:
    body=response.read(80_000_001)
    if len(body)>80_000_000:raise ValueError('Download exceeds 80MB inspection limit')
    metadata={'url':args.url,'resolved_url':response.url,'status':response.status,'content_type':response.headers.get('Content-Type'),'bytes':len(body),'sha256':hashlib.sha256(body).hexdigest(),'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(root/args.name).write_bytes(body);(root/(args.name+'.provenance.json')).write_text(json.dumps(metadata,indent=2));print(json.dumps(metadata))
