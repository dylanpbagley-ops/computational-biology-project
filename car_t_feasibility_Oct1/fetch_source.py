"""Retrieve primary source material and record checksums; no remote code execution."""
import argparse, hashlib, json, pathlib, urllib.request, datetime

p = argparse.ArgumentParser()
p.add_argument('url')
p.add_argument('name')
p.add_argument('--range', dest='byte_range')
args = p.parse_args()
root = pathlib.Path(__file__).resolve().parent
name = pathlib.Path(args.name)
if name.name != args.name:
    raise ValueError('Use a single filename inside the investigation directory')
headers = {'User-Agent': 'Mozilla/5.0 scientific-data-feasibility/1.0'}
if args.byte_range:
    headers['Range'] = 'bytes=' + args.byte_range
request = urllib.request.Request(args.url, headers=headers)
with urllib.request.urlopen(request, timeout=45) as response:
    body = response.read(75_000_001)
    if len(body) > 75_000_000:
        raise ValueError('Download exceeded the 75 MB limit')
    metadata = {'requested_url': args.url, 'resolved_url': response.url,
                'status': response.status, 'headers': dict(response.headers),
                'bytes': len(body), 'sha256': hashlib.sha256(body).hexdigest(),
                'retrieved_utc': datetime.datetime.now(datetime.timezone.utc).isoformat()}
(root / args.name).write_bytes(body)
(root / (args.name + '.provenance.json')).write_text(json.dumps(metadata, indent=2), encoding='utf-8')
print(json.dumps({k: metadata[k] for k in ('requested_url','status','bytes','sha256')}))
