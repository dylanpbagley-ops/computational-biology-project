import io, json, pathlib, zipfile
root=pathlib.Path(__file__).resolve().parent
start=2437400000
body=(root/'aml_data_zip_tail.bin').read_bytes()

class Tail(io.RawIOBase):
    def __init__(self): self.position=start
    def seekable(self): return True
    def seek(self,offset,whence=0):
        self.position=offset if whence==0 else self.position+offset if whence==1 else start+len(body)+offset
        return self.position
    def tell(self): return self.position
    def read(self,size=-1):
        if self.position<start: raise ValueError('Requested bytes outside fetched archive tail')
        begin=self.position-start
        out=body[begin:] if size<0 else body[begin:begin+size]
        self.position+=len(out)
        return out
with zipfile.ZipFile(Tail()) as archive:
    inventory=[{'name':i.filename,'bytes':i.file_size,'compressed_bytes':i.compress_size,'header_offset':i.header_offset,'compression':i.compress_type,'crc':i.CRC} for i in archive.infolist()]
    for member in ('LT6.csv','2023-09-19_kept_chains.jld2'):
        (root/member).write_bytes(archive.read(member))
(root/'aml_data_remote_inventory.json').write_text(json.dumps(inventory,indent=2),encoding='utf-8')
print(json.dumps({'entries':len(inventory),'total_uncompressed':sum(i['bytes'] for i in inventory),'small_data_files':[i for i in inventory if i['name'].lower().endswith(('.csv','.data','.txt','.xlsx','.fcs'))]},indent=2))
