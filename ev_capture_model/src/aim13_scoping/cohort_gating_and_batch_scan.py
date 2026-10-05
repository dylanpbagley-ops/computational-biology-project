import sys, os, glob, json, numpy as np
SRC = os.path.expanduser('~/mnt/outputs/ev_capture_model/src')
sys.path.insert(0, SRC)
import aim12_single_cell_cytof_reordering as a12
RAW = os.path.expanduser('~/mnt/outputs/ev_capture_data/raw_cytof')
out = {}
for f in sorted(glob.glob(os.path.join(RAW,'*.fcs'))):
    sid = os.path.basename(f).replace('.fcs','')
    raw = open(f,'rb').read()
    h = a12.parse_fcs3_header(raw); tk, stains = a12.parse_fcs3_text_segment(raw, h)
    tot=int(tk['TOT']); npar=int(tk['PAR'])
    arr = np.frombuffer(raw, dtype='>f4', count=tot*npar, offset=h['data_start']).reshape(tot,npar).astype(np.float32)
    def col(m):
        i = a12.find_parameter_index_by_stain_substring(stains, m)
        return arr[:, i].astype(float)
    dna = col('193Ir-DNA'); cd3=col('CD3'); cd19=col('CD19'); cd20=col('CD20')
    r = {'n':tot,'batch':sid.split('_')[0],'pid':sid.split('_')[-1],
         'dna_med':float(np.median(dna)),
         'cd19_mean_all':float(cd19.mean()),'cd19_med_all':float(np.median(cd19)),
         'cd19_zero_all':float((cd19<=0).mean())}
    # B gate: CD19 or CD20 clearly positive AND CD3 low
    t_hi = np.percentile(cd3, 75)
    B = ((cd19 > 10) | (cd20 > 10)) & (cd3 < 10)
    T = cd3 > 10
    r['bfrac']=float(B.mean()); r['tfrac']=float(T.mean())
    if B.sum()>200:
        r['cd19_mean_B']=float(cd19[B].mean()); r['cd19_med_B']=float(np.median(cd19[B]))
        r['cd19_zero_B']=float((cd19[B]<=0).mean())
        r['dna_med_B']=float(np.median(dna[B]))
    out[sid]=r
    print('%-22s n=%7d dna=%7.1f  cd19_all=%8.2f  cd19_B=%8.2f  B%%=%.2f T%%=%.2f' % (
        sid,tot,r['dna_med'],r['cd19_mean_all'],r.get('cd19_mean_B',-1),r['bfrac'],r['tfrac']), flush=True)
json.dump(out, open(os.path.expanduser('~/scratch/cohort_scan.json'),'w'), indent=1)
print('DONE')
