"""
Re-run Aim 12's patient-reordering finding on properly gated tumor cells.

Builds, per patient, three CD19 populations:
  ungated  -- every event in the file (what Aim 12 actually used)
  Bgated   -- QC-gated CD19/CD20+ CD3- B cells
  clone    -- light-chain-restricted malignant cells only (where resolvable)

Then reports the reordering metric three ways, plus the exact continuum version made
possible by the identity escape_fraction == F(n*).
Resumable: run repeatedly until it prints DONE.
"""
import os, sys, glob, json, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE,'..')))
import aim12_single_cell_cytof_reordering as a12
DATA=os.path.abspath(os.path.join(HERE,'..','..','..','ev_capture_data'))
OUT=os.path.join(DATA,'cd19_populations'); os.makedirs(OUT, exist_ok=True)
lc=json.load(open(os.path.join(HERE,'light_chain_calls.json'))); calls=lc['calls']

def build(path):
    raw=open(path,'rb').read(); h=a12.parse_fcs3_header(raw); tk,st=a12.parse_fcs3_text_segment(raw,h)
    tot,npar=int(tk['TOT']),int(tk['PAR'])
    arr=np.frombuffer(raw,dtype='>f4',count=tot*npar,offset=h['data_start']).reshape(tot,npar)
    g=lambda m: arr[:,a12.find_parameter_index_by_stain_substring(st,m)].astype(np.float32)
    dna,via,bd,cd45=g('193Ir-DNA'),g('195Pt-Viability'),g('140Ce-Beads'),g('CD45')
    cd3,cd19,cd20=g('CD3'),g('CD19'),g('CD20'); ik,il=g('IgKap'),g('IgLam')
    qc=((dna>np.percentile(dna,5))&(dna<np.percentile(dna,99.5))&(via<np.percentile(via,95))
        &(bd<np.percentile(bd,99))&(cd45>np.percentile(cd45,2)))
    c3,c19,c20=np.percentile(cd3[qc],50),np.percentile(cd19[qc],50),np.percentile(cd20[qc],50)
    B=qc&((cd19>c19)|(cd20>c20))&(cd3<=c3); T=qc&(cd3>c3)&(cd19<=c19)&(cd20<=c20)
    out={'ungated':cd19.astype(np.float32),'Bgated':cd19[B].astype(np.float32)}
    key=os.path.basename(path).replace('.fcs','')
    call=calls.get(key,{}).get('call','')
    if call in ('kappa-restricted','lambda-restricted') and T.sum()>500:
        kc,lcut=np.percentile(ik[T],99.0),np.percentile(il[T],99.0)
        kp=(ik>kc)&(il<=lcut)&B; lp=(il>lcut)&(ik<=kc)&B
        m = kp if call=='kappa-restricted' else lp
        if m.sum()>=500: out['clone']=cd19[m].astype(np.float32)
    return out

todo=[f for f in sorted(glob.glob(os.path.join(DATA,'raw_cytof','*.fcs')))
      if not os.path.exists(os.path.join(OUT, os.path.basename(f).replace('.fcs','.npz')))]
print('remaining %d' % len(todo))
lim=int(sys.argv[1]) if len(sys.argv)>1 else 14
for f in todo[:lim]:
    o=build(f); np.savez_compressed(os.path.join(OUT, os.path.basename(f).replace('.fcs','.npz')), **o)
    print('  %-22s ungated=%7d Bgated=%7d clone=%s' % (os.path.basename(f)[:-4], len(o['ungated']), len(o['Bgated']),
          len(o['clone']) if 'clone' in o else '-'), flush=True)
print('DONE' if len(todo)<=lim else 'MORE: %d' % (len(todo)-lim))
