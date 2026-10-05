"""
Step 9 (replication in a SECOND disease): load per-cell surface-protein (ADT) counts
from the GSE234261 multiple-myeloma BCMA CAR-T CITE-seq dataset.

Why this dataset: it measures 63 surface proteins on every cell -- the same physical
quantity the DLBCL CyTOF data measures -- in a different cancer (multiple myeloma),
from a different lab, on a different technology (CITE-seq, not mass cytometry), and it
includes the real approved myeloma targets: BCMA (CD269), CD38, SLAMF7 (CD319),
FcRH5 (CD307) and CD138. If the ceiling framework is measuring biology rather than an
artifact of one dataset, it should behave consistently here.

Only the 63 antibody rows are read from each Matrix Market file (streamed, so the
36,602 gene rows are never materialised). Resumable: rerun until it prints DONE.
"""
import os, sys, gzip, glob, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
SRC=os.path.abspath(os.path.join(HERE,'..','..','..','GSE234261_BCMA_CAR-T_data'))
OUT=os.path.abspath(os.path.join(HERE,'..','..','..','ev_capture_data','myeloma_adt'))
os.makedirs(OUT, exist_ok=True)

def load_sample(sid):
    feat=os.path.join(SRC,f'GSE234261_MXMERZ002A_{sid}_features.tsv.gz')
    mtx =os.path.join(SRC,f'GSE234261_MXMERZ002A_{sid}_matrix.mtx.gz')
    names=[]; adt_rows={}
    with gzip.open(feat,'rt') as fh:
        for i,line in enumerate(fh,1):
            p=line.rstrip('\n').split('\t')
            if len(p)>=3 and p[2]=='Antibody Capture':
                adt_rows[i]=len(names); names.append(p[1].replace('_proteona',''))
    with gzip.open(mtx,'rt') as fh:
        for line in fh:
            if line.startswith('%'): continue
            nrow,ncol,nnz=(int(x) for x in line.split()); break
        M=np.zeros((ncol,len(names)),dtype=np.float32)
        for line in fh:
            r,c,v=line.split()
            r=int(r)
            j=adt_rows.get(r)
            if j is not None: M[int(c)-1,j]=float(v)
    return names, M

sids=sorted({os.path.basename(p).split('_')[2] for p in glob.glob(os.path.join(SRC,'*_matrix.mtx.gz'))})
todo=[s for s in sids if not os.path.exists(os.path.join(OUT,f'{s}.npz'))]
print('samples %d, remaining %d' % (len(sids), len(todo)))
lim=int(sys.argv[1]) if len(sys.argv)>1 else 4
for sid in todo[:lim]:
    names,M=load_sample(sid)
    np.savez_compressed(os.path.join(OUT,f'{sid}.npz'), adt=M, names=np.array(names))
    nz=(M.sum(1)>0).mean()
    print('  sample %s: %d cells x %d proteins, %.1f%% cells with any ADT signal, median total ADT/cell %.0f'
          % (sid, M.shape[0], M.shape[1], 100*nz, np.median(M.sum(1))), flush=True)
print('DONE' if len(todo)<=lim else 'MORE: %d' % (len(todo)-lim))
