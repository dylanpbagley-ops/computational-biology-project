import sys, os, glob, json, numpy as np
SRC = os.path.expanduser('~/mnt/outputs/ev_capture_model/src'); sys.path.insert(0, SRC)
import aim12_single_cell_cytof_reordering as a12
RAW = os.path.expanduser('~/mnt/outputs/ev_capture_data/raw_cytof')
rng = np.random.default_rng(0)

ANTIGENS = ['CD19','CD20','CD22','CD79B','CD38','CD40','CD24','CD21','CD72','CD48','CD84',
            'CD43','CD44','CD27','CD25','CD23','CD80','CD86','CD83','CXCR4','CXCR5','PD-L1',
            'FasR','TRAIL-R1','IL-4R','CD49D','CD10','HLA-DR','HLA-ABC','HLA-E','CCR4','BTLA',
            'CD45RA','IgM','IgD','IgKap','IgLam']

def auc_rank(x, y):
    """P(X > Y) via rank statistic; ties get 0.5 credit."""
    n1, n2 = len(x), len(y)
    z = np.concatenate([x, y])
    r = np.empty(len(z))
    order = np.argsort(z, kind='mergesort')
    zs = z[order]
    # average ranks for ties
    i = 0; rank = np.empty(len(z))
    while i < len(z):
        j = i
        while j+1 < len(z) and zs[j+1] == zs[i]: j += 1
        rank[i:j+1] = (i + j)/2.0 + 1
        i = j+1
    r[order] = rank
    return (r[:n1].sum() - n1*(n1+1)/2.0) / (n1*n2)

rows = {}
for f in sorted(glob.glob(os.path.join(RAW,'*.fcs'))):
    sid = os.path.basename(f).replace('.fcs','')
    pid = sid.split('_')[-1]
    raw = open(f,'rb').read()
    h = a12.parse_fcs3_header(raw); tk, stains = a12.parse_fcs3_text_segment(raw, h)
    tot=int(tk['TOT']); npar=int(tk['PAR'])
    arr = np.frombuffer(raw, dtype='>f4', count=tot*npar, offset=h['data_start']).reshape(tot,npar)
    def col(m):
        return arr[:, a12.find_parameter_index_by_stain_substring(stains, m)].astype(np.float32)
    dna = col('193Ir-DNA'); via = col('195Pt-Viability'); bead = col('140Ce-Beads')
    cd3 = col('CD3'); cd19 = col('CD19'); cd20 = col('CD20'); cd45 = col('CD45')
    # batch-invariant QC gates (percentile based)
    intact = (dna > np.percentile(dna,5)) & (dna < np.percentile(dna,99.5))
    live   = via < np.percentile(via, 95)
    nobead = bead < np.percentile(bead, 99)
    leuko  = cd45 > np.percentile(cd45, 2)
    qc = intact & live & nobead & leuko
    # lineage gates using per-file positivity cutoffs
    c3 = np.percentile(cd3[qc], 50); c19 = np.percentile(cd19[qc], 50); c20 = np.percentile(cd20[qc], 50)
    B = qc & ((cd19 > c19) | (cd20 > c20)) & (cd3 <= c3)
    T = qc & (cd3 > c3) & (cd19 <= c19) & (cd20 <= c20)
    if B.sum() < 500 or T.sum() < 500:
        print('SKIP', sid, B.sum(), T.sum()); continue
    bi = np.flatnonzero(B); ti = np.flatnonzero(T)
    if len(bi) > 20000: bi = rng.choice(bi, 20000, replace=False)
    if len(ti) > 20000: ti = rng.choice(ti, 20000, replace=False)
    rec = {'n_B': int(B.sum()), 'n_T': int(T.sum()), 'batch': sid.split('_')[0]}
    for m in ANTIGENS:
        try: v = col(m)
        except Exception: continue
        rec[m] = {'auc': float(auc_rank(v[bi].astype(float), v[ti].astype(float))),
                  'mean_B': float(v[B].mean()), 'med_B': float(np.median(v[B]))}
    rows[pid] = rec
    print('%-8s B=%6d T=%6d  CD19 auc=%.3f meanB=%7.1f | HLA-DR auc=%.3f meanB=%7.1f | CD38 auc=%.3f' % (
        pid, rec['n_B'], rec['n_T'], rec['CD19']['auc'], rec['CD19']['mean_B'],
        rec['HLA-DR']['auc'], rec['HLA-DR']['mean_B'], rec['CD38']['auc']), flush=True)

json.dump(rows, open(os.path.expanduser('~/scratch/pilot_roc.json'),'w'))
print('DONE', len(rows))
