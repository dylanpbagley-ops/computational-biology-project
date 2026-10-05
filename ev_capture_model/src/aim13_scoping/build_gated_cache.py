"""
Build a per-file cache of GATED single-cell antigen values for the whole cohort
(41 DLBCL tumor files + 36 reactive-lymph-node control files).

Why a cache: each FCS file is 15-130 MB and reparsing them for every analysis is
slow. This extracts, once, the QC-gated B-cell and T-cell populations for the
antigen panel and stores them as compressed npz.

Gating is deliberately percentile-based (not absolute-intensity based) because
per-batch signal scale varies ~10x across this dataset; percentile gates are
invariant to that scaling. See docs/Aim13_Therapeutic_Index_Ceiling_PROPOSAL.md.

Resumable: already-cached files are skipped. Run repeatedly until it prints DONE.
"""
import os, sys, glob, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..')))
import aim12_single_cell_cytof_reordering as a12

DATA = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'ev_capture_data'))
CACHE = os.path.join(DATA, 'gated_cache')
os.makedirs(CACHE, exist_ok=True)

ANTIGENS = ['CD19','CD20','CD22','CD79B','CD38','CD40','CD24','CD21','CD72','CD48','CD84',
            'CD43','CD44','CD27','CD25','CD23','CD80','CD86','CD83','CXCR4','CXCR5','PD-L1',
            'FasR','TRAIL-R1','IL-4R','CD49D','CD10','HLA-DR','HLA-ABC','HLA-E','CCR4','BTLA',
            'CD45RA','IgM','IgD','IgKap','IgLam']
MAX_CELLS = 25000
rng = np.random.default_rng(11)

def process(path):
    raw = open(path, 'rb').read()
    h = a12.parse_fcs3_header(raw)
    tk, stains = a12.parse_fcs3_text_segment(raw, h)
    tot, npar = int(tk['TOT']), int(tk['PAR'])
    arr = np.frombuffer(raw, dtype='>f4', count=tot*npar, offset=h['data_start']).reshape(tot, npar)
    idx = {m: a12.find_parameter_index_by_stain_substring(stains, m) for m in ANTIGENS + ['CD3','CD45']}
    col = lambda m: arr[:, idx[m]].astype(np.float32)
    dna = arr[:, a12.find_parameter_index_by_stain_substring(stains, '193Ir-DNA')].astype(np.float32)
    via = arr[:, a12.find_parameter_index_by_stain_substring(stains, '195Pt-Viability')].astype(np.float32)
    bead = arr[:, a12.find_parameter_index_by_stain_substring(stains, '140Ce-Beads')].astype(np.float32)
    cd45 = col('CD45')
    cd3, cd19, cd20 = col('CD3'), col('CD19'), col('CD20')
    qc = ((dna > np.percentile(dna,5)) & (dna < np.percentile(dna,99.5))
          & (via < np.percentile(via,95)) & (bead < np.percentile(bead,99))
          & (cd45 > np.percentile(cd45,2)))
    c3  = np.percentile(cd3[qc],50); c19 = np.percentile(cd19[qc],50); c20 = np.percentile(cd20[qc],50)
    B = qc & ((cd19 > c19) | (cd20 > c20)) & (cd3 <= c3)
    T = qc & (cd3 > c3) & (cd19 <= c19) & (cd20 <= c20)
    # third normal compartment: CD45+ leukocytes that are neither B nor T
    # (myeloid, NK, and other non-lymphoid bystanders present in the same tissue)
    O = qc & (~B) & (~T)
    out = {}
    for lab, mask in (('B', B), ('T', T), ('O', O)):
        i = np.flatnonzero(mask)
        if len(i) > MAX_CELLS: i = rng.choice(i, MAX_CELLS, replace=False)
        out[lab] = np.stack([col(m)[i] for m in ANTIGENS], axis=1) if len(i) else np.zeros((0,len(ANTIGENS)),np.float32)
    return out, dict(n_total=tot, n_B=int(B.sum()), n_T=int(T.sum()), n_O=int(O.sum()),
                     dna_med=float(np.median(dna)), qc_frac=float(qc.mean()))

files = (sorted(glob.glob(os.path.join(DATA,'raw_cytof','*.fcs')))
         + sorted(glob.glob(os.path.join(DATA,'raw_cytof_rln','*.fcs'))))
todo = [f for f in files if not os.path.exists(os.path.join(CACHE, os.path.basename(f).replace('.fcs','.npz')))]
print('total files %d, already cached %d, remaining %d' % (len(files), len(files)-len(todo), len(todo)))
limit = int(sys.argv[1]) if len(sys.argv) > 1 else 8
for f in todo[:limit]:
    name = os.path.basename(f).replace('.fcs','')
    o, meta = process(f)
    np.savez_compressed(os.path.join(CACHE, name+'.npz'), B=o['B'], T=o['T'], O=o['O'],
                        antigens=np.array(ANTIGENS), **{k: np.array(v) for k, v in meta.items()})
    print('cached %-24s n=%7d B=%6d T=%6d dna_med=%7.1f' % (name, meta['n_total'], meta['n_B'], meta['n_T'], meta['dna_med']), flush=True)
print('DONE' if len(todo) <= limit else 'MORE REMAINING: %d' % (len(todo)-limit))
