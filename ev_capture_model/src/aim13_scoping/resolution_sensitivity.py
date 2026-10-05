"""
DOES THE MEASURE NEED CELL-TYPE RESOLUTION? A controlled test.

Four references have now been run, and the results line up with how finely each one resolves
normal tissue rather than with whether it measures RNA or protein:

    reference                          groups   separation      p
    RNA, single-cell consensus            154        0.868   0.0003
    PROTEIN, IHC, tissue x cell type       76        0.785   0.0104
    RNA, bulk tissue consensus             51        0.726   0.0244
    PROTEIN, mass spectrometry             20        0.650   0.1148   (H5 NOT MET)

The obvious self-serving story is "the MS reference is just too coarse". That is exactly the
kind of excuse a pre-registration exists to stop, so it is not asserted here -- it is tested.

THE TEST. Take the 154-cell-type RNA reference, the one where the result is strongest, and
COARSEN IT ON PURPOSE: randomly merge the 154 cell types into K groups, summing expression
within each group, then recompute the separation exactly as before. Repeat many times per K.
Nothing about the data changes except resolution.

If separation decays with K and lands near the observed value at K = 20, then resolution
explains the pattern and the MS result is a statement about that reference's granularity.
If separation stays high at K = 20, resolution does NOT explain it, the excuse is dead, and
the honest conclusion is that the finding does not transfer to quantitative protein.
"""
import os, io, json, zipfile, collections
import numpy as np
import target_labels as TL
import wholebody_v2 as WB

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'ev_capture_data', 'normal_tissue_atlas'))
GENES = {t[0] for t in TL.TARGETS}
LABEL = {t[0]: t[3] for t in TL.TARGETS}
N_REPEATS = 60

z = zipfile.ZipFile(os.path.join(D, 'rna_single_cell_type.tsv.zip'))
vals = collections.defaultdict(dict)
with z.open(z.namelist()[0]) as f:
    fh = io.TextIOWrapper(f, 'utf-8')
    hdr = fh.readline().rstrip('\n').split('\t')
    gi, ci, vi = hdr.index('Gene name'), hdr.index('Cell type'), hdr.index('nCPM')
    for line in fh:
        p = line.rstrip('\n').split('\t')
        if p[gi] in GENES:
            try:
                vals[p[gi]][p[ci]] = float(p[vi])
            except ValueError:
                pass

cell_types = sorted(next(iter(vals.values())))
n_ct = len(cell_types)
genes = [g for g in vals if LABEL[g] in ('safe', 'toxic')]
M = np.array([[vals[g].get(c, 0.0) for c in cell_types] for g in genes], float)
print('reference: %d cell types, %d labelled targets' % (n_ct, len(genes)))

OBSERVED = [(154, 0.868, 'RNA single-cell (the primary reference)'),
            (76, 0.785, 'PROTEIN immunohistochemistry'),
            (51, 0.726, 'RNA bulk tissue'),
            (20, 0.650, 'PROTEIN mass spectrometry')]

print('\ncoarsening the 154-cell-type reference by random merging, %d repeats per K:\n' % N_REPEATS)
print('  %5s  %-28s %-28s' % ('K', 'separation  mean [min, max]', 'observed at this K'))
out = {}
for K in (154, 120, 100, 76, 60, 51, 40, 30, 20, 14, 10, 6):
    seps = []
    for r in range(N_REPEATS if K < n_ct else 1):
        rng = np.random.default_rng(1000 + r)
        assign = rng.integers(0, K, n_ct)
        # guarantee every group is used so K is the true group count
        assign[rng.permutation(n_ct)[:K]] = np.arange(K)
        G = np.zeros((len(genes), K))
        np.add.at(G.T, assign, M.T)
        tot = G.sum(1)
        top = G.max(1)
        share = np.where(tot > 0, top / np.maximum(tot, 1e-12), 0.0)
        sc = {g: float(s) for g, s in zip(genes, share)}
        seps.append(WB.sep([sc[g] for g in sc if LABEL[g] == 'safe'],
                           [sc[g] for g in sc if LABEL[g] == 'toxic']))
    obs = next((f'{v:.3f}  {n}' for k, v, n in OBSERVED if k == K), '')
    print('  %5d  %.3f  [%.3f, %.3f]          %s' % (K, np.mean(seps), min(seps), max(seps), obs))
    out[K] = {'mean': float(np.mean(seps)), 'min': float(min(seps)), 'max': float(max(seps))}

k20 = out[20]
print("""
READING.""")
if k20['min'] <= 0.650 <= k20['max'] or abs(k20['mean'] - 0.650) < 0.06:
    print("""  At K = 20 the coarsened RNA reference gives %.3f [%.3f, %.3f], and the mass-spectrometry
  reference -- 20 tissues -- gave 0.650. The same data, measured the same way, loses the same
  amount of separation purely from being pooled into 20 groups. So the MS failure is consistent
  with GRANULARITY, not with RNA-versus-protein.

  This does NOT rescue H5, and H5 is still reported as NOT MET on that reference. What it says
  is narrower and checkable: the measure needs cell-type resolution to work, and any reference
  that pools normal tissue into ~20 bulk samples cannot test it. That is now a stated scope
  condition of the framework rather than a hopeful excuse.""" % (k20['mean'], k20['min'], k20['max']))
else:
    print("""  At K = 20 the coarsened RNA reference still gives %.3f [%.3f, %.3f], well above the 0.650
  the mass-spectrometry reference produced. Granularity does NOT explain the failure. The
  honest conclusion is that the finding does not transfer to quantitative protein measurement,
  and that limitation goes in the write-up as stated.""" % (k20['mean'], k20['min'], k20['max']))

json.dump({'coarsening': out, 'observed': OBSERVED},
          open(os.path.join(HERE, 'resolution_sensitivity.json'), 'w'), indent=1)
print('\nwrote resolution_sensitivity.json')
