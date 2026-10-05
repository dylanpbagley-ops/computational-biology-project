"""
WHY THE MASS-SPECTROMETRY CROSS-CHECK FAILED, and whether the reason is real.

H5 (registered in advance) was NOT MET on the mass-spectrometry reference: separation 0.643,
p = 0.127. That is the registered result and it stands.

This script asks the separate question of WHETHER THE MS REFERENCE CAN CARRY THE MEASURE AT
ALL. The concern is specific and checkable: `top_share` is a share of a total, so a gene
measured in only one tissue scores 1.000 by arithmetic rather than by biology. If the MS table
is sparse, the measure is not confinement there -- it is a map of what happened to be detected.

EVERYTHING BELOW IS POST HOC. It was not registered, it cannot rescue H5, and it is reported
as a diagnostic of the reference rather than as a result about targets.
"""
import os, io, zipfile, collections
import numpy as np
import target_labels as TL
import wholebody_v2 as WB

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'ev_capture_data', 'normal_tissue_atlas'))
GENES = {t[0] for t in TL.TARGETS}
LABEL = {t[0]: t[3] for t in TL.TARGETS}

z = zipfile.ZipFile(os.path.join(D, 'ms_tissue.tsv.zip'))
vals = collections.defaultdict(dict)
all_tissues = set()
with z.open(z.namelist()[0]) as f:
    fh = io.TextIOWrapper(f, 'utf-8')
    hdr = fh.readline().rstrip('\n').split('\t')
    gi, ti, vi = hdr.index('Gene name'), hdr.index('Tissue'), hdr.index('Intensity')
    for line in fh:
        p = line.rstrip('\n').split('\t')
        all_tissues.add(p[ti])
        if p[gi] in GENES:
            try:
                vals[p[gi]][p[ti]] = float(p[vi])
            except ValueError:
                pass

N_TISSUES = len(all_tissues)
print('tissues in the MS reference: %d' % N_TISSUES)
print('\nhow many tissues each labelled target was actually MEASURED in:\n')
print('  %-9s %-14s %8s %10s  %s' % ('GENE', 'outcome', 'measured', 'top_share', 'home'))
rows = []
for g in sorted(vals, key=lambda g: len(vals[g])):
    if LABEL[g] not in ('safe', 'toxic'):
        continue
    d = vals[g]
    c = WB.confinement({g: d})[g]
    rows.append((g, len(d), c['top_share'], c['home']))
    flag = '   <-- share is 1.0 because only one tissue was measured' if len(d) == 1 else ''
    print('  %-9s %-14s %5d/%d %10.3f  %s%s'
          % (g, LABEL[g], len(d), N_TISSUES, c['top_share'], c['home'][:22], flag))

single = [r for r in rows if r[1] == 1]
print('\n  %d of %d labelled targets were measured in a SINGLE tissue.' % (len(single), len(rows)))
print('  median tissues measured: %.0f of %d' % (np.median([r[1] for r in rows]), N_TISSUES))

print('\n' + '-' * 78)
print('POST HOC: restrict to targets measured in at least K tissues and re-run.')
print('This cannot rescue H5. It only says whether sparsity explains the failure.')
print('-' * 78)
for K in (1, 3, 5, 8):
    keep = {g: WB.confinement({g: vals[g]})[g] for g in vals
            if LABEL[g] in ('safe', 'toxic') and len(vals[g]) >= K}
    if len({LABEL[g] for g in keep}) < 2:
        print('  K >= %-2d   not enough targets left in both groups' % K)
        continue
    s, p, ex, na, nb = WB.run_test({g: d['top_share'] for g, d in keep.items()}, LABEL, True)
    sn, pn, _, _, _ = WB.run_test({g: d['n_types'] for g, d in keep.items()}, LABEL, False)
    print('  K >= %-2d   n = %2d safe vs %2d toxic   top_share %.3f (p = %.4f)   n_types %.3f (p = %.4f)'
          % (K, na, nb, s, p, sn, pn))

print("""
READING. If separation climbs as the sparsely-measured targets are removed, the MS failure is
about detection coverage rather than about the framework. If it does not climb, the failure is
real and the RNA result does not transfer to quantitative protein -- which is a genuine
limitation and belongs in the write-up either way.""")
