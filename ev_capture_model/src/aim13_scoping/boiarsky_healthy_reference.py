"""
Myeloma-side healthy-reference check, structurally paralleling section 16's
CD19/DLBCL healthy-donor-PBMC test but for BCMA/myeloma using Boiarsky et al.'s
independent, external NBM (normal bone marrow) scRNA-seq cohort.

Question: in genuinely healthy, external, independent human tissue -- not this
project's own within-patient comparisons -- is BCMA (TNFRSF17) restricted to
plasma cells, or is it detectably expressed on lineages a therapy would need to
spare (T, myeloid, NK, B)?

Data: Boiarsky et al. 2022 (GSE193531), extracted_markers.json (built by
extract_boiarsky_markers.py from the raw GEO matrix). This matrix is CD138-sorted
plasma-cell-enriched, not whole bone marrow -- confirmed below by checking that
canonical non-plasma-cell markers (CD3D/CD3E T-cell, CD14 myeloid, MS4A1/CD19
B-cell, NCAM1 NK) are nearly silent even in cells nominally in this file, while
plasma cell markers (CD38, SDC1/CD138) are near-ubiquitous. This is an important
scope limitation, stated plainly: this dataset can show whether BCMA is confined
to the plasma-cell compartment among plasma-cell-like/sorted cells, but it is NOT
a whole-marrow dataset that can independently re-derive the T/myeloid/NK/B
"other lineage" negative controls the way the CD19/PBMC comparison in section 16
could (that dataset was whole blood, unsorted). Reported honestly either way.
"""
import json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'ev_capture_data', 'healthy_bone_marrow_boiarsky'))
IN = os.path.join(DATA, 'extracted_markers.json')

with open(IN) as f:
    D = json.load(f)

expr = {g: np.array(v) for g, v in D['expr'].items()}
disease_stage = np.array(D['disease_stage'])
cluster = np.array(D['cluster'])
norm_neo = np.array(D['normal_or_neoplastic'])
n_cells = len(D['cell_ids'])
print('Total cells in file: %d' % n_cells)
print('Genes extracted: %s' % sorted(expr.keys()))
print()

# ---- Step 1: confirm this is plasma-cell-sorted, not whole marrow ----
print('=== Step 1: confirm plasma-cell-sorted composition (all %d cells) ===' % n_cells)
pc_markers = ['CD38', 'SDC1']
other_markers = ['CD3D', 'CD3E', 'CD14', 'MS4A1', 'CD19', 'NCAM1']
for g in pc_markers + other_markers:
    v = expr[g]
    pct_pos = 100.0 * np.mean(v > 0)
    tag = '(plasma-cell marker)' if g in pc_markers else '(non-plasma-cell lineage marker)'
    print('  %-10s %-28s %%>0 = %6.2f%%   median(all)=%.2f  median(>0)=%.2f' % (
        g, tag, pct_pos, np.median(v), np.median(v[v > 0]) if pct_pos > 0 else float('nan')))
print()

# ---- Step 2: restrict to NBM (genuinely healthy donor) cells ----
is_nbm = disease_stage == 'NBM'
n_nbm = int(is_nbm.sum())
print('=== Step 2: restrict to NBM (healthy-donor) cells: n=%d ===' % n_nbm)
print('  NBM cluster label breakdown:')
for lab, cnt in zip(*np.unique(cluster[is_nbm], return_counts=True)):
    print('    cluster=%-10s n=%d' % (lab, cnt))
print('  NBM normal_or_neoplastic breakdown:')
for lab, cnt in zip(*np.unique(norm_neo[is_nbm], return_counts=True)):
    print('    %-12s n=%d' % (lab, cnt))
print()

# ---- Step 3: BCMA (TNFRSF17) expression in healthy-donor plasma cells vs everything else in the file ----
print('=== Step 3: TNFRSF17 (BCMA) in NBM cells specifically ===')
bcma = expr['TNFRSF17']
bcma_nbm = bcma[is_nbm]
print('  NBM cells: %%>0 = %.2f%%, median = %.2f, mean = %.3f, p90 = %.2f, max = %.0f' % (
    100.0*np.mean(bcma_nbm>0), np.median(bcma_nbm), np.mean(bcma_nbm), np.percentile(bcma_nbm,90), bcma_nbm.max()))
print('  cf. CD38 in same NBM cells: %%>0 = %.2f%%, median = %.2f' % (
    100.0*np.mean(expr['CD38'][is_nbm]>0), np.median(expr['CD38'][is_nbm])))
print('  cf. SDC1/CD138 in same NBM cells: %%>0 = %.2f%%, median = %.2f' % (
    100.0*np.mean(expr['SDC1'][is_nbm]>0), np.median(expr['SDC1'][is_nbm])))
print()

print('=== Step 4: BCMA in the (few) non-plasma-cell-marker-positive NBM cells, if any exist ===')
non_pc_like = is_nbm & (expr['CD38'] == 0) & (expr['SDC1'] == 0)
n_non_pc = int(non_pc_like.sum())
print('  NBM cells with CD38==0 AND SDC1==0 (i.e. NOT plasma-cell-marker-positive): n=%d (%.2f%% of NBM)' % (
    n_non_pc, 100.0*n_non_pc/n_nbm))
if n_non_pc > 0:
    bcma_non_pc = bcma[non_pc_like]
    print('    BCMA in these cells: %%>0 = %.2f%%, median = %.2f, max = %.0f' % (
        100.0*np.mean(bcma_non_pc>0), np.median(bcma_non_pc), bcma_non_pc.max()))
    for g in other_markers:
        v = expr[g][non_pc_like]
        print('    %-8s in these cells: %%>0 = %.2f%%, median = %.2f' % (g, 100.0*np.mean(v>0), np.median(v)))
print()

print('=== Step 5: GPRC5D and SLAMF7 (other myeloma antigens) in NBM, for context ===')
for g in ['GPRC5D', 'SLAMF7']:
    v = expr[g][is_nbm]
    print('  %-8s in NBM: %%>0 = %.2f%%, median = %.2f, p90 = %.2f' % (g, 100.0*np.mean(v>0), np.median(v), np.percentile(v,90)))

print()
print('=== Step 6: is the CD38-/SDC1- subset just low-depth plasma cells (dropout), or a real other lineage? ===')
# reload n_counts/n_genes from metadata to check sequencing depth in that subset
import gzip, csv
META = os.path.join(DATA, 'GSE193531_cell-level-metadata.csv.gz')
with gzip.open(META, 'rt') as f:
    r = csv.reader(f)
    header = next(r)
    rows = list(r)
ncounts = np.array([float(row[5]) for row in rows])
ngenes = np.array([float(row[3]) for row in rows])
cell_ids_meta = [row[0] for row in rows]
assert cell_ids_meta == D['cell_ids']

nbm_ncounts = ncounts[is_nbm]
nbm_ngenes = ngenes[is_nbm]
print('  NBM overall: median n_counts=%.0f, median n_genes=%.0f' % (np.median(nbm_ncounts), np.median(nbm_ngenes)))
print('  NBM CD38+/SDC1+ (n=%d): median n_counts=%.0f, median n_genes=%.0f' % (
    int((~non_pc_like[is_nbm]).sum()) if False else int((is_nbm & ~non_pc_like).sum()),
    np.median(ncounts[is_nbm & ~non_pc_like]), np.median(ngenes[is_nbm & ~non_pc_like])))
print('  NBM CD38-/SDC1- (n=%d): median n_counts=%.0f, median n_genes=%.0f' % (
    n_non_pc, np.median(ncounts[non_pc_like]), np.median(ngenes[non_pc_like])))
