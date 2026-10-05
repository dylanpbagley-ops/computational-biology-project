"""
Extract target-gene rows from the Boiarsky et al. 2022 (GSE193531) NBM/MGUS/SMM/MM
UMI count matrix and save a compact per-cell table (genes x cells -> cells x genes)
for the myeloma healthy-reference check. This is a one-time extraction step; downstream
analysis scripts should read the output of this script, not the raw 55MB matrix.

Boiarsky et al. 2022, Blood Cancer J. (doi:10.1038/s41408-022-00666-1), GEO GSE193531.
CD138-sorted (plasma-cell-enriched) scRNA-seq across NBM (normal bone marrow, n=9 donors),
MGUS, SMM, and MM samples. This matrix has 29,387 total cells and only 8 metadata columns
with no separate cell-type-cluster field beyond disease_stage/cluster/normal_or_neoplastic --
consistent with the whole file already being plasma-cell-sorted/filtered, not whole bone
marrow requiring T/myeloid/NK gating (unlike this project's own CITE-seq myeloma_ceiling.py
gating, which starts from unsorted marrow).
"""
import gzip, csv, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'ev_capture_data', 'healthy_bone_marrow_boiarsky'))
MATRIX = os.path.join(DATA, 'GSE193531_umi-count-matrix.csv.gz')
META = os.path.join(DATA, 'GSE193531_cell-level-metadata.csv.gz')
OUT = os.path.join(DATA, 'extracted_markers.json')

TARGET_GENES = ['TNFRSF17','CD38','SDC1','SLAMF7','GPRC5D','CD19','CD3D','CD3E','CD14','NCAM1','MS4A1']

print('Reading metadata...')
with gzip.open(META, 'rt') as f:
    r = csv.reader(f)
    header = next(r)
    meta_rows = list(r)
print('  %d metadata rows, columns: %s' % (len(meta_rows), header))

cell_ids = [row[0] for row in meta_rows]
disease_stage = {row[0]: row[2] for row in meta_rows}
cluster = {row[0]: row[6] for row in meta_rows}
norm_neo = {row[0]: row[7] for row in meta_rows}

print('Scanning matrix for target genes (single pass over %s)...' % MATRIX)
gene_rows = {}
header_cells = None
n_lines = 0
with gzip.open(MATRIX, 'rt') as f:
    r = csv.reader(f)
    header_cells = next(r)[1:]  # drop leading blank col name
    remaining = set(TARGET_GENES)
    for row in r:
        n_lines += 1
        g = row[0]
        if g in remaining:
            gene_rows[g] = [float(x) for x in row[1:]]
            remaining.discard(g)
            if not remaining:
                break
    print('  scanned %d gene rows; found %d/%d target genes' % (n_lines, len(gene_rows), len(TARGET_GENES)))
    missing = set(TARGET_GENES) - set(gene_rows.keys())
    if missing:
        print('  MISSING:', missing)

assert header_cells == cell_ids, 'matrix column order does not match metadata row order!'

print('Assembling per-cell records...')
out = {
    'genes': list(gene_rows.keys()),
    'cell_ids': cell_ids,
    'disease_stage': [disease_stage[c] for c in cell_ids],
    'cluster': [cluster[c] for c in cell_ids],
    'normal_or_neoplastic': [norm_neo[c] for c in cell_ids],
    'expr': {g: gene_rows[g] for g in gene_rows},
}
with open(OUT, 'w') as f:
    json.dump(out, f)
sz = os.path.getsize(OUT)
print('Wrote %s (%.1f MB)' % (OUT, sz/1e6))
