"""
Compute the confinement measure for EVERY gene in the Human Protein Atlas single-cell
consensus, not just the 40-target panel. This is the artifact the tool ships: a genome-wide
table of how confined each human gene is across 154 normal cell types, plus which cell type
it calls home and whether that home is an expendable (haematopoietic) lineage.

Outputs:
  confinement_table.tsv   full table, one row per gene (human-readable, for the repository)
  confinement_web.json    compact form for the interactive tool (cell types stored once,
                          referenced by index)
"""
import os, io, json, zipfile, collections
import numpy as np
import target_labels as TL

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'ev_capture_data', 'normal_tissue_atlas'))
Z = os.path.join(D, 'rna_single_cell_type.tsv.zip')

z = zipfile.ZipFile(Z)
name = z.namelist()[0]

cell_types = []
ct_index = {}
tot = collections.defaultdict(float)
top = collections.defaultdict(lambda: (-1.0, -1))
by_gene = collections.defaultdict(list)
ens = {}

with z.open(name) as f:
    fh = io.TextIOWrapper(f, 'utf-8')
    hdr = fh.readline().rstrip('\n').split('\t')
    gi, ci, vi = hdr.index('Gene name'), hdr.index('Cell type'), hdr.index('nCPM')
    ei = hdr.index('Gene') if 'Gene' in hdr else None
    for line in fh:
        p = line.rstrip('\n').split('\t')
        try:
            v = float(p[vi])
        except (ValueError, IndexError):
            continue
        g = p[gi]
        c = p[ci]
        j = ct_index.get(c)
        if j is None:
            j = ct_index[c] = len(cell_types)
            cell_types.append(c)
        tot[g] += v
        by_gene[g].append(v)
        if v > top[g][0]:
            top[g] = (v, j)
        if ei is not None and g not in ens:
            ens[g] = p[ei]

print('cell types: %d ; genes: %d' % (len(cell_types), len(tot)))

rows = []
for g, t in tot.items():
    if t <= 0:
        continue
    mx, j = top[g]
    v = np.asarray(by_gene[g], float)
    rows.append((g, ens.get(g, ''), mx / t, int((v > 0.10 * mx).sum()), j))
rows.sort(key=lambda r: -r[2])

with io.open(os.path.join(HERE, 'confinement_table.tsv'), 'w', encoding='utf-8', newline='\n') as fh:
    fh.write('gene\tensembl\ttop_share\tn_cell_types\thome_cell_type\thome_is_expendable_lineage\n')
    for g, e, s, n, j in rows:
        fh.write('%s\t%s\t%.6f\t%d\t%s\t%s\n'
                 % (g, e, s, n, cell_types[j], 'yes' if TL.haematopoietic(cell_types[j]) else 'no'))

json.dump({'cell_types': cell_types,
           'haem': [bool(TL.haematopoietic(c)) for c in cell_types],
           'genes': [[g, round(s, 4), n, j] for g, e, s, n, j in rows]},
          io.open(os.path.join(HERE, 'confinement_web.json'), 'w', encoding='utf-8'),
          separators=(',', ':'))

print('rank 1  : %s  %.3f  (%s)' % (rows[0][0], rows[0][2], cell_types[rows[0][4]]))
for probe in ('CD19', 'MS4A1', 'CD47', 'ERBB2', 'CLDN18'):
    for i, r in enumerate(rows):
        if r[0] == probe:
            print('%-8s rank %5d of %d   confinement %.3f   %d cell types   home %s'
                  % (probe, i + 1, len(rows), r[2], r[3], cell_types[r[4]]))
            break
for f in ('confinement_table.tsv', 'confinement_web.json'):
    print('%-24s %8.1f KB' % (f, os.path.getsize(os.path.join(HERE, f)) / 1024))
