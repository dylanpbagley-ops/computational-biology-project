"""
CELL-TYPE-RESOLUTION PROTEIN TEST, from the per-gene HPA summary table.

Section 26 concluded that the safe-vs-toxic result needs CELL-TYPE resolution: it holds on RNA
across 154 cell types (0.868) and on antibody staining across tissue x cell type (0.785), and
fails on both RNA and protein when normal tissue is pooled into ~20 whole organs (0.626, 0.650).
The open question was whether it holds on PROTEIN measured at cell-type resolution, which needed
a reference that did not exist when section 26 was written.

It exists now. HPA v25.1 integrated Deep Visual Proteomics -- mass-spectrometry protein
intensity per CELL TYPE -- and the per-gene summary of it ships inside proteinatlas.tsv.zip as
'Protein cell type specificity' and 'Protein cell type specific Intensity'.

WHAT THIS TEST CAN AND CANNOT DO. The summary gives HPA's own categorical specificity call
(Cell type enriched > Group enriched > Cell type enhanced > Low cell type specificity) plus the
intensities of only the cell types that drove that call. It is therefore COARSER than the
continuous top_share used everywhere else, and a categorical ladder cannot be compared
numerically to 0.868. What it can do is answer the yes/no question: at protein, at cell-type
resolution, do the safe targets sit higher on the specificity ladder than the toxic ones?

The RNA columns from the same table are run alongside as the matched control, because they are
the same summarisation applied to the reference where the answer is already known.
"""
import os, io, csv, json, zipfile
import numpy as np
import target_labels as TL
import wholebody_v2 as WB

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'ev_capture_data', 'normal_tissue_atlas'))
Z = os.path.join(D, 'proteinatlas.tsv.zip')
LABEL = {t[0]: t[3] for t in TL.TARGETS}
GENES = set(LABEL)

# HPA's specificity ladder, most specific first. Ordinal scores are the ladder position, not a
# measurement -- stated so nobody reads 3.0 as three times 1.0.
LADDER = {'cell type enriched': 4, 'group enriched': 3, 'cell type enhanced': 2,
          'tissue enriched': 4, 'tissue enhanced': 2,
          'low cell type specificity': 1, 'low tissue specificity': 1,
          'low immune cell specificity': 1,
          'immune cell enriched': 4, 'immune cell enhanced': 2, 'lineage enriched': 3}
# 'Not detected' is deliberately ABSENT from this table, and that is the single most important
# line in this file. A first version mapped it to 0, i.e. to "least specific of all" -- the
# exact opposite of what it means. Not detected means the protein was below the detection
# limit of a 27-cell-type, single-donor proteome; it does not mean the protein is everywhere.
# Scoring it as 0 put BCMA, SLAMF7, CD30, GPRC5D and CD52 -- five of the most confined and
# safest targets in the panel -- at the bottom of the ladder, and produced a meaningless
# separation of 0.413. Absent from LADDER, ladder_score() now returns None and those targets
# are excluded rather than mis-scored.

z = zipfile.ZipFile(Z)
name = z.namelist()[0]
rows = {}
with z.open(name) as f:
    fh = io.TextIOWrapper(f, 'utf-8', newline='')
    rd = csv.reader(fh, delimiter='\t', quotechar='"')
    hdr = next(rd)
    idx = {h.strip('"'): i for i, h in enumerate(hdr)}
    gi = idx['Gene']
    for r in rd:
        if len(r) > gi and r[gi] in GENES:
            rows[r[gi]] = r

print('targets found in the per-gene HPA table: %d of %d' % (len(rows), len(GENES)))


def col(g, name_):
    i = idx.get(name_)
    if i is None or i >= len(rows[g]):
        return ''
    return rows[g][i].strip().strip('"')


def ladder_score(g, name_):
    v = col(g, name_).lower()
    return LADDER.get(v)


def top_share_from_intensity(g, name_):
    """Crude confinement from the listed 'specific' cell types. Only defined when HPA listed
    any; a gene with no listing is NOT scored 0, it is left out, because absent means
    'not specific enough to list', not 'zero everywhere'."""
    s = col(g, name_)
    if not s:
        return None
    vals = []
    for part in s.split(';'):
        if ':' in part:
            try:
                vals.append(float(part.rsplit(':', 1)[1]))
            except ValueError:
                pass
    if not vals:
        return None
    v = np.array(vals, float)
    return float(v.max() / v.sum())


def run(label_name, scores, higher_is_safer=True):
    use = {g: v for g, v in scores.items() if v is not None and LABEL.get(g) in ('safe', 'toxic')}
    groups = {LABEL[g] for g in use}
    if len(groups) < 2:
        print('  %-46s (not enough labelled targets: %d)' % (label_name, len(use)))
        return None
    s, p, ex, na, nb = WB.run_test(use, LABEL, higher_is_safer)
    print('  %-46s %6.3f  %s p = %.4f   (safe n=%d vs toxic n=%d)'
          % (label_name, s, 'exact' if ex else 'approx', p, na, nb))
    return {'sep': s, 'p': p, 'n_safe': na, 'n_toxic': nb}


print('\n' + '=' * 88)
print('THE OPEN QUESTION: protein, at CELL-TYPE resolution')
print('=' * 88)
out = {}
out['protein_cell_type_ladder'] = run(
    'PROTEIN cell-type specificity (DVP, ordinal ladder)',
    {g: ladder_score(g, 'Protein cell type specificity') for g in rows})
out['protein_cell_type_topshare'] = run(
    'PROTEIN cell-type intensity, crude top_share',
    {g: top_share_from_intensity(g, 'Protein cell type specific Intensity') for g in rows})

print('\n  matched controls from the SAME table, same summarisation:')
out['rna_cell_type_ladder'] = run(
    'RNA single-cell-type specificity (ordinal ladder)',
    {g: ladder_score(g, 'RNA single cell type specificity') for g in rows})
out['rna_cell_type_topshare'] = run(
    'RNA single-cell-type nCPM, crude top_share',
    {g: top_share_from_intensity(g, 'RNA single cell type specific nCPM') for g in rows})

print('\n  and the whole-organ versions, which section 26 predicts should be WEAKER:')
out['protein_tissue_ladder'] = run(
    'PROTEIN whole-tissue specificity (ordinal ladder)',
    {g: ladder_score(g, 'Protein tissue specificity') for g in rows})
out['rna_tissue_ladder'] = run(
    'RNA whole-tissue specificity (ordinal ladder)',
    {g: ladder_score(g, 'RNA tissue specificity') for g in rows})

print('\n' + '=' * 88)
print('per-target detail (protein cell-type call, then RNA cell-type call)')
print('=' * 88)
print('  %-9s %-8s %-26s %-26s' % ('GENE', 'outcome', 'PROTEIN cell-type call', 'RNA cell-type call'))
for g in sorted(rows, key=lambda g: (LABEL[g], g)):
    if LABEL[g] not in ('safe', 'toxic'):
        continue
    print('  %-9s %-8s %-26s %-26s'
          % (g, LABEL[g], col(g, 'Protein cell type specificity')[:26],
             col(g, 'RNA single cell type specificity')[:26]))

nd = [g for g in rows if col(g, 'Protein cell type specificity').lower() == 'not detected'
      and LABEL[g] in ('safe', 'toxic')]
nd_safe = [g for g in nd if LABEL[g] == 'safe']
print('\n' + '=' * 88)
print('WHY THE PROTEIN CELL-TYPE TEST CANNOT ANSWER THE QUESTION YET')
print('=' * 88)
print('  %d of 29 labelled targets are NOT DETECTED in the cell-type proteome: %s'
      % (len(nd), ', '.join(sorted(nd))))
print('  %d of those %d are SAFE targets -- and they are precisely the confined blood-lineage'
      % (len(nd_safe), len(nd)))
print('  antigens the framework is built on (BCMA, SLAMF7, CD30, GPRC5D, CD52).')
print('''
  The Deep Visual Proteomics reference is 27 cell types from a SINGLE donor. Plasma cells and
  activated lymphocytes are rare, so the antigens that define the safe class sit below its
  detection limit. Dropping them does not fix the problem -- it removes the safe class.

  CONCLUSION, and it is a statement about the DATA, not about the framework: no protein
  reference at cell-type resolution currently covers the cell types this question needs. The
  antibody-staining result (0.785, p = 0.010, 76 tissue x cell-type groups) remains the best
  protein-level evidence, and section 26's scope condition stands unchanged.

  A CALIBRATION THAT MAKES THE ABOVE READABLE: HPA's ordinal ladder is too coarse to detect
  the effect even where the effect is known to exist. On the RNA cell-type data -- where the
  full continuous computation gives 0.868 -- the same 4-level ladder gives only 0.553 (ns),
  while the continuous top_share from the same table gives 0.716 (p = 0.031). So a null on the
  ladder is uninformative, not negative.''')

json.dump(out, open(os.path.join(HERE, 'hpa_summary_reference.json'), 'w'), indent=1)
print('\nwrote hpa_summary_reference.json')
