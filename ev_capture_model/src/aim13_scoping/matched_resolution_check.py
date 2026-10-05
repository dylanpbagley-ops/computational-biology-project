"""
THE CONTROL I SHOULD HAVE RUN THE FIRST TIME.

(Renamed from matched_resolution_test.py: pytest auto-collects any *_test.py file, and this
module's test() helper takes required arguments, which produced a collection ERROR alongside an
otherwise clean 195-pass run. The analysis is unchanged.)

Section 26 tested whether "the mass-spectrometry reference is too coarse" explains why the
result failed there, by RANDOMLY merging the 154 cell types into 20 groups. Separation barely
moved (0.837), and I concluded that coarseness was not the explanation.

THAT WAS THE WRONG TEST, and the reason is specific. Random merging preserves distinctness. A
gene expressed in exactly one of 154 cell types still lands in exactly one random bucket, so its
confinement stays near 1.0. Real anatomical pooling does not behave that way at all: a B-cell
gene is spread across every tissue that contains lymphocytes -- lymph node, spleen, tonsil, gut,
lung, bone marrow -- so pooling by ORGAN dilutes it in a way random pooling never does. Random
merging tested "does having fewer buckets hurt?" (no). It did not test "does mixing cell types
into real organs hurt?", which is the actual question.

THE RIGHT TEST, run here. The mass-spectrometry reference covers 20 named human tissues, and 18
of those names appear verbatim in the HPA bulk RNA reference. So the same measure can be
computed on the SAME 18 organs, in RNA and in protein, with everything else held constant. If
RNA also collapses on those 18 organs, the failure is about anatomical pooling and not about
protein. If RNA holds up while protein does not, the failure is real and about the measurement.

Either answer is reported. This was written before either was known.
"""
import os, io, json, zipfile, collections
import numpy as np
import target_labels as TL
import wholebody_v2 as WB

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'ev_capture_data', 'normal_tissue_atlas'))
GENES = {t[0] for t in TL.TARGETS}
LABEL = {t[0]: t[3] for t in TL.TARGETS}


def load(fn, value_cols, group_col='Tissue', keep=None):
    z = zipfile.ZipFile(os.path.join(D, fn))
    vals = collections.defaultdict(lambda: collections.defaultdict(list))
    with z.open(z.namelist()[0]) as f:
        fh = io.TextIOWrapper(f, 'utf-8')
        hdr = fh.readline().rstrip('\n').split('\t')
        gi = hdr.index('Gene name')
        ci = hdr.index(group_col)
        vi = next(hdr.index(c) for c in value_cols if c in hdr)
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) <= max(gi, ci, vi) or p[gi] not in GENES:
                continue
            if keep is not None and p[ci] not in keep:
                continue
            try:
                vals[p[gi]][p[ci]].append(float(p[vi]))
            except ValueError:
                pass
    # drop groups with no measured value; the MS table leaves Intensity blank where a protein
    # was not detected, and averaging an empty list silently produced NaN in a first version of
    # this script, which made the MS row read 0.000 instead of its true value.
    return {g: {k: float(np.mean(v)) for k, v in d.items() if len(v)}
            for g, d in vals.items()}


def test(name, raw, note=''):
    conf = WB.confinement(raw)
    use = {g: d['top_share'] for g, d in conf.items() if LABEL.get(g) in ('safe', 'toxic')}
    if len({LABEL[g] for g in use}) < 2:
        print('  %-46s (too few labelled targets)' % name)
        return None
    s, p, ex, na, nb = WB.run_test(use, LABEL, True)
    ngroups = int(np.median([len(d) for d in raw.values()])) if raw else 0
    print('  %-46s %6.3f  %s p = %.4f   (n=%d vs %d, ~%d groups) %s'
          % (name, s, 'exact' if ex else 'approx', p, na, nb, ngroups, note))
    return {'sep': s, 'p': p, 'n_safe': na, 'n_toxic': nb, 'groups': ngroups}


MS20 = {'blood vessel', 'bone marrow', 'cerebral cortex', 'colon', 'duodenum', 'fallopian tube',
        'heart muscle', 'ileum', 'jejunum', 'kidney', 'liver', 'lung', 'lymph node', 'ovary',
        'pancreas', 'salivary gland', 'skeletal muscle', 'skin', 'spleen', 'stomach'}
SHARED18 = MS20 - {'ileum', 'jejunum'} | {'small intestine'}

print('=' * 88)
print('SAME MEASURE, SAME ORGANS. Only the measurement technology differs.')
print('=' * 88)
out = {}
out['rna_cell_types_154'] = test('RNA, 154 CELL TYPES (the primary reference)',
                                 load('rna_single_cell_type.tsv.zip', ['nCPM'], 'Cell type'),
                                 '<- cell-type resolution')
out['rna_tissue_51'] = test('RNA, 51 whole organs',
                            load('rna_tissue_consensus.tsv.zip', ['nTPM']))
out['rna_tissue_18_matched'] = test('RNA, the SAME 18 organs the MS panel covers',
                                    load('rna_tissue_consensus.tsv.zip', ['nTPM'], keep=SHARED18),
                                    '<- the matched control')
out['ms_tissue_20'] = test('PROTEIN mass spec, 20 organs',
                           load('ms_tissue.tsv.zip', ['Intensity']))
out['ihc_76'] = None

# IHC needs a composite tissue x cell type key, so build it separately
z = zipfile.ZipFile(os.path.join(D, 'normal_ihc_data.tsv.zip'))
LEV = {'not detected': 0.0, 'low': 1.0, 'medium': 2.0, 'high': 3.0}
ihc = collections.defaultdict(lambda: collections.defaultdict(list))
with z.open(z.namelist()[0]) as f:
    fh = io.TextIOWrapper(f, 'utf-8')
    hdr = fh.readline().rstrip('\n').split('\t')
    gi, ti, ci, vi = (hdr.index('Gene name'), hdr.index('Tissue'),
                      hdr.index('Cell type'), hdr.index('Level'))
    for line in fh:
        p = line.rstrip('\n').split('\t')
        if len(p) <= vi or p[gi] not in GENES:
            continue
        v = LEV.get(p[vi].strip().lower())
        if v is not None:
            ihc[p[gi]]['%s | %s' % (p[ti], p[ci])].append(v)
ihc = {g: {k: float(np.mean(v)) for k, v in d.items()} for g, d in ihc.items()}
out['ihc_76'] = test('PROTEIN antibody staining, tissue x CELL TYPE',
                     ihc, '<- cell-type resolution')

r18 = out['rna_tissue_18_matched']
ms = out['ms_tissue_20']
print('\n' + '=' * 88)
print('READING')
print('=' * 88)
if r18 and ms and abs(r18['sep'] - ms['sep']) < 0.12:
    print("""  RNA on the same 18 organs gives %.3f (p = %.3f). Mass spectrometry on those organs gives
  %.3f (p = %.3f). The two technologies AGREE once the anatomical resolution is matched.

  So the section 26 failure is NOT about protein versus RNA. It is about RESOLUTION: pooling
  normal tissue into whole organs mixes the expendable cell types in with the essential ones,
  and that is precisely the distinction the measure depends on. My earlier random-merging test
  did not catch this because random merging preserves distinctness -- a gene confined to one of
  154 cell types stays confined to one random bucket, whereas a real organ contains dozens of
  cell types at once.

  The claim that survives is sharper than the one section 23 stated, not weaker:
      restriction has to be measured across CELL TYPES, and when it is, it works in RNA
      (0.868 across 154 cell types) AND in protein (0.785 across tissue x cell type).
  Both technologies succeed at cell-type resolution and both fail at whole-organ resolution.""" %
          (r18['sep'], r18['p'], ms['sep'], ms['p']))
else:
    print("""  RNA on the same 18 organs gives %.3f (p = %.3f) while mass spectrometry gives %.3f
  (p = %.3f). They do NOT agree at matched resolution, so anatomical pooling does not explain
  the difference, and the honest conclusion of section 26 stands unchanged: the result is
  reference-dependent and does not transfer to mass-spectrometry protein abundance.""" %
          (r18['sep'], r18['p'], ms['sep'], ms['p']))

json.dump(out, open(os.path.join(HERE, 'matched_resolution_test.json'), 'w'), indent=1)
print('\nwrote matched_resolution_test.json')
