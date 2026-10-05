"""
PROTEIN-LEVEL CROSS-CHECK of the whole-body restriction result (results doc section 23).

The stated limitation of section 23 is that the reference is RNA. This script repeats the
identical analysis on Human Protein Atlas PROTEIN references and reports whether the result
transfers. The hypotheses were registered in the pre-registration addendum BEFORE this ran:

  H5  top_share separates safe from toxic with separation > 0.70 and p < 0.05.
      Below that, the RNA result is reported as reference-dependent.
  H6  the two-factor rule beats the one-factor rule. If not, the expendability factor is
      reported as not replicating.

Two protein references are used, and they are NOT equivalent:

  MS   ms_tissue_sample_data.tsv.zip / ms_tissue.tsv.zip -- mass spectrometry intensity.
       QUANTITATIVE, so top_share transfers unchanged. This is the primary protein test.
  IHC  normal_ihc_data.tsv.zip -- antibody staining scored Not detected / Low / Medium / High.
       ORDINAL with four levels, so a "share of total" is a much coarser object here; a gene
       scored High in three tissues and Not detected everywhere else has top_share 1/3 by
       construction. Reported second, with that caveat attached, never as the primary result.

Usage:  python protein_cross_check.py            # runs whichever files are present
        python protein_cross_check.py --headers  # just print each file's columns
"""
import os, io, sys, json, zipfile, collections
import numpy as np
import target_labels as TL
import wholebody_v2 as WB

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'ev_capture_data', 'normal_tissue_atlas'))

IHC_LEVEL = {'not detected': 0.0, 'low': 1.0, 'medium': 2.0, 'high': 3.0}

GENE_COLS = ('Gene name', 'Gene_name', 'Gene')
GROUP_COLS = ('Cell type', 'Tissue', 'Sample', 'Cell_type')
VALUE_COLS = ('nMS', 'Intensity', 'intensity', 'nTPM', 'pTPM', 'nCPM', 'Level',
              'Protein intensity', 'Normalized intensity', 'Value')


def open_first(zf):
    return zf.namelist()[0]


def headers(path):
    z = zipfile.ZipFile(path)
    with z.open(open_first(z)) as f:
        fh = io.TextIOWrapper(f, 'utf-8')
        hdr = fh.readline().rstrip('\n').split('\t')
        sample = fh.readline().rstrip('\n').split('\t')
    return hdr, sample


def load(path, genes, prefer_group=None, composite=None, reliability=None):
    """Return {gene: {group: value}}, plus the column names actually used.

    composite   list of column names to join with ' | ' as the group key. The IHC file is
                indexed by (Tissue, Cell type) and cell-type names repeat across tissues
                ('glandular cells' appears in a dozen organs), so grouping on 'Cell type'
                alone would silently merge unrelated tissues.
    reliability keep only rows whose Reliability is in this set. HPA marks some antibody
                results 'Uncertain'; dropping them is the stricter reading and is run as a
                separate sensitivity analysis, never as a quiet default.
    """
    z = zipfile.ZipFile(path)
    vals = collections.defaultdict(lambda: collections.defaultdict(float))
    counts = collections.defaultdict(lambda: collections.defaultdict(int))
    with z.open(open_first(z)) as f:
        fh = io.TextIOWrapper(f, 'utf-8')
        hdr = fh.readline().rstrip('\n').split('\t')
        gi = next((hdr.index(c) for c in GENE_COLS if c in hdr), None)
        gcands = [c for c in GROUP_COLS if c in hdr]
        if prefer_group and prefer_group in hdr:
            gcands = [prefer_group] + [c for c in gcands if c != prefer_group]
        ci = hdr.index(gcands[0]) if gcands else None
        vi = next((hdr.index(c) for c in VALUE_COLS if c in hdr), None)
        if gi is None or ci is None or vi is None:
            raise SystemExit('cannot identify columns in %s: %s' % (os.path.basename(path), hdr))
        is_level = hdr[vi] == 'Level'
        comp = [hdr.index(c) for c in (composite or []) if c in hdr]
        ri = hdr.index('Reliability') if ('Reliability' in hdr and reliability) else None
        n_seen = n_blank = 0
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) <= max(gi, ci, vi):
                continue
            g = p[gi]
            if g not in genes:
                continue
            n_seen += 1
            if ri is not None and p[ri].strip() not in reliability:
                continue
            key = ' | '.join(p[i] for i in comp) if comp else p[ci]
            raw = p[vi]
            if is_level:
                v = IHC_LEVEL.get(raw.strip().lower())
                if v is None:
                    n_blank += 1
                    continue
            else:
                try:
                    v = float(raw)
                except ValueError:
                    n_blank += 1      # the MS files leave Intensity empty where not measured
                    continue
            # average replicates/samples that share a group label
            vals[g][key] += v
            counts[g][key] += 1
    out = {g: {k: vals[g][k] / counts[g][k] for k in vals[g]} for g in vals}
    gcol = hdr[gi]
    ccol = ' | '.join(composite) if comp else hdr[ci]
    coverage = 1.0 - (n_blank / n_seen if n_seen else 0.0)
    return out, gcol, ccol, hdr[vi], coverage


BLOOD = {'CD19', 'MS4A1', 'CD22', 'CD79B', 'TNFRSF17', 'SLAMF7', 'CD38', 'TNFRSF8',
         'CD33', 'GPRC5D', 'CD52', 'PTPRC', 'IL3RA', 'FCER2', 'CD40', 'CD80',
         'HLA-DRA', 'SDC1', 'CD74', 'CD47'}


def analyse(name, path, label, prefer_group=None, composite=None, reliability=None):
    genes = {t[0] for t in TL.TARGETS}
    raw, gcol, ccol, vcol, coverage = load(path, genes, prefer_group, composite, reliability)
    conf = WB.confinement(raw)
    ngroups = len(next(iter(raw.values()))) if raw else 0
    print('\n' + '=' * 78)
    print('%s   (%s)' % (name, os.path.basename(path)))
    print('=' * 78)
    print('  columns used: gene=%r  group=%r  value=%r' % (gcol, ccol, vcol))
    print('  measured coverage: %.1f%% of rows for these genes carry a value%s'
          % (100 * coverage, '' if coverage > 0.5 else '   <-- SPARSE, read with care'))
    print('  genes resolved: %d of %d ; groups in reference: ~%d'
          % (len(conf), len(genes), ngroups))
    missing = sorted(genes - set(conf))
    if missing:
        print('  NOT FOUND: %s' % ', '.join(missing))

    use = {g: d for g, d in conf.items() if label.get(g) in ('safe', 'toxic')}
    if len({label[g] for g in use}) < 2:
        print('  too few labelled targets resolved to run the test.')
        return None

    res = {}
    s1, p1, ex1, na, nb = WB.run_test({g: d['top_share'] for g, d in use.items()}, label, True)
    print('\n  H5  one factor (top_share)   separation %.3f  %s p = %.5f   (safe n=%d vs toxic n=%d)'
          % (s1, 'exact' if ex1 else 'approx', p1, na, nb))
    print('      registered threshold was separation > 0.70 and p < 0.05  ->  %s'
          % ('MET' if (s1 > 0.70 and p1 < 0.05) else 'NOT MET'))
    sn, pn, _, _, _ = WB.run_test({g: d['n_types'] for g, d in use.items()}, label, False)
    print('      n_types                    separation %.3f  p = %.5f' % (sn, pn))

    sc2 = {g: d['top_share'] + (1.0 if TL.haematopoietic(d['home']) else 0.0)
           for g, d in use.items()}
    n_haem = sum(1 for g, d in use.items() if TL.haematopoietic(d['home']))
    s2, p2, ex2, _, _ = WB.run_test(sc2, label, True)
    print('\n  H6  two factors              separation %.3f  %s p = %.5f   (change %+.3f)'
          % (s2, 'exact' if ex2 else 'approx', p2, s2 - s1))
    print('      the expendability flag fires on %d of %d labelled targets' % (n_haem, len(use)))
    if n_haem == 0:
        print('      NOTE: this reference is indexed by TISSUE, not cell type, so no home can be')
        print('      a blood lineage and the two-factor rule collapses onto the one-factor rule.')
        print('      H6 is NOT TESTABLE on this reference. Reported as such, not as a failure.')

    useS = {g: d for g, d in use.items() if g not in BLOOD}
    if len({label[g] for g in useS}) == 2:
        sS, pS, exS, nsa, nsb = WB.run_test({g: d['top_share'] for g, d in useS.items()}, label, True)
        print('\n      solid tumours only         separation %.3f  %s p = %.5f   (n=%d vs %d)'
              % (sS, 'exact' if exS else 'approx', pS, nsa, nsb))
        res['solid'] = {'sep': sS, 'p': pS, 'n_safe': nsa, 'n_toxic': nsb}

    print('\n  most and least confined in this reference:')
    order = sorted(conf.items(), key=lambda kv: -kv[1]['top_share'])
    for g, d in order[:4] + order[-4:]:
        print('    %-9s %6.3f  %4d groups  home %-28s [%s]'
              % (g, d['top_share'], d['n_types'], d['home'][:28], label.get(g, '-')))

    res.update({'one_factor': s1, 'p_one': p1, 'two_factor': s2, 'p_two': p2,
                'n_groups': ngroups, 'n_haem_homes': n_haem,
                'columns': {'gene': gcol, 'group': ccol, 'value': vcol},
                'per_gene': conf})
    return res


def main():
    label = {t[0]: t[3] for t in TL.TARGETS}
    files = [
        ('PROTEIN, mass spectrometry, 60 samples across 20 tissues (PRIMARY)',
         'ms_tissue_sample_data.tsv.zip', None, None, None),
        ('PROTEIN, mass spectrometry, summarised to 20 tissues',
         'ms_tissue.tsv.zip', None, None, None),
        ('PROTEIN, immunohistochemistry, tissue x cell type (SECONDARY, ordinal)',
         'normal_ihc_data.tsv.zip', None, ['Tissue', 'Cell type'], None),
        ('PROTEIN, immunohistochemistry, uncertain antibodies dropped (sensitivity)',
         'normal_ihc_data.tsv.zip', None, ['Tissue', 'Cell type'],
         {'Enhanced', 'Supported', 'Approved'}),
    ]
    if '--headers' in sys.argv:
        for row in files:
            fn = row[1]
            p = os.path.join(D, fn)
            if os.path.exists(p):
                h, srow = headers(p)
                print('%-34s %s\n%-34s %s\n' % (fn, h, '', srow))
            else:
                print('%-34s (not present)' % fn)
        return

    out = {}
    any_found = False
    for name, fn, prefer, composite, reliability in files:
        p = os.path.join(D, fn)
        if not os.path.exists(p):
            print('\n(skipping %s — not present)' % fn)
            continue
        any_found = True
        r = analyse(name, p, label, prefer, composite, reliability)
        if r:
            out[name] = r
    if not any_found:
        raise SystemExit('No protein reference files found in %s' % D)

    json.dump(out, open(os.path.join(HERE, 'protein_cross_check.json'), 'w'), indent=1)
    print('\nwrote protein_cross_check.json')


if __name__ == '__main__':
    main()
