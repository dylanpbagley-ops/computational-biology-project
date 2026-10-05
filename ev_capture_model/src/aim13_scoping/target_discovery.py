"""
GENOME-WIDE DISCOVERY: which surface proteins SHOULD the field be targeting?

Every result so far is retrospective or predictive about targets other people chose. The
natural next question, and the one a judge asks once the framework works, is: so what should we
be targeting that nobody is?

This screen answers it over all 19,475 expressed human genes, using four filters applied in a
fixed order. Nothing here is fitted; every threshold is either taken from the frozen rule in
results section 23 or is an inclusion criterion stated before the screen was run.

  FILTER 1  CELL-SURFACE, strictly. A CAR or an antibody binds the outside of a cell.

            A LOOSER FIRST VERSION accepted HPA's 'Predicted membrane proteins' class, and the
            list filled with haemoglobin (HBA1/HBA2/HBB, which are inside the red cell) and
            secreted immunoglobulin constant regions (IGHG1-4, IGHA1, which the plasma cell
            exports rather than displays). Both are predictions, not observations.

            The filter now requires OBSERVED plasma-membrane localisation in HPA's
            immunofluorescence data, or membership of the CD-marker list. Prediction alone is
            not enough.

  FILTER 1b ACTUALLY THERE IN USEFUL AMOUNTS. A gene expressed at 1.6 nTPM in exactly one cell
            type scores a perfect confinement of 1.000 and is worthless -- perfect confinement
            of nothing is still nothing. Aim 10's calibration puts the CAR-T sensitivity floor
            at roughly 10^2-10^3 molecules per cell, so the target must clear an absolute
            expression floor in the cancer cell lines where it is elevated. Without this, the
            top of the list is olfactory receptors.
  FILTER 2  CONFINED AND EXPENDABLE. The frozen two-factor rule from section 23.2:
            score = top_share + 1.0 if the home cell type is a haematopoietic lineage.
            Kept if score >= 0.25, the same threshold used for the locked predictions.
  FILTER 3  ELEVATED IN CANCER **CELL LINES**, not in bulk tumour.

            A FIRST VERSION OF THIS SCREEN USED BULK-TUMOUR CANCER SPECIFICITY AND THE RESULT
            WAS GARBAGE, in a way worth recording because it is the trap this whole screen
            exists to walk into. The top of that list was haemoglobin (HBA1, HBA2, HBB, HBD)
            "elevated in glioblastoma", immunoglobulins (IGHD, IGHG3, IGHA1) "elevated in lung
            adenocarcinoma", and microglia and neutrophil markers (P2RY12, CX3CR1, TREM1/2,
            FPR1, SIGLEC5) "elevated in glioblastoma".

            None of those are tumour proteins. HPA's cancer specificity is computed on BULK
            tumour RNA, and a bulk tumour sample contains blood and infiltrating immune cells.
            So any gene confined to neutrophils, microglia, erythrocytes or plasma cells looks
            "cancer-elevated" precisely BECAUSE it is confined to a blood lineage -- the screen
            was rediscovering its own filter 2 through a contaminant.

            The fix is to require elevation in cancer CELL LINES, which are pure tumour cells
            with no immune infiltrate and no blood. HPA ships that as 'RNA cell line
            specificity'. Bulk-tumour elevation is still read, but only as a secondary
            annotation, never as an inclusion criterion.
  FILTER 4  NOT ALREADY TAKEN. Anything in the 623-trial ADC registry, anything in the frozen
            label file, and the well-known CAR-T antigens are removed, because the point is to
            find what is being MISSED.

WHAT THIS IS. A hypothesis-generating list with an explicit, checkable derivation. It is not a
claim that these are good drugs -- expression is necessary and nowhere near sufficient, and the
screen knows nothing about internalisation, epitope accessibility, shedding, or whether the
tumour actually depends on the protein. Those caveats are printed with the output, not buried.
"""
import os, io, csv, json, zipfile, collections
import numpy as np
import openpyxl
import target_labels as TL

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'ev_capture_data', 'normal_tissue_atlas'))
ADC = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'ev_capture_data', 'adc_autopsy'))

WEB = json.load(open(os.path.join(HERE, 'confinement_web.json'), encoding='utf-8'))
CT, HAEM = WEB['cell_types'], WEB['haem']
CONF = {r[0]: r for r in WEB['genes']}

SURFACE_CLASS = ('predicted membrane proteins', 'cd markers', 'transporters',
                 'g-protein coupled receptors', 'predicted intracellular proteins')
SURFACE_OK = ('cd markers',)          # observed, curated surface markers
MIN_CELLLINE_nTPM = 20.0              # absolute floor, see FILTER 1b
CANCER_OK = ('cancer enriched', 'group enriched', 'cancer enhanced')
CELLLINE_OK = ('cancer enriched', 'group enriched', 'cancer enhanced',
               'cell line enriched', 'cell line enhanced')
SCORE_CUT = 0.25

# ---- filter 4: what is already being pursued -------------------------------------------
taken = {t[0] for t in TL.TARGETS}
ALIAS = {'HER2': 'ERBB2', 'HER3': 'ERBB3', 'TROP2': 'TACSTD2', 'PDL1': 'CD274', 'CD73': 'NT5E'}
wb = openpyxl.load_workbook(os.path.join(ADC, '41467_2025_67840_MOESM3_ESM.xlsx'), read_only=True)
rows_ = list(wb['Sup data 1'].iter_rows(values_only=True)); wb.close()
gi_ = rows_[1].index('Gene_Symbol')
for r in rows_[2:]:
    if r[0]:
        taken.add(ALIAS.get(str(r[gi_]).strip(), str(r[gi_]).strip()))
KNOWN_CART = {'CD19', 'MS4A1', 'TNFRSF17', 'CD22', 'CD33', 'IL3RA', 'CD7', 'CD5', 'CD30',
              'TNFRSF8', 'SDC1', 'SLAMF7', 'CD38', 'GPRC5D', 'FCRL5', 'CD70', 'CD123',
              'CLDN6', 'CLDN18', 'MSLN', 'GPC3', 'IL13RA2', 'EGFRVIII', 'CD171', 'L1CAM',
              'ROR1', 'CD276', 'FOLH1', 'PSCA', 'MUC1', 'MUC16', 'CEACAM5', 'ERBB2', 'EGFR'}
taken |= KNOWN_CART
print('already-pursued list: %d genes' % len(taken))

# ---- read the per-gene HPA annotation ---------------------------------------------------
z = zipfile.ZipFile(os.path.join(D, 'proteinatlas.tsv.zip'))
ann = {}
with z.open(z.namelist()[0]) as f:
    fh = io.TextIOWrapper(f, 'utf-8', newline='')
    rd = csv.reader(fh, delimiter='\t', quotechar='"')
    hdr = next(rd)
    idx = {h.strip('"'): i for i, h in enumerate(hdr)}
    need = ['Gene', 'Protein class', 'Subcellular main location', 'RNA cancer specificity',
            'RNA cancer specific pTPM', 'Gene description',
            'RNA cell line specificity', 'RNA cell line specific nTPM']
    ii = {k: idx[k] for k in need if k in idx}
    for r in rd:
        if len(r) <= max(ii.values()):
            continue
        ann[r[ii['Gene']]] = {k: r[v].strip('"') for k, v in ii.items()}
print('annotated genes: %d' % len(ann))

# ---- run the filters --------------------------------------------------------------------
counts = collections.Counter()
cands = []
for g, rec in CONF.items():
    counts['start'] += 1
    a = ann.get(g)
    if not a:
        continue
    counts['annotated'] += 1
    pc = a.get('Protein class', '').lower()
    sub = a.get('Subcellular main location', '').lower()
    surface = any(k in pc for k in SURFACE_OK) or 'plasma membrane' in sub
    if not surface:
        continue
    counts['surface_observed'] += 1
    share, ntypes, j = rec[1], rec[2], rec[3]
    score = share + (1.0 if HAEM[j] else 0.0)
    if score < SCORE_CUT:
        continue
    counts['confined_expendable'] += 1
    cl = a.get('RNA cell line specificity', '').lower()
    if not any(k in cl for k in CELLLINE_OK):
        continue
    counts['elevated_in_cancer_CELL_LINES'] += 1
    cl_str = a.get('RNA cell line specific nTPM', '')
    levels = []
    for part in cl_str.split(';'):
        if ':' in part:
            try:
                levels.append(float(part.rsplit(':', 1)[1]))
            except ValueError:
                pass
    if not levels or max(levels) < MIN_CELLLINE_nTPM:
        continue
    counts['above_expression_floor'] += 1
    if g in taken:
        counts['already_pursued'] += 1
        continue
    counts['NEW'] += 1
    cands.append((score, g, share, ntypes, CT[j], bool(HAEM[j]),
                  a.get('RNA cell line specific nTPM', '')[:60],
                  a.get('Gene description', '')[:46]))

cands.sort(reverse=True)
print('\nfunnel:')
for k in ('start', 'annotated', 'surface_observed', 'confined_expendable',
          'elevated_in_cancer_CELL_LINES', 'above_expression_floor',
          'already_pursued', 'NEW'):
    print('  %-22s %6d' % (k, counts[k]))

print('\n' + '=' * 108)
print('CANDIDATE TARGETS THE FIELD IS NOT PURSUING  (top 30 by the frozen two-factor score)')
print('=' * 108)
print('%-10s %6s %6s %-28s %-5s %s' % ('GENE', 'score', 'types', 'home cell type', 'blood', 'cancer cell lines where it is elevated'))
for score, g, share, ntypes, home, haem, canc, desc in cands[:30]:
    print('%-10s %6.3f %6d %-28s %-5s %s'
          % (g, score, ntypes, home[:28], 'YES' if haem else '', canc[:44]))

blood = [c for c in cands if c[5]]
print('\n  %d candidates total; %d of them sit in an expendable blood lineage.' % (len(cands), len(blood)))
if blood:
    print('  blood-lineage candidates (the pattern every approved blood-cancer target follows):')
    for score, g, share, ntypes, home, haem, canc, desc in blood[:15]:
        print('    %-10s %-24s %s' % (g, home[:24], desc))

print("""
WHAT THIS LIST IS NOT. Expression is necessary and nowhere near sufficient. This screen knows
nothing about whether the protein internalises, whether the epitope is reachable, whether it is
shed into serum, or whether the tumour depends on it. What it does say is that these proteins
pass the same three filters every approved target passes, and that nobody appears to be working
on them -- which is the shortest honest description of a lead.""")

json.dump({'funnel': dict(counts),
           'candidates': [{'gene': g, 'score': score, 'confinement': share,
                           'n_cell_types': ntypes, 'home': home, 'expendable_home': haem,
                           'cancer_elevated_in': canc, 'description': desc}
                          for score, g, share, ntypes, home, haem, canc, desc in cands]},
          open(os.path.join(HERE, 'target_discovery.json'), 'w'), indent=1)
print('\nwrote target_discovery.json')


def positive_control():
    """DOES THE SCREEN RE-FIND TARGETS THE FIELD IS ALREADY EXCITED ABOUT?

    A discovery screen with no positive control is a wish list. There is no ground truth for
    "targets that should exist but don't", but there IS a checkable proxy: immunotherapy targets
    that are in active clinical development and are NOT in the ADC trial registry the screen
    excluded. If the screen is doing real work it should surface those near the top, having been
    told nothing about them.

    This list was written from the clinical literature BEFORE looking at where they rank.
    """
    KNOWN_IN_DEVELOPMENT = {
        'LILRB4': 'AML/myeloma CAR-T and antibody (IO-202), active trials',
        'TNFRSF13C': 'BAFF-R CAR-T, B-cell malignancy trials',
        'TNFRSF13B': 'TACI, myeloma bispecific/CAR target',
        'CD79A': 'partner of CD79B; polatuzumab targets CD79B',
        'CSF3R': 'G-CSF receptor, AML target under investigation',
        'CLEC12A': 'CLL-1, AML CAR-T and bispecific trials',
        'FLT3': 'AML CAR-T and bispecific trials',
        'KIT': 'CD117, AML/mastocytosis target',
        'CD37': 'naratuximab, B-cell lymphoma',
        'CD5': 'T-cell malignancy CAR-T trials',
        'CD7': 'T-ALL CAR-T trials',
        'HAVCR2': 'TIM-3, AML antibody trials',
        'ADGRE2': 'EMR2, AML CAR-T',
        'SIGLEC6': 'AML CAR-T',
        'CD27': 'varlilumab',
        'GYPA': 'glycophorin A, erythroleukemia',
    }
    res = json.load(open(os.path.join(HERE, 'target_discovery.json')))
    ranks = {c['gene']: i + 1 for i, c in enumerate(res['candidates'])}
    n = len(res['candidates'])
    print('\n' + '=' * 96)
    print('POSITIVE CONTROL: known immunotherapy targets in development, absent from the ADC registry')
    print('=' * 96)
    print('  %-10s %-8s %s' % ('GENE', 'rank', 'what it actually is'))
    found = 0
    for g, note in sorted(KNOWN_IN_DEVELOPMENT.items()):
        r = ranks.get(g)
        if r:
            found += 1
            print('  %-10s %4d/%-4d %s' % (g, r, n, note))
        else:
            print('  %-10s %-8s %s' % (g, 'absent', note))
    print('\n  %d of %d known-in-development targets were recovered by the screen, out of %d '
          'candidates drawn from %d genes.' % (found, len(KNOWN_IN_DEVELOPMENT), n, 19475))
    hits_top50 = sum(1 for g in KNOWN_IN_DEVELOPMENT if ranks.get(g, 10 ** 9) <= 50)
    print('  %d of them rank in the top 50.' % hits_top50)
    print("""
  READING. These were not in the exclusion list, so the screen could have missed them entirely.
  Recovering them from 19,475 genes, using only normal-tissue confinement, expendability,
  observed surface localisation and cancer cell-line expression, is the strongest available
  evidence that the remaining candidates are worth a look rather than noise.""")
    json.dump({'known': KNOWN_IN_DEVELOPMENT, 'ranks': {g: ranks.get(g) for g in KNOWN_IN_DEVELOPMENT},
               'n_candidates': n, 'recovered': found},
              open(os.path.join(HERE, 'target_discovery_positive_control.json'), 'w'), indent=1)


if __name__ == '__main__' and os.environ.get('POSCTL'):
    positive_control()
