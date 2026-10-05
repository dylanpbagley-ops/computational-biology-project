"""The set of antigens already under therapeutic development -- the screen's exclusion set.

WHY THIS EXISTS
---------------
`target_discovery.py` originally excluded only the genes in the clinical label panel
(`taken = {t[0] for t in TL.TARGETS}`, 44 genes). That meant the screen had NO knowledge of
any target in development that did not happen to be in the panel, so its "214 NEW
candidates" overstated novelty by an unknown amount.

This was found concretely while building the DCT dossier: TYRP1 is in active CAR-T
development (targeting its surface expression in melanoma), TYRP1 is not in the label panel,
and nothing in the pipeline would have flagged it had it reached the candidate list.

WHAT IS IN HERE
---------------
Three tiers, kept separate so a reader can see where each exclusion comes from and disagree
with any tier independently.

  PANEL      the 44 genes in target_labels.py -- outcomes already adjudicated.
  REGISTRY   the 72 ADC targets compiled by De Schepper et al. 2025 (Nat Commun), taken
             from that paper's own supplementary data, already local in
             ev_capture_data/adc_autopsy/. Compiled by an independent group, which is why
             it is the strongest tier: it is not this project's own judgement.
  CURATED    additional antigens under CAR-T / bispecific / TCR development that appear in
             neither of the above. Each carries a one-line reason. This tier IS this
             project's own judgement and is the one to challenge first.

HONEST LIMITS
-------------
This is still not exhaustive. There is no single public registry of every antigen under
development across every modality, and the CURATED tier is assembled by hand. The correct
claim after this change is "not in the label panel, the De Schepper ADC registry, or a
curated list of known CAR-T/bispecific targets" -- NOT "novel". A future session with access
to a trials database could widen this further, and should.
"""
import os

import target_labels as TL

# ---- tier 2: De Schepper et al. 2025 ADC trial registry ------------------------------
# Extracted from 41467_2025_67840_MOESM13_ESM.xlsx, sheet "Fig 6B - sup data 7", column
# 'gene' (72 entries). Aliases in that sheet are mapped to HGNC symbols below.
REGISTRY_RAW = [
    'ADAM9', 'ALCAM', 'AXL', 'B3GALT5', 'CA9', 'CD276', 'CD38', 'CD44', 'CD46', 'CD70',
    'CD73', 'CDH3', 'CDH6', 'CEACAM5', 'CEACAM6', 'CFC1B', 'CLDN18', 'CLDN6', 'DLK1',
    'DLL3', 'DPEP3', 'EDNRB', 'EFNA4', 'EGFR', 'ENPP3', 'EPHA2', 'EPHA5', 'F3', 'FAP',
    'FGFR2', 'FGFR3', 'FN1', 'FOLH1', 'FOLR1', 'FUT2', 'GPNMB', 'GUCY2C', 'HAVCR1', 'HER2',
    'HER3', 'IGF1R', 'ITGAV', 'ITGB6', 'KIT', 'LAMP1', 'LRRC15', 'LY6E', 'LY75', 'LYPD3',
    'MET', 'MSLN', 'MUC1', 'MUC16', 'NCAM1', 'NECTIN4', 'NOTCH3', 'PDL1', 'PTK7', 'ROR1',
    'ROR2', 'SEZ6', 'SLC34A2', 'SLC39A6', 'SLC44A4', 'SLITRK6', 'ST8SIA1', 'STEAP1', 'TFRC',
    'TNFRSF8', 'TPBG', 'TROP2', 'VTCN1',
]

# the registry uses common names for four entries; map to the symbols the HPA data uses
ALIASES = {
    'HER2': 'ERBB2',
    'HER3': 'ERBB3',
    'PDL1': 'CD274',
    'TROP2': 'TACSTD2',
    'CD73': 'NT5E',
}

# ---- tier 3: curated, and therefore the tier to challenge ----------------------------
CURATED = {
    'TYRP1': 'CAR-T targeting surface TYRP1 in cutaneous and rare melanoma (CIRM-funded)',
    'CD19': 'foundational CAR-T target',
    'CD20': 'rituximab and CD20-directed bispecifics',
    'MS4A1': 'CD20, HGNC symbol',
    'CD22': 'inotuzumab; CD22 CAR-T',
    'TNFRSF17': 'BCMA -- ide-cel, cilta-cel, teclistamab',
    'GPRC5D': 'talquetamab',
    'FCRL5': 'FcRH5 -- cevostamab, in trials',
    'SLAMF7': 'elotuzumab',
    'CD33': 'gemtuzumab; CD33 CAR-T in AML',
    'IL3RA': 'CD123 -- tagraxofusp; CD123 CAR-T',
    'CLEC12A': 'CLL-1 -- AML CAR-T programmes',
    'LILRB4': 'IO-202, AML/myeloma',
    'TNFRSF13C': 'BAFF-R CAR-T, B-cell malignancy trials',
    'TNFRSF13B': 'TACI -- myeloma bispecific/CAR target',
    'CD27': 'varlilumab',
    'CD79A': 'partner of the polatuzumab target; CAR programmes',
    'CD5': 'CD5 CAR-T, T-cell malignancy',
    'CD7': 'CD7 CAR-T, T-cell malignancy',
    'FLT3': 'FLT3 CAR-T / bispecifics in AML',
    'CD30': 'brentuximab; CD30 CAR-T',
    'TNFRSF8': 'CD30, HGNC symbol',
    'GPC3': 'GPC3 CAR-T in hepatocellular carcinoma',
    'CD47': 'magrolimab and successors',
    'EPCAM': 'catumaxomab; EpCAM CAR-T',
    'PSCA': 'PSCA CAR-T in prostate/pancreatic',
    'CEACAM5': 'CEA CAR-T and ADCs',
    'MAGEA4': 'MAGE-A4 TCR therapy (afamitresgene)',
    'NY-ESO-1': 'NY-ESO-1 TCR therapy',
    'CTAG1B': 'NY-ESO-1, HGNC symbol',
    'KLK2': 'KLK2-directed programmes in prostate',
    'DLL3': 'tarlatamab; DLL3 CAR-T',
    'SLC34A2': 'NaPi2b ADCs',
    'B7H3': 'B7-H3 CAR-T (CD276)',
    'CD276': 'B7-H3, HGNC symbol',
}


def registry_symbols():
    """The De Schepper registry, with aliases resolved to HGNC symbols."""
    return {ALIASES.get(g, g) for g in REGISTRY_RAW}


def panel_symbols():
    return {t[0] for t in TL.TARGETS}


def pursued(include_curated=True):
    """Every gene the screen should NOT call novel."""
    s = panel_symbols() | registry_symbols()
    if include_curated:
        s |= set(CURATED)
    return s


def tiers():
    p, r, c = panel_symbols(), registry_symbols(), set(CURATED)
    return {
        'panel': sorted(p),
        'registry': sorted(r),
        'curated': sorted(c),
        'n_panel': len(p), 'n_registry': len(r), 'n_curated': len(c),
        'n_union': len(p | r | c),
        'registry_beyond_panel': sorted(r - p),
        'curated_beyond_panel_and_registry': sorted(c - p - r),
    }


if __name__ == '__main__':
    t = tiers()
    print('exclusion set tiers')
    print('  panel    %3d' % t['n_panel'])
    print('  registry %3d  (+%d beyond the panel)' % (t['n_registry'], len(t['registry_beyond_panel'])))
    print('  curated  %3d  (+%d beyond panel and registry)' % (t['n_curated'], len(t['curated_beyond_panel_and_registry'])))
    print('  UNION    %3d' % t['n_union'])
    print()
    print('registry adds:', ', '.join(t['registry_beyond_panel']))
    print()
    print('curated adds:', ', '.join(t['curated_beyond_panel_and_registry']))
