"""Expendability weights from an EXTERNAL measured dataset: cell lifespans (Sender & Milo, 2021).

Source: Sender R, Milo R. "The distribution of cellular turnover in the human body." Nat Med 2021;27:45-48.
Data file: ev_capture_data/cell_turnover_sender_milo/Summary.xlsx  (github.com/milo-lab/cellular_turnover).
Lifespans below are copied from that file (human 'lifespan' in days; where the human value is blank the file's
'lifespan in rodents' value is used, marked R).

IDEA: a cell type that turns over fast is replaced by the body, so losing it is expendable; a long-lived,
non-renewing cell type is not.  w = log10(lifespan_days) / log10(300000), clipped to [0, 1]
(300000 d = neurons, the longest in the file -> 1.0; lifespans <= 1 d -> 0).

The mapping from the 154 HPA cell types to the ~26 categories in the file was written from cell biology
(lineage / tissue family) only; it is fixed here and is not edited after any score is computed.
HPA types with no reasonable category are UNMAPPED and receive a neutral weight (median of the mapped types).
"""
import math

LIFESPAN_DAYS = {  # category: (days, 'H' human value or 'R' rodent value)
    'Erythrocytes': (119, 'H'), 'Neutrophils': (6.6, 'H'), 'Thymocytes': (3.8, 'R'),
    'Mature T cells': (320, 'H'), 'B cells progenitors': (0.83, 'R'), 'Mature B cells': (63, 'H'),
    'Monocytes': (3.5, 'H'), 'Stomach epithelia': (13, 'H'), 'Small intestine epithelia': (5.3, 'H'),
    'Colon epithelia': (3.4, 'H'), 'Epidermal cells': (31, 'H'), 'Dermal fibroblasts': (140, 'R'),
    'Endothelial cells': (2300, 'H'), 'Alveolar epithelial cells': (32, 'R'), 'Alveolar macrophages': (165, 'H'),
    'Bronchial epithelial cells': (139, 'H'), 'Adipocytes': (4000, 'H'), 'Myocytes': (18000, 'H'),
    'Cardiomyocytes': (56000, 'H'), 'Neurons': (300000, 'H'), 'Glial cells': (21000, 'H'),
    'Hepatocytes': (360, 'R'), 'Kupffer cells': (400, 'R'), 'Stellate cells': (30, 'R'),
}
MAX_DAYS = 300000.0

_FAMILIES = {
    'Adipocytes': ['adipocytes'],
    'Alveolar epithelial cells': ['alveolar cells type 1', 'alveolar cells type 2', 'transitional alveolar cells'],
    'Glial cells': ['astrocytes', 'bergmann glia', 'microglia', 'oligodendrocytes', 'oligodendrocyte progenitor cells',
                    'müller glia', 'ependymal cells', 'schwann cells', 'pituicytes/fscs'],
    'Neurons': ['brain excitatory neurons', 'brain inhibitory neurons', 'other brain neurons', 'retinal amacrine cells',
                'retinal bipolar cells', 'retinal ganglion cells', 'retinal horizontal cells', 'rod photoreceptor cells',
                'cone photoreceptor cells'],
    'Mature B cells': ['b-cells', 'plasma cells'],
    'Mature T cells': ['t-cells', 'nk-cells', 'innate lymphoid cells'],
    'Thymocytes': ['thymocytes'],
    'Neutrophils': ['neutrophils'],
    'B cells progenitors': ['neutrophil progenitors', 'monocyte progenitors', 'megakaryocyte progenitors',
                            'megakaryocyte-erythroid progenitors', 'erythrocyte progenitors'],
    'Monocytes': ['monocytes'],
    'Alveolar macrophages': ['macrophages', 'hofbauer cells'],
    'Kupffer cells': ['kupffer cells'],
    'Erythrocytes': ['erythrocytes'],
    'Cardiomyocytes': ['cardiomyocytes'],
    'Myocytes': ['myonuclei', 'myosatellite cells'],
    'Hepatocytes': ['hepatocytes', 'cholangiocytes'],
    'Stellate cells': ['hepatic stellate cells'],
    'Endothelial cells': ['vascular endothelial cells', 'lymphatic endothelial cells'],
    'Dermal fibroblasts': ['fibroblasts', 'fibro-adipogenic progenitors', 'endometrial stromal cells',
                           'decidual stromal cells', 'ovarian stromal cells'],
    'Epidermal cells': ['basal keratinocytes', 'suprabasal keratinocytes', 'esophageal apical cells',
                        'esophageal basal cells', 'esophageal suprabasal cells', 'ocular epithelial cells',
                        'conjunctival goblet cells'],
    'Bronchial epithelial cells': ['respiratory basal cells', 'respiratory ciliated cells',
                                   'respiratory deuterosomal cells', 'respiratory ionocytes',
                                   'respiratory secretory cells', 'submucosal glandular cells'],
    'Stomach epithelia': ['foveolar cells', 'gastric chief cells', 'gastric progenitor cells', 'mucous neck cells',
                          'parietal cells'],
    'Small intestine epithelia': ['enterocytes', 'goblet cells', 'paneth cells', 'tuft cells', 'enteric stem cells',
                                  'enteric transient amplifying cells'],
    'Colon epithelia': ['colonocytes'],
}
MAP = {hpa: cat for cat, lst in _FAMILIES.items() for hpa in lst}
assert set(_FAMILIES) <= set(LIFESPAN_DAYS)


def lifespan(cell_type):
    cat = MAP.get(cell_type)
    return None if cat is None else LIFESPAN_DAYS[cat][0]


def raw_weight(cell_type):
    L = lifespan(cell_type)
    if L is None:
        return None
    return min(1.0, max(0.0, math.log10(L) / math.log10(MAX_DAYS)))


def neutral_weight(all_types):
    """Median weight over the MAPPED types among all_types (used for unmapped types)."""
    ws = sorted(raw_weight(c) for c in all_types if c in MAP)
    n = len(ws)
    return ws[n // 2] if n % 2 else 0.5 * (ws[n // 2 - 1] + ws[n // 2])


def weight_fun(all_types):
    nw = neutral_weight(all_types)
    return lambda c: (nw if raw_weight(c) is None else raw_weight(c))


def flag_fun(threshold_days):
    """1 if the cell type is mapped and its lifespan <= threshold (expendable by turnover), else 0."""
    return lambda c: 1.0 if (lifespan(c) is not None and lifespan(c) <= threshold_days) else 0.0
