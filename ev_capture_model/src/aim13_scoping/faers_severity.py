"""Measured organ-system severity from FDA adverse-event reports (openFDA / FAERS), used as expendability weights.

WHY: hand-assigned and lifespan-based weights are judgement or a proxy. This module defines a MEASURED weight: for each
organ system, among all FAERS reports that record a loss-of-function reaction of that system, what fraction record
the patient's death?  Severity = P(death | report of that organ-system reaction). It builds in real-world supportive
care (dialysis, insulin, transfusion, G-CSF, immunoglobulin) because it is measured after the fact.
Fixed BEFORE any count was fetched and before any target score was computed (see faers_precommitment.md).

CAVEATS (stated up front): FAERS is spontaneous reporting; fatal reports are over-reported; patients are often
seriously ill from their underlying disease; a PT-to-system assignment is a human choice. The weight is a
rank-informative estimate of severity, not a probability of dying from losing the organ.
"""
import json, os, re

HERE = os.path.dirname(os.path.abspath(__file__))
COUNTS = os.path.join(HERE, 'faers_counts.json')

SYSTEM_PTS = {
    'liver': ['HEPATIC FAILURE', 'ACUTE HEPATIC FAILURE', 'HEPATOTOXICITY', 'DRUG-INDUCED LIVER INJURY', 'HEPATIC NECROSIS'],
    'biliary': ['CHOLESTASIS', 'CHOLANGITIS', 'CHOLESTATIC LIVER INJURY'],
    'kidney_tubule': ['ACUTE KIDNEY INJURY', 'RENAL FAILURE', 'RENAL TUBULAR NECROSIS', 'RENAL IMPAIRMENT'],
    'kidney_glomerulus': ['NEPHROTIC SYNDROME', 'GLOMERULONEPHRITIS', 'PROTEINURIA'],
    'heart_muscle': ['CARDIAC FAILURE', 'MYOCARDITIS', 'CARDIOMYOPATHY', 'CARDIOTOXICITY', 'LEFT VENTRICULAR DYSFUNCTION'],
    'lung': ['PNEUMONITIS', 'INTERSTITIAL LUNG DISEASE', 'ACUTE RESPIRATORY DISTRESS SYNDROME', 'RESPIRATORY FAILURE', 'PULMONARY FIBROSIS'],
    'brain': ['ENCEPHALOPATHY', 'TOXIC ENCEPHALOPATHY', 'NEUROTOXICITY', 'LEUKOENCEPHALOPATHY'],
    'peripheral_nerve': ['PERIPHERAL NEUROPATHY', 'POLYNEUROPATHY', 'NEUROPATHY PERIPHERAL', 'GUILLAIN-BARRE SYNDROME'],
    'colon_gut': ['COLITIS', 'ENTEROCOLITIS', 'INTESTINAL PERFORATION', 'GASTROINTESTINAL PERFORATION'],
    'upper_gi': ['OESOPHAGITIS', 'GASTRITIS', 'GASTRIC ULCER', 'STOMATITIS'],
    'pancreas_exocrine': ['PANCREATITIS', 'PANCREATITIS ACUTE'],
    'skin': ['TOXIC EPIDERMAL NECROLYSIS', 'STEVENS-JOHNSON SYNDROME', 'DERMATITIS EXFOLIATIVE',
             'DRUG REACTION WITH EOSINOPHILIA AND SYSTEMIC SYMPTOMS'],
    'melanocyte': ['VITILIGO'],
    'neutrophil': ['NEUTROPENIA', 'FEBRILE NEUTROPENIA', 'AGRANULOCYTOSIS'],
    'red_cell': ['ANAEMIA', 'APLASIA PURE RED CELL'],
    'platelet': ['THROMBOCYTOPENIA', 'IMMUNE THROMBOCYTOPENIA'],
    'marrow_stem': ['APLASTIC ANAEMIA', 'PANCYTOPENIA', 'BONE MARROW FAILURE', 'MYELOSUPPRESSION'],
    'lymphoid': ['HYPOGAMMAGLOBULINAEMIA', 'B-CELL APLASIA', 'LYMPHOPENIA', 'LYMPHOCYTE COUNT DECREASED'],
    'thyroid': ['HYPOTHYROIDISM', 'THYROIDITIS'],
    'adrenal': ['ADRENAL INSUFFICIENCY', 'ADRENOCORTICAL INSUFFICIENCY ACUTE'],
    'islet': ['TYPE 1 DIABETES MELLITUS', 'DIABETIC KETOACIDOSIS'],
    'pituitary': ['HYPOPITUITARISM', 'HYPOPHYSITIS'],
    'eye_surface': ['KERATITIS', 'DRY EYE', 'ULCERATIVE KERATITIS', 'CORNEAL DISORDER'],
    'retina_nerve': ['RETINOPATHY', 'OPTIC NEURITIS', 'VISUAL IMPAIRMENT', 'BLINDNESS'],
    'reproductive': ['INFERTILITY', 'AZOOSPERMIA', 'PREMATURE MENOPAUSE', 'OVARIAN FAILURE'],
    'salivary': ['DRY MOUTH', 'SIALOADENITIS'],
    'vascular': ['CAPILLARY LEAK SYNDROME', 'VASCULITIS', 'THROMBOTIC MICROANGIOPATHY'],
    'skeletal_muscle': ['RHABDOMYOLYSIS', 'MYOSITIS', 'MYOPATHY'],
}

_C = {
    'lymphoid': ['b-cells', 'plasma cells', 't-cells', 'thymocytes', 'nk-cells', 'innate lymphoid cells', 'cdc', 'pdcs',
                 'macrophages', 'monocytes', 'mast cells', 'kupffer cells', 'hofbauer cells', 'microglia',
                 'medullary thymic epithelial cells', 'thymic myoid cells'],
    'neutrophil': ['neutrophils', 'neutrophil progenitors'],
    'red_cell': ['erythrocytes', 'erythrocyte progenitors'],
    'platelet': ['platelets', 'megakaryocytes', 'megakaryocyte progenitors', 'megakaryocyte-erythroid progenitors'],
    'marrow_stem': ['hematopoietic stem cells', 'monocyte progenitors'],
    'vascular': ['vascular endothelial cells', 'lymphatic endothelial cells', 'pericytes', 'vascular smooth muscle cells'],
    'adrenal': ['adrenal cortex cells', 'adrenal medulla cells'],
    'lung': ['alveolar cells type 1', 'alveolar cells type 2', 'transitional alveolar cells', 'respiratory basal cells',
             'respiratory ciliated cells', 'respiratory deuterosomal cells', 'respiratory ionocytes',
             'respiratory secretory cells', 'submucosal glandular cells'],
    'brain': ['astrocytes', 'bergmann glia', 'brain excitatory neurons', 'brain inhibitory neurons', 'other brain neurons',
              'ependymal cells', 'choroid plexus epithelial cells', 'oligodendrocytes', 'oligodendrocyte progenitor cells'],
    'peripheral_nerve': ['schwann cells'],
    'retina_nerve': ['cone photoreceptor cells', 'rod photoreceptor cells', 'retinal amacrine cells', 'retinal bipolar cells',
                     'retinal ganglion cells', 'retinal horizontal cells', 'müller glia', 'retinal pigment epithelial cells'],
    'eye_surface': ['conjunctival goblet cells', 'lacrimal acinar cells', 'ocular epithelial cells'],
    'skin': ['basal keratinocytes', 'suprabasal keratinocytes'],
    'melanocyte': ['melanocytes'],
    'heart_muscle': ['cardiomyocytes', 'epicardial cells'],
    'skeletal_muscle': ['myonuclei', 'myosatellite cells'],
    'liver': ['hepatocytes', 'hepatic stellate cells'],
    'biliary': ['cholangiocytes'],
    'colon_gut': ['colonocytes', 'enterocytes', 'enteric stem cells', 'enteric transient amplifying cells', 'goblet cells',
                  'paneth cells', 'tuft cells'],
    'upper_gi': ['foveolar cells', 'gastric chief cells', 'gastric progenitor cells', 'mucous neck cells', 'parietal cells',
                 'esophageal apical cells', 'esophageal basal cells', 'esophageal suprabasal cells'],
    'pancreas_exocrine': ['pancreatic acinar cells', 'pancreatic duct cells'],
    'islet': ['pancreatic islet cells'],
    'pituitary': ['corticotrophs', 'gonadotrophs', 'lactotrophs', 'somatotrophs', 'thyrotrophs', 'pituitary stem cells', 'pituicytes/fscs'],
    'kidney_tubule': ['distal convoluted tubule cells', 'loop of henle epithelial cells', 'proximal tubule cells',
                      'renal collecting duct intercalated cells', 'renal collecting duct principal cells',
                      'renal connecting tubule cells', 'papillary tip epithelial cells'],
    'kidney_glomerulus': ['podocytes'],
    'salivary': ['salivary acinar cells', 'salivary basal cells', 'salivary duct cells', 'salivary ionocytes', 'salivary myoepithelial cells'],
    'reproductive': ['basal prostatic cells', 'prostatic club cells', 'prostatic glandular cells', 'prostatic hillock cells',
                     'breast hormone-responsive cells', 'breast lactating cells', 'breast myoepithelial cells', 'breast secretory cells',
                     'cytotrophoblasts', 'extravillous trophoblasts', 'migrating cytotrophoblasts', 'syncytiotrophoblasts',
                     'decidual stromal cells', 'differentiating spermatogonia', 'undifferentiated spermatogonia',
                     'early primary spermatocytes', 'late primary spermatocytes', 'early spermatids', 'late spermatids', 'oocytes',
                     'sertoli cells', 'leydig cells', 'granulosa cells', 'ovarian stromal cells', 'peritubular myoid cells',
                     'epididymal basal cells', 'epididymal clear cells', 'epididymal efferent duct absorptive cells',
                     'epididymal efferent duct ciliated cells', 'epididymal principal cells',
                     'endometrial ciliated cells', 'endometrial glandular cells', 'endometrial luminal cells',
                     'endometrial secretory cells', 'endometrial stromal cells', 'fallopian secretory cells', 'fallopian tube ciliated cells'],
}
# No thyroid follicular cell exists among the 154 types. UNMAPPED types get the median measured weight (sensitivity: min and max).
UNMAPPED_TYPES = ['adipocytes', 'fibroblasts', 'fibro-adipogenic progenitors', 'mesothelial cells', 'smooth muscle cells',
                  'neuroendocrine cells', 'urothelial cells']
CELL_SYSTEM = {c: s for s, lst in _C.items() for c in lst}
assert len(CELL_SYSTEM) == sum(len(v) for v in _C.values()), 'a cell type is listed twice'
assert not set(UNMAPPED_TYPES) & set(CELL_SYSTEM)


def check_cell_types(all_types):
    """Every one of the reference's cell types must be mapped or explicitly unmapped, and every name must exist."""
    missing = [t for t in all_types if t not in CELL_SYSTEM and t not in UNMAPPED_TYPES]
    bad = [t for t in list(CELL_SYSTEM) + UNMAPPED_TYPES if t not in all_types]
    assert not missing and not bad, (missing, bad)


def ihc_system(tissue, cell_type):
    """Map an HPA IHC row (tissue, cell type) to a system by identity of the tissue/cell type; None = unmapped."""
    t, c = tissue.lower(), cell_type.lower()
    if 'endothelial' in c or 'vascular mural' in c:
        return 'vascular'
    lymph_c = ('germinal', 'lymphoid tissue', 'lymphocytes', 'langerhans', 'macrophages', 'cells in red pulp', 'cells in white pulp',
               'cortical cells', 'medullary cells')
    if any(k in c for k in lymph_c) and t in ('appendix', 'lymph node', 'spleen', 'thymus', 'tonsil', 'skin', 'endometrium', 'lung', 'colon', 'rectum'):
        return 'lymphoid'
    if 'peripheral nerve' in c:
        return 'peripheral_nerve'
    if t == 'liver':
        return 'biliary' if 'cholangio' in c else 'liver'
    if t == 'gallbladder':
        return 'biliary'
    if t == 'kidney':
        return 'kidney_glomerulus' if ('glomerul' in c or 'bowman' in c) else 'kidney_tubule'
    if t == 'pancreas':
        return 'islet' if 'endocrine' in c else 'pancreas_exocrine'
    if t == 'skin' or t in ('hair', 'sole of foot'):
        if 'melanocyte' in c:
            return 'melanocyte'
        return None if any(k in c for k in ('fibro', 'extracellular')) else 'skin'
    simple = {'heart muscle': 'heart_muscle', 'skeletal muscle': 'skeletal_muscle', 'lung': 'lung', 'bronchus': 'lung',
              'caudate': 'brain', 'cerebellum': 'brain', 'cerebral cortex': 'brain', 'dorsal raphe': 'brain', 'hippocampus': 'brain',
              'hypothalamus': 'brain', 'substantia nigra': 'brain', 'choroid plexus': 'brain',
              'colon': 'colon_gut', 'rectum': 'colon_gut', 'duodenum': 'colon_gut', 'small intestine': 'colon_gut', 'appendix': 'colon_gut',
              'esophagus': 'upper_gi', 'oral mucosa': 'upper_gi', 'stomach': 'upper_gi', 'tonsil': 'upper_gi',
              'eye': 'eye_surface', 'retina': 'retina_nerve', 'bone marrow': 'marrow_stem',
              'lymph node': 'lymphoid', 'spleen': 'lymphoid', 'thymus': 'lymphoid',
              'adrenal gland': 'adrenal', 'thyroid gland': 'thyroid', 'pituitary gland': 'pituitary', 'salivary gland': 'salivary'}
    if t in simple:
        return simple[t]
    if t in ('breast', 'lactating breast', 'cervix', 'endometrium', 'epididymis', 'efferent ducts', 'fallopian tube', 'ovary', 'placenta',
             'prostate', 'seminal vesicle', 'testis', 'vagina'):
        return 'reproductive'
    return None    # adipose, brown adipose, cartilage, soft tissue, smooth muscle, urinary bladder, parathyroid, nasopharynx


def severity_table():
    """{system: {'reports': n, 'deaths': d, 'death_frac': d/n}} from faers_counts.json (written by the fetch step)."""
    c = json.load(open(COUNTS))
    out = {}
    for s in SYSTEM_PTS:
        n, d = c['systems'][s]['reports'], c['systems'][s]['deaths']
        out[s] = {'reports': n, 'deaths': d, 'death_frac': d / n}
    return out


def measured_weights():
    tab = severity_table()
    mx = max(v['death_frac'] for v in tab.values())
    return {s: v['death_frac'] / mx for s, v in tab.items()}
