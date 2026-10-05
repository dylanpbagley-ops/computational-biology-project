"""Coarse lineage families for the 154 HPA single-cell cell types (biology only; no labels, no scores).

Used to define a HOME LINEAGE (the family with the largest aggregate share of a gene's expression) instead of one
arbitrary home cell type. 'blood' is exactly target_labels.haematopoietic (the frozen list); HSC, platelets and
innate lymphoid cells are not on that list, so they sit in 'blood_other' and are NOT counted as blood.
Everything not listed below is 'epithelial' (the biggest natural family: hepatocytes, tubules, keratinocytes,
airway, gut, glands, ducts...).  Fixed before any score using it was computed.
"""
import target_labels as TL

_LISTS = {
    'blood_other': ['hematopoietic stem cells', 'platelets', 'innate lymphoid cells'],
    'germline': ['differentiating spermatogonia', 'undifferentiated spermatogonia', 'early primary spermatocytes',
                 'late primary spermatocytes', 'early spermatids', 'late spermatids', 'oocytes'],
    'placental': ['syncytiotrophoblasts', 'cytotrophoblasts', 'extravillous trophoblasts', 'migrating cytotrophoblasts'],
    'neural': ['astrocytes', 'bergmann glia', 'oligodendrocytes', 'oligodendrocyte progenitor cells', 'müller glia',
               'ependymal cells', 'schwann cells', 'pituicytes/fscs', 'brain excitatory neurons',
               'brain inhibitory neurons', 'other brain neurons', 'retinal amacrine cells', 'retinal bipolar cells',
               'retinal ganglion cells', 'retinal horizontal cells', 'rod photoreceptor cells',
               'cone photoreceptor cells'],
    'endocrine': ['adrenal cortex cells', 'adrenal medulla cells', 'corticotrophs', 'gonadotrophs', 'lactotrophs',
                  'somatotrophs', 'thyrotrophs', 'pituitary stem cells', 'pancreatic islet cells',
                  'neuroendocrine cells', 'sertoli cells', 'leydig cells', 'granulosa cells'],
    'endothelial': ['vascular endothelial cells', 'lymphatic endothelial cells'],
    'stromal_muscle': ['fibroblasts', 'fibro-adipogenic progenitors', 'endometrial stromal cells',
                       'decidual stromal cells', 'ovarian stromal cells', 'smooth muscle cells',
                       'vascular smooth muscle cells', 'pericytes', 'myonuclei', 'myosatellite cells',
                       'cardiomyocytes', 'adipocytes', 'thymic myoid cells', 'peritubular myoid cells',
                       'mesothelial cells', 'epicardial cells', 'hepatic stellate cells'],
    'other': ['melanocytes'],
}
_BY_TYPE = {t: k for k, lst in _LISTS.items() for t in lst}
assert len(_BY_TYPE) == sum(len(v) for v in _LISTS.values()), 'a cell type is listed twice'


def lineage(cell_type):
    if TL.haematopoietic(cell_type):
        return 'blood'
    return _BY_TYPE.get(cell_type, 'epithelial')


def check(all_types):
    """Every listed name must exist in the reference (guards against typos)."""
    missing = [t for t in _BY_TYPE if t not in all_types]
    assert not missing, missing
