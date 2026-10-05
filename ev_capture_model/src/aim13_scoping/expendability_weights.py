"""Continuous expendability weights, replacing the binary haematopoietic flag.

WHY THIS EXISTS
---------------
Section 29 proves the ceiling theorem for a CONTINUOUS weight w_c on each normal cell type
(0 = the body replaces it or can live without it, 1 = loss is catastrophic and permanent).
The project then uses a BINARY approximation: +1.0 if the home cell type is haematopoietic,
0 otherwise. Two costs follow.

1. The project proves something more general than it uses. The generality is decorative.
2. The binary flag is substantively wrong in both directions, and demonstrably so:
     * It treats T cells, microglia and B cells as equally expendable. T-cell aplasia is
       severe immunodeficiency; microglia are CNS-resident and do not repopulate from
       marrow; B-cell aplasia is managed routinely with immunoglobulin replacement.
     * It treats gastric mucosa, melanocytes and cardiomyocytes as equally essential.
       Claudin-18.2 (gastric) is an approved target; melanocyte-directed therapy causes
       vitiligo and is tolerated; cardiomyocyte loss is irreversible.
   In the discovery screen this inverts the ranking: AKAP10 (expressed across 137 cell
   types, confinement 0.026) outranks DCT (confinement 0.948, one cell type) purely because
   AKAP10's top cell type happens to be a neutrophil.

WHAT A WEIGHT MEANS HERE
------------------------
w_c is the fraction of harm that counts. Three stated criteria, applied in this order:
  (a) REGENERATION  -- does the tissue continuously renew from a stem/progenitor pool?
  (b) REPLACEABILITY -- is there an established clinical replacement for its function?
  (c) SURVIVABILITY  -- is acute loss survivable while regeneration happens?
Criterion (c) is why colonocytes are NOT scored as freely expendable despite renewing
rapidly: barrier function fails during the gap, which is exactly the CEACAM5 colitis case.

HONESTY ABOUT WHAT THIS IS
--------------------------
These weights are assignments made from stated biology, with a rationale recorded for every
entry. That is still judgement -- it is a transparent, per-cell-type, individually
challengeable judgement instead of an opaque binary one, which is an improvement, not an
escape. A reader who disagrees with any single weight can change it and re-run.

Accordingly the deliverable is NOT one number. It is a SENSITIVITY ANALYSIS across several
weighting schemes, including the binary scheme as baseline and a shuffled scheme as a null.
If the conclusion only survives under one hand-built scheme, that is a finding and must be
reported as such.
"""

# w, rationale. w in [0, 1]: 0 = fully expendable, 1 = catastrophic irreversible loss.
WEIGHTS = {
    # ---- fully or near-fully expendable ------------------------------------------------
    'b-cells': (0.00, 'B-cell aplasia managed routinely with immunoglobulin replacement'),
    'plasma cells': (0.05, 'as B cells; residual loss of long-lived antibody memory'),
    'syncytiotrophoblasts': (0.00, 'placental; absent outside pregnancy'),
    'extravillous trophoblasts': (0.00, 'placental; absent outside pregnancy'),
    'hofbauer cells': (0.00, 'placental macrophages; absent outside pregnancy'),
    'decidual stromal cells': (0.05, 'endometrial/placental; cyclically regenerated'),
    'oocytes': (0.10, 'fertility loss; not life-threatening'),
    'early spermatids': (0.10, 'fertility loss; not life-threatening'),
    'late spermatids': (0.10, 'fertility loss; not life-threatening'),
    'early primary spermatocytes': (0.10, 'fertility loss; not life-threatening'),
    'differentiating spermatogonia': (0.10, 'fertility loss; not life-threatening'),
    'epididymal principal cells': (0.10, 'fertility loss; not life-threatening'),
    'epididymal efferent duct absorptive cells': (0.10, 'fertility loss only'),
    'prostatic glandular cells': (0.10, 'organ dispensable and frequently already resected'),
    'breast lactating cells': (0.15, 'function confined to lactation'),
    'melanocytes': (0.15, 'loss presents as vitiligo, which is tolerated clinically'),
    'endometrial glandular cells': (0.20, 'cyclically shed and regenerated'),
    'endometrial secretory cells': (0.20, 'cyclically shed and regenerated'),

    # ---- regenerating epithelium: renews, but barrier/function gap matters -------------
    'foveolar cells': (0.30, 'gastric mucosa renews; basis of approved claudin-18.2'),
    'tuft cells': (0.30, 'rare chemosensory epithelial population; renews'),
    'esophageal apical cells': (0.35, 'stratified epithelium renews from basal layer'),
    'basal keratinocytes': (0.35, 'skin renews; barrier loss survivable with care'),
    'salivary duct cells': (0.40, 'renews; xerostomia is morbid but survivable'),
    'colonocytes': (0.50, 'renews fast BUT barrier fails during the gap -- CEACAM5 colitis'),
    'enterocytes': (0.50, 'as colonocytes; barrier failure is the binding constraint'),
    'cholangiocytes': (0.70, 'limited regeneration; biliary injury poorly tolerated'),

    # ---- haematopoietic, but NOT uniformly expendable ----------------------------------
    'erythrocyte progenitors': (0.40, 'marrow regenerates; transfusion support exists'),
    'megakaryocyte-erythroid progenitors': (0.45, 'as above; bleeding risk adds'),
    'megakaryocyte progenitors': (0.45, 'platelet transfusion exists but bleeding is acute'),
    'megakaryocytes': (0.45, 'platelet transfusion exists but bleeding is acute'),
    'platelets': (0.45, 'transfusible; acute bleeding risk'),
    'monocyte progenitors': (0.50, 'regenerates; transient innate deficit'),
    'monocytes': (0.50, 'regenerates; transient innate deficit'),
    'neutrophils': (0.60, 'neutropenia is dangerous; G-CSF and recovery mitigate'),
    'neutrophil progenitors': (0.60, 'as neutrophils'),
    'mast cells': (0.30, 'loss well tolerated'),
    'nk-cells': (0.40, 'regenerates; partial innate deficit'),
    'innate lymphoid cells': (0.40, 'regenerates; partial innate deficit'),
    'cdc': (0.40, 'conventional dendritic cells; regenerate from marrow'),
    'pdcs': (0.40, 'plasmacytoid dendritic cells; regenerate from marrow'),
    'kupffer cells': (0.40, 'liver macrophages repopulate'),
    'thymocytes': (0.60, 'developing T cells; depletion impairs T-cell reconstitution'),
    't-cells': (0.70, 'T-cell aplasia is severe immunodeficiency with no routine replacement'),
    'thymic myoid cells': (0.50, 'rare thymic stromal population'),
    'microglia': (0.90, 'CNS-resident; do NOT repopulate from marrow; loss is permanent'),

    # ---- non-regenerating or functionally critical --------------------------------------
    'hepatocytes': (0.50, 'liver has substantial regenerative capacity'),
    'fibroblasts': (0.50, 'broadly distributed; partial regeneration'),
    'pancreatic acinar cells': (0.60, 'exocrine insufficiency is replaceable but morbid'),
    'fibro-adipogenic progenitors': (0.60, 'marrow/muscle stroma; supports other lineages'),
    'myosatellite cells': (0.60, 'muscle stem pool; loss impairs repair'),
    'myonuclei': (0.60, 'post-mitotic muscle'),
    'pituitary stem cells': (0.50, 'endocrine output largely replaceable by hormone therapy'),
    'lactotrophs': (0.40, 'prolactin loss is tolerable and replaceable'),
    'proximal tubule cells': (0.70, 'partial regeneration; renal failure is life-limiting'),
    'epicardial cells': (0.70, 'cardiac; minimal regeneration'),
    'ocular epithelial cells': (0.70, 'corneal damage impairs vision; partial regeneration'),
    'conjunctival goblet cells': (0.70, 'ocular surface integrity; poor regeneration'),
    'lymphatic endothelial cells': (0.70, 'poor regeneration; lymphoedema is permanent'),
    'vascular smooth muscle cells': (0.80, 'vascular integrity; no meaningful regeneration'),
    'choroid plexus epithelial cells': (0.80, 'CSF production; no regeneration'),
    'alveolar cells type 1': (0.90, 'gas exchange surface; minimal regeneration'),
    'podocytes': (0.90, 'terminally differentiated; loss drives irreversible renal failure'),
    'pancreatic islet cells': (0.90, 'insulin-producing; loss is permanent diabetes'),
    'oligodendrocytes': (0.90, 'CNS myelination; very limited regeneration'),
    'oligodendrocyte progenitor cells': (0.85, 'some remyelination capacity'),
    'retinal pigment epithelial cells': (0.95, 'vision loss is permanent'),
    'retinal amacrine cells': (0.95, 'neural retina; no regeneration'),
    'retinal bipolar cells': (0.95, 'neural retina; no regeneration'),
    'cone photoreceptor cells': (0.95, 'vision loss is permanent'),
    'brain inhibitory neurons': (1.00, 'CNS neurons; no regeneration'),
    'cardiomyocytes': (1.00, 'negligible regeneration; loss is permanent'),
}

DEFAULT_WEIGHT = 0.70  # unlisted cell type: assume non-expendable rather than flatter it


def weight(cell_type, scheme='continuous'):
    """Expendability weight for a normal cell type under the named scheme."""
    if scheme == 'continuous':
        return WEIGHTS.get(cell_type, (DEFAULT_WEIGHT, 'unlisted'))[0]
    if scheme == 'binary':
        import target_labels as TL
        return 0.0 if TL.haematopoietic(cell_type) else 1.0
    if scheme == 'regeneration_only':
        # coarser: only asks whether the tissue renews, ignoring replaceability/survivability
        w, _ = WEIGHTS.get(cell_type, (DEFAULT_WEIGHT, ''))
        return 0.0 if w <= 0.35 else 1.0
    raise ValueError(f'unknown scheme {scheme!r}')


def score(top_share, home, scheme='continuous'):
    """Two-factor score: confinement, penalised by how much the home tissue is needed.

    Under 'binary' this reduces EXACTLY to the section 23 rule (+1.0 for a haematopoietic
    home), so the baseline is reproduced rather than approximated.
    """
    return top_share + (1.0 - weight(home, scheme))


def coverage(cell_types):
    """How many of the supplied cell types have an explicit weight (not the default)."""
    known = sum(1 for c in cell_types if c in WEIGHTS)
    return known, len(set(cell_types))
