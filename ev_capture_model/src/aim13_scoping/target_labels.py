"""
FROZEN LABEL TABLE for the whole-body restriction analyses.

Every label here is drug-safety history from the published clinical record. NOTHING about
expression enters a label. This file is the object that an external reviewer checks; it is
deliberately separated from every analysis script so that the analysis cannot quietly edit it.

Columns per entry:
  gene      HGNC symbol (the key used against the Human Protein Atlas reference)
  common    the name clinicians use
  agent     the specific drug/construct the verdict is about
  outcome   'safe'  = reached approval or documented clinical tolerability, with no
                      dose-limiting on-target normal-tissue toxicity
            'toxic' = documented SEVERE, DOSE-LIMITING or FATAL on-target toxicity in a
                      normal tissue
            'efficacy_fail' = the programme failed, but for lack of efficacy, NOT toxicity.
                      EXCLUDED from the toxicity test (it is not evidence either way) and
                      used only as a separate control group.
  evidence  one-line citation of the clinical fact
  home_expected  the normal tissue the toxicity did/would fall on (documentation only;
                 never used as a predictor)

PRE-STATED EXCLUSIONS (fixed before any number was computed):
  - non-protein antigens (GD2, Lewis-Y, gangliosides): no gene, no RNA reference.
  - targets whose only verdict is 'still in early trials': no outcome to predict.
  - CD123/IL3RA: BOTH an approval (tagraxofusp, BPDCN) and a toxicity failure
    (talacotuzumab) exist for the same gene. Ambiguous by our own rule -> excluded.
  - PTPRC (CD45): no agent ever reached patients as a cytotoxic target. Retained ONLY as a
    negative control for the restriction measure, excluded from the outcome test.
"""
import re

TARGETS = [
 # gene, common, agent, outcome, evidence, home_expected
 ('CD19','CD19','tisagenlecleucel / axi-cel','safe','Approved CAR-T; on-target B-cell aplasia is managed with immunoglobulin replacement','B cells'),
 ('MS4A1','CD20','rituximab','safe','Approved 1997; B-cell depletion clinically managed','B cells'),
 ('CD22','CD22','inotuzumab ozogamicin','safe','FDA approved 2017 (ALL); B-cell depletion is on-target and clinically managed','B cells'),
 ('CD79B','CD79b','polatuzumab vedotin','safe','FDA approved 2019 (DLBCL); no dose-limiting normal-tissue on-target toxicity','B cells'),
 ('TNFRSF17','BCMA','ide-cel / cilta-cel','safe','Approved 2021; plasma-cell aplasia managed','plasma cells'),
 ('SLAMF7','CS1','elotuzumab','safe','FDA approved 2015 (myeloma); well tolerated in combination regimens','plasma cells / NK'),
 ('CD38','CD38','daratumumab','safe','FDA approved 2015 (myeloma); on-target plasma-cell loss is clinically managed','plasma cells'),
 ('TNFRSF8','CD30','brentuximab vedotin','safe','FDA approved 2011 (Hodgkin, ALCL); neuropathy is payload-driven, not on-target','activated lymphocytes'),
 ('CD33','CD33','gemtuzumab ozogamicin','safe','Approved 2000/2017; myelosuppression is on-target but marrow regenerates','myeloid'),
 ('GPRC5D','GPRC5D','talquetamab','safe','Approved 2023; on-target skin/nail/taste effects are common but not dose-limiting','plasma cells, keratinised tissue'),
 ('FOLH1','PSMA','177Lu-PSMA-617','safe','FDA approved 2022 (mCRPC radioligand); salivary uptake is on-target but not dose-limiting','prostate epithelium'),
 ('TACSTD2','TROP2','sacituzumab govitecan','safe','FDA approved 2020 (triple-negative breast); toxicity is payload-driven neutropenia','epithelium'),
 ('NECTIN4','Nectin-4','enfortumab vedotin','safe','Approved 2019 (skin effects on-target but not dose-limiting)','urothelium/skin'),
 ('CLDN18','claudin-18.2','zolbetuximab','safe','Approved 2024; nausea/vomiting on-target on gastric mucosa, manageable','gastric mucosa'),
 ('DLL3','DLL3','tarlatamab','safe','FDA approved 2024 (small-cell lung); DLL3 is essentially absent from adult normal tissue','neuroendocrine'),
 ('MET','c-MET','telisotuzumab vedotin','safe','FDA approved 2025 (c-MET-overexpressing non-squamous NSCLC); no dose-limiting on-target normal toxicity','epithelium'),
 ('MSLN','mesothelin','anti-MSLN CAR-T / anetumab','safe','Repeatedly tolerable in trials; serositis only at high intrapleural dose','mesothelium'),
 ('GPC3','glypican-3','GPC3 CAR-T','safe','Tolerable in hepatocellular trials','fetal liver (absent in adult)'),
 ('CD52','CD52','alemtuzumab','safe','Approved; profound but managed lymphodepletion','lymphocytes'),

 ('ERBB2','HER2','HER2 CAR-T (Morgan 2010)','toxic','Fatal pulmonary event within hours of infusion; trastuzumab also causes on-target cardiotoxicity','lung epithelium, cardiomyocytes'),
 ('EGFR','EGFR','anti-EGFR CAR-T / cetuximab','toxic','Dose-limiting on-target skin and pulmonary toxicity','keratinocytes, lung'),
 ('CEACAM5','CEA','CEA CAR-T (Parkhurst 2011)','toxic','Severe transient colitis, dose-limiting','colonocytes'),
 ('CA9','CAIX','CAIX CAR-T (Lamers 2006)','toxic','Grade 2-4 on-target liver enzyme toxicity, trial halted','bile duct epithelium'),
 ('CD47','CD47','magrolimab','toxic','On-target anaemia; programme discontinued 2024','erythrocytes, ubiquitous'),
 ('CD44','CD44v6','CD44v6 CAR-T / bivatuzumab','toxic','Fatal on-target skin toxicity with bivatuzumab mertansine','keratinocytes'),
 ('NCAM1','CD56','lorvotuzumab mertansine','toxic','Dose-limiting on-target peripheral neuropathy; programme halted','neurons'),
 ('EPCAM','EpCAM','catumaxomab / adecatumumab','toxic','On-target hepatotoxicity; catumaxomab withdrawn','hepatocytes/biliary'),
 ('FOLR1','FRalpha','mirvetuximab soravtansine','toxic','Dose-limiting on-target ocular keratopathy (boxed warning) despite approval','corneal epithelium'),
 ('FAP','FAP','FAP CAR-T','toxic','On-target toxicity to FAP+ bone-marrow stromal cells; cachexia and marrow failure','marrow stroma'),
 ('F3','tissue factor','tisotumab vedotin (Tivdak)','toxic','FDA approved 2021 (accelerated) / 2024 (full, cervical cancer); label carries a BOXED WARNING for ocular toxicity -- up to 60%% of patients affected (conjunctivitis, dry eye, corneal ulceration/damage, blepharitis), dose-limiting and requiring a dedicated eye-care REMS. Approved despite the toxicity, not because it is absent -- same pattern as FOLR1/mirvetuximab.','corneal/conjunctival epithelium'),

 ('FCER2','CD23','lumiliximab','efficacy_fail','Phase 3 stopped for futility, not toxicity','B cells'),
 ('GPNMB','GPNMB','glembatumumab vedotin (CDX-011)','efficacy_fail','METRIC phase 2b randomized trial in gpNMB-overexpressing triple-negative breast cancer explicitly did not meet its primary endpoint (PFS); development not advanced to approval. No dose-limiting on-target normal-tissue toxicity signal reported -- a clean efficacy failure, not a safety failure.','melanocytes / osteoclasts (low-level normal expression)'),
 ('CD40','CD40','dacetuzumab','efficacy_fail','Stopped for lack of efficacy','B cells / APC'),
 ('CD80','CD80','galiximab','efficacy_fail','Stopped for lack of efficacy','B cells / APC'),
 ('HLA-DRA','HLA-DR','apolizumab','efficacy_fail','Stopped for lack of efficacy','all APC'),
 ('SDC1','CD138','indatuximab ravtansine','efficacy_fail','Development discontinued, efficacy','plasma cells / epithelium'),
 ('MUC1','MUC1','several','efficacy_fail','No agent has reached approval; failures attributed to efficacy/shedding','glandular epithelium'),
 ('MUC16','CA125','oregovomab','efficacy_fail','Phase 3 negative for efficacy; no on-target normal-tissue toxicity signal reported','coelomic epithelium'),
 ('SLC34A2','NaPi2b','lifastuzumab vedotin','efficacy_fail','Did not beat comparator. NOTE (added Sept 2026): a second NaPi2b agent, upifitamab rilsodotin, was independently discontinued for efficacy (UPLIFT trial missed its pre-specified response-rate threshold, 15.6%% ORR) after an earlier FDA partial clinical hold for 5 grade-5 bleeding events (4 treatment-related) in ~560 patients. The bleeding mechanism was never publicly resolved as on-target NaPi2b toxicity, payload toxicity, or otherwise -- every source found states the cause was \"under evaluation\" and the story moves on to the efficacy result without returning to it. Kept as efficacy_fail, not relabeled toxic, because the program\'s own stated ending was efficacy and the toxicity mechanism was never confirmed on-target -- but flagged here so this label is not read as a clean safety record.','lung/thyroid'),
 ('CD74','CD74','milatuzumab','efficacy_fail','Development discontinued','B cells / APC'),
 ('ERBB3','HER3','patritumab deruxtecan','efficacy_fail','FDA BLA voluntarily withdrawn 2026 after the confirmatory phase 3 HERTHENA-Lung02 trial (EGFR-mutated NSCLC) missed its overall-survival endpoint despite an earlier positive PFS signal; Daiichi Sankyo/Merck attributed the withdrawal to the OS miss, not a safety finding. Two ILD-related deaths occurred in the treatment arm but were not cited as the reason for withdrawal -- a clean efficacy failure by the sponsors\' own account, same pattern as GPNMB.','epithelium (broadly expressed, low level)'),
 ('AXL','AXL','enapotamab vedotin','efficacy_fail','Genmab discontinued development in 2020 after early-phase expansion cohorts; company statement: did not meet the "stringent criteria for proof-of-concept," with "some evidence of clinical activity" that could not be improved by alternate dosing or biomarker selection. Explicitly an efficacy-driven decision, no toxicity concern cited.','epithelium / mesenchymal (broadly expressed)'),

 ('PTPRC','CD45','(none)','excluded_control','No cytotoxic agent ever reached patients; retained as a restriction negative control','all leukocytes'),
 ('IL3RA','CD123','tagraxofusp vs talacotuzumab','excluded_ambiguous','One approval and one toxicity failure on the same gene','myeloid/pDC'),
]

# --- objective, NON-CIRCULAR lineage class of the normal home cell type -------------------
# Assigned from the identity of the cell type alone. A haematopoietic home means the cell
# type is continuously regenerated from marrow stem cells AND its loss has an established
# clinical replacement (immunoglobulin, transfusion, growth factor). No drug outcome is used.
HAEM_PATTERNS = (
    r'\bb[- ]cells?\b', r'\bplasma cells?\b', r'\bt[- ]cells?\b', r'\bnk[- ]cells?\b',
    r'\bmonocytes?\b', r'\bmacrophages?\b', r'\bgranulocytes?\b', r'\bneutrophils?\b',
    r'\beosinophils?\b', r'\bbasophils?\b', r'\bdendritic\b', r'\bcdc\b', r'\bpdcs?\b',
    r'\berythroid\b', r'\berythrocytes?\b', r'\bmegakaryocytes?\b', r'\blangerhans cells?\b',
    r'\bkupffer cells?\b', r'\bmicroglia\b', r'\bmast cells?\b', r'\bmyeloid\b',
    r'\blymphocytes?\b', r'\bosteoclasts?\b', r'\bhofbauer cells?\b', r'\bthymocytes?\b',
)
_HAEM_RE = re.compile('|'.join(HAEM_PATTERNS), re.I)


def haematopoietic(cell_type_name):
    """True when the normal home cell type is a blood/marrow-derived lineage.

    Matching is word-boundary anchored on purpose: a plain substring test calls
    'goblet cells' haematopoietic because it contains 't cell'. That bug was found and
    fixed; the guard test below keeps it fixed.
    """
    return bool(_HAEM_RE.search(cell_type_name))


_GUARD = {'goblet cells': False, 'conjunctival goblet cells': False, 'b-cells': True,
          'plasma cells': True, 'kupffer cells': True, 'basal cells': False,
          'mast cells': True, 'alveolar cells type 1': False, 'colonocytes': False,
          'neutrophils': True, 'cdc': True, 'oligodendrocytes': False,
          'epicardial cells': False, 'foveolar cells': False}
for _k, _v in _GUARD.items():
    assert haematopoietic(_k) == _v, ('haematopoietic() guard failed on %r' % _k)


def by_outcome(*wanted):
    return [t for t in TARGETS if t[3] in wanted]
