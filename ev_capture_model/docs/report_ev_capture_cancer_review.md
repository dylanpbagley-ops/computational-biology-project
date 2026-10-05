# Extracellular Vesicle Capture in Cancer Detection and Treatment: A Literature-Grounded Landscape and a Portfolio of Computational Research Topics

**Prepared for:** Dylan Bagley · **Date:** 2026-08-26
**Purpose:** Map the research landscape at the intersection of extracellular vesicle (EV) capture, cancer detection, and cancer treatment; identify where the evidence is strong, contested, or thin; and derive a portfolio of novel thesis-scale research topics that are executable with a **mainly computational workflow**.

---

## 1. Scope and method

This review synthesizes 258 unique records retrieved through 15 targeted literature searches spanning: microfluidic capture; label-free physical separation (acoustofluidics, deterministic lateral displacement, nanowires, dielectrophoresis); aptamer-based affinity capture; glycan/lectin capture; computational fluid dynamics (CFD) and transport modeling of capture; EV cargo biomarkers; machine-learning (ML) diagnostics; single-EV analysis; multi-analyte liquid-biopsy integration; EV PD-L1 and immunotherapy monitoring; minimal residual disease (MRD) and longitudinal monitoring; therapeutic EV removal/scavenging; engineered EVs for therapy; computational aptamer design; and EV databases/standards. The window was ~2020–2026 with foundational older work retained (earliest record 2008 [169]).

Thirteen pivotal papers were read in full text ([12, 47, 54, 57, 60, 63, 130, 133, 156, 157, 162, 216]); all other records are represented by abstracts and structured metadata. A machine-readable evidence table with per-record facet, study type, citation count, and evidence source accompanies this report (`evidence_table_ev_capture_cancer.csv`). Inline numbers [N] refer to that table's `ref_index`.

---

## 2. Executive summary

- **Capture technology is prolific but empirically designed.** Dozens of microfluidic, affinity, and label-free devices report high capture efficiencies, but only one study in our corpus uses computational optimization of capture geometry — and it showed that intuition (maximize surface area) is wrong [216]. Rational, simulation-driven device design is a near-empty niche.
- **ML diagnostics work impressively in retrospective pilots and almost nowhere else.** EV-protein panels [47], EVP proteomics [46], and SERS+deep learning [54, 39] report AUCs of 0.95–0.97, but the authors themselves flag small, single-site, ethnically homogeneous cohorts and absent external validation [47, 54]. Benchmarking and validation methodology is a genuine gap.
- **Single-EV measurement has outrun single-EV computation.** Platforms such as SEVEN [57], SEVEN-UP [68], and NOBEL-SPA [63] now resolve per-vesicle protein/miRNA colocalization, and SVAtlas aggregates single-EV omics [109] — yet methods for subpopulation deconvolution, integration, and statistical testing on such data are immature [65].
- **Therapeutic EV removal has crossed from concept to clinic without quantitative modeling.** Hemofiltration concepts from 2012 [157] are now in a safety/dose-finding trial combined with anti-PD-1 [167, 168], and in vivo dual-removal devices exist preclinically [162] — but no pharmacokinetic-style models of EV clearance, rebound, or dosing exist anywhere in the corpus.
- **Computational aptamer design is mature as a method but untouched for tumor-EV targets.** Diffusion models, transformers, and ML-guided refinement design aptamers in silico [143, 148, 152, 156], yet essentially none target EV surface proteins for capture reagents.
- **The data substrate for computational work is rich and underused.** exoRBase 3.0 [111], Vesiclepedia 2024 [119], SVAtlas [109], and EV-TRACK/MISEV metadata [122, 123, 125] are public, large, and — because isolation-method bias is documented [47, 110] — in need of harmonization methods.

Ten computational thesis topics derived from these gaps are presented in Section 5.

---

## 3. The landscape

### 3.1 Capture technologies

#### 3.1.1 Microfluidic and immunoaffinity capture

Microfluidics is the most densely populated capture modality, repeatedly reviewed as enabling selective, low-volume, integrated EV isolation and detection [70, 71, 72, 73, 85, 86]. Representative performance claims are strong: a phosphatidylserine-binding ExoChip reported ~90% capture of cancer-associated EVs, outperforming ultracentrifugation [78]; a graphene-oxide/polydopamine nano-IMEX chip achieved ultrasensitive detection [84]; the ExoSearch chip enabled continuous-flow multiplexed immunomagnetic isolation for ovarian cancer diagnostics [218]; the ExoProfile chip reported a limit of detection of 21 exosomes/mL with multiplexed immunophenotyping [220]; and a 3D-nanopatterned chip achieved ultrasensitive circulating-exosome detection in plasma [211]. Recent engineering emphasizes automation and release efficiency: the ExoCPR magnetic chip reports 75.8% capture and 62.7% release [82], the DSPE-Ti-Chip reports 86.9% recovery from microvolume biofluids [80], and the SpinEx automated disc isolates and labels 16 EV proteins from 150 µL of whole blood in under 75 minutes, reporting 90% cancer-vs-non-cancer accuracy and 96% five-tumor-type classification [87].

**Agreement:** immunoaffinity microfluidics can deliver high capture efficiency and multiplexed readout from small volumes. **Friction:** efficiency metrics are defined inconsistently across studies (spike-in recovery vs. endogenous yield vs. purity), which the MISEV guidelines and EV-TRACK framework explicitly flag as a field-wide reproducibility problem [125, 124, 122, 123].

#### 3.1.2 Label-free physical separation

Label-free modalities separate EVs by size, density, or polarizability without capture reagents: acoustofluidic isolation of exosomes from whole blood [237] and EV/lipoprotein separation [227]; deterministic lateral displacement (DLD) including nano-targeted DLD [240] and non-uniform DLD arrays enhanced with dielectrophoresis reporting >90% purity [224]; insulator-based dielectrophoresis for rapid small-EV purification [231]; viscoelastic [219] and elasto-inertial [207] flow fractionation; oxide-nanowire size-based manipulation [236]; and antibody-cocktail magnetic nanowires giving ~3× yield over conventional isolation [241]. Reviews position these as biocompatible, continuous, and contact-free [230, 232, 234, 228, 229, 235], with microfluidic size-exclusion chromatography emerging as an automatable option [79].

**Thin zone:** these papers are dominated by physics demonstrations; very few connect separated fractions to downstream clinical classification, and design choices (pillar geometry, acoustic frequency, flow rates) are made by intuition or one-factor-at-a-time experiments rather than optimization.

#### 3.1.3 Aptamer-based affinity capture

Aptamers are used as capture reagents against tetraspanins and tumor-associated surface proteins: dual-aptamer systems for direct tumor-exosome quantification [20]; apta-magnetic lab-on-a-chip pre-concentration for colorectal cancer exosomes [21]; on-chip isolation with protein analysis [22]; tetrahedral-DNA microelectrode sensors [23]; electrochemical [24] and SERS [25] aptasensors; target-specific aptamer microfluidics capturing ~10^7–10^8 exosomes/mL within 20 minutes [222]; and multivalent aptamers achieving ~45% isolation of spiked vesicles in 40 minutes with an 84.7% reduction of albumin contamination [37]. Target-specific aptamers include CD63 [29], breast-cancer-EV-selective ex-50.T [27], and prostate-cancer-EV Seq25 [28]. Reviews cover aptamer EV theranostics [26], aptamer liquid biopsy [30], aptamer-microfluidic biosensors [31, 35], EV isolation/analysis/therapeutics [32], and — critically for downstream cargo work — methods for eluting intact EVs from aptamer affinity surfaces [36]. Notably, aptamers also appear on the therapeutic side: the A8 peptide aptamer blocks HSP70 on tumor-derived exosomes to restore anticancer immune responses [33].

**Gap:** every capture aptamer in this corpus was discovered by experimental SELEX; none was designed or optimized computationally, despite a mature in silico aptamer literature (Section 3.1.5's sibling, below).

#### 3.1.4 Glycan/lectin-mediated capture

A distinct capture axis exploits cancer-associated glycosylation rather than protein epitopes. Lectins are reviewed as EV biomarker-discovery tools [191, 194, 203]. Demonstrations include: EVLET, a lectin-based thermophoretic glycan-profiling assay for triple-negative breast cancer reporting 91% diagnostic and 96% response-monitoring accuracy [192]; lectin-affinity glycosylation-pattern analysis distinguishing gastric cancer stages [193]; lectin microarrays linking exosomal glycan patterns to colorectal cancer subtypes and enabling dynamic monitoring [195]; lectin-SERS multiplex glycan profiling [196]; O-glycan-altered EVs as a pancreatic cancer serum marker [198]; Helix pomatia agglutinin binding elevated in breast cancer EVs [199]; consistent LGALS3BP enrichment on cancer-derived EVs [200]; GlyExo-Capture, LCA-immobilized magnetic enrichment of fucosylated EVs capturing ~60% of target vesicles, with a downstream 4-miRNA signature for hepatocellular carcinoma [201]; distinct glycosignatures in high-grade serous ovarian carcinoma EVs [202]; and 756 intact N-glycopeptides (344 differentially expressed) in urinary EVs in HCC glycoproteomics [205]. Dielectrophoretic capture of lectin-bound Janus nanoparticles extends glycan capture to devices [190].

**Gap:** glycan signatures are generated study-by-study, one cancer type at a time; no cross-cancer computational mining of EV glycan data exists, and glycan capture has never been integrated with the ML-diagnostic pipeline mainstream.

#### 3.1.5 Computational modeling of capture devices

This is the thinnest and most revealing facet. The corpus contains transport/CFD work adjacent to EV capture — predictive nanoparticle capture on nanoporous membranes [206], elasto-inertial nanoparticle separation simulation [207], Euler–Lagrange reactive-separation modeling [210], coarse-grained MD of exosome squeezing through nanopores for drug loading [209, 223], transvascular nanoparticle transport simulation [215], and CFD reviews for nanoparticle–tumor interaction [208] — but **only one study optimizes an EV capture device computationally**: Hisey et al. coupled Python parallel pattern search with COMSOL particle tracing on HPC to optimize triangular micropillar arrays, finding counterintuitively that *larger pillars with wider spacing* outperform maximum-pillar-density designs (streamlines follow pillar contours), yielding +140% predicted capture and +102%/+62% experimentally validated improvement (fluorescence/NTA) of best vs. worst geometries [216].

**Implication:** the entire capture-device field optimizes by trial and error; a single proof-of-concept shows simulation-driven optimization is feasible, HPC-tractable, and capable of overturning design intuition. This is arguably the most under-exploited computational niche in the landscape.

### 3.2 Detection and diagnostics

#### 3.2.1 EV cargo biomarkers

Tumor-derived EV cargo — miRNA, mRNA, protein, DNA — is extensively reviewed as a liquid-biopsy substrate [2, 11, 12, 13, 14, 17, 19, 34]. Specific signatures include a four-miRNA EV panel for early diagnosis [1]; colorectal EV miR-182-5p/miR-196b-5p/miR-429 [9]; plasma small-EV miRNA signatures for early lung cancer [16]; blood-derived exosomal tumor RNA signatures for multi-cancer detection [18]; CRISPR/Cas13a EV miRNA detection [3]; and encoded fusion-mediated multiplex miRNA profiling [10]. The landmark EVP proteomics study (Hoshino et al., Cell 2020; read at abstract level) reported that EV/particle protein biomarkers define multiple human cancers with 90% sensitivity and 95% specificity [46]. Reviews consistently caution that purification quality and cargo attribution remain confounders [13].

#### 3.2.2 Machine-learning diagnostics

ML on EV measurements is the fastest-moving diagnostic thread. Read in full text:

- **Hinestrosa et al. 2022** [47]: alternating current electrokinetic (ACE/Verita) EV isolation (240 µL plasma, ~2 h, no immunoaffinity) plus a 13-EV-protein panel and logistic regression with recursive feature elimination, resampled 100× on 67/33 splits; N=323 (139 stage I–II cases across pancreatic/ovarian/bladder; 184 controls). **AUC 0.95, 71.2% sensitivity at 99.5% specificity**; stage-I sensitivity was 95.5% (pancreatic), 74.4% (ovarian), 43.8% (bladder); EV-associated proteins outperformed plasma free proteins (AUC 0.95 vs 0.87). The authors state the limitations plainly: pilot scale, ethnically homogeneous cohort, no external validation [47].
- **Shin et al. 2023** [54]: Exosome-SERS-AI — SEC isolation, AuNP SERS, CNN with multiple-instance learning; 520 held-out test samples across six cancer types; **AUC 0.970 for cancer presence, mean AUC 0.945 for tissue-of-origin, 90.2% sensitivity / 94.4% specificity, 88.1% early-stage sensitivity**. Stated limitations: training-set size, need for prospective validation, spectral interpretability, chip manufacturability [54].

Around these anchors: deep-learning SERS diagnosed early lung cancer at ~95% accuracy [39]; label-free SERS+ML for lung cancer [40]; ML classification of pancreatic neoplasia from plasma EVs [41]; deep-learning AFM nanomechanical classification of NSCLC EVs [42]; nanoaperture optical tweezers + deep learning for label-free classification [44]; EV microbiome markers + ML for pancreatic cancer [45]; plasma EV metabolomics with RF/SVM for ovarian cancer [50]; a four-protein EV panel (ECM1, MBL2, BTD, RAB5C) for TNBC recurrence [51]; an interpretable deep neural network on EV RNA for liver cancer (AUC 0.8877) [52]; ML on exosomal surface sialic acid SERS spectra (93% accuracy, ovarian) [53]; and the metaEVchip Kolmogorov–Arnold-network nanoplasmonic sensor reporting an **external-validation AUC of 0.99** for pancreatic cancer — a rare instance of true external validation [49]. Reviews map the ML-EV omics space [43, 48, 179].

**Conflict/tension:** reported performance is uniformly high, but nearly every study is retrospective, single-site, and small; the two full-text anchors both concede the absence of external validation [47, 54], while [49] shows external validation is achievable. The field lacks shared benchmarks, leakage-aware evaluation protocols, and cross-study replication — a methodological vacuum that is itself a research opportunity.

#### 3.2.3 Single-EV analysis

Single-EV platforms now resolve heterogeneity invisible to bulk assays. Read in full text: **SEVEN** (single-EV nanoscopy) quantified per-vesicle size and tetraspanin content in pancreatic cancer plasma, finding CD9-enriched EVs of 58–59 nm in PDAC vs 60 nm in healthy donors (MRPS: 67–72 nm vs 70 nm), with stromal CD9 linked to poor prognosis [57]; **NOBEL-SPA** colocalized protein and miRNA markers on single EVs to reveal early-breast-cancer-associated subpopulations that bulk analysis misses entirely [63]; and a Nature Nanotechnology review frames heterogeneity exploitation as the field's diagnostic/therapeutic frontier [60]. Supporting work: deep-learning multi-miRNA profiling in single EVs [38]; multiomics of large oncosomes [56]; cryo-EM architectural heterogeneity with computational analysis [61]; single-cell-secreted EV tissue fingerprints [62]; mutant KRAS/BRAF effects on single-vesicle uptake [67]; and SEVEN-UP, an affordable universal SRRF-based single-EV imaging protocol [68]. Critically, a 2025 perspective argues the field still lacks a standardized definition of "single EV research" and rigorous heterogeneity quantification [65], and SVAtlas now aggregates single-EV omics into a public resource [109].

**Gap:** measurement platforms are producing high-dimensional per-vesicle data faster than computational methods exist to analyze it — subpopulation deconvolution, batch correction, multi-view (protein+miRNA) integration, and statistical testing at the single-vesicle level are all open.

#### 3.2.4 Multi-analyte integration

Combining EVs with cfDNA, CTCs, and protein markers is actively explored: an EV miRNA/mRNA + cfDNA + CA19-9 panel improved pancreatic cancer diagnosis and staging [174]; tri-modal CTC+exosome+cfDNA analysis [180]; the ONCE single-aliquot EV-RNA+cfDNA approach improved HER2-positivity detection [182]; EV-based multianalyte biopsy is reviewed as improving early-stage sensitivity [178]; and multi-omics liquid biopsy with ML is surveyed broadly [172, 173, 176, 181, 183, 185, 187]. **Direct conflict:** Thakur et al. found EV-DNA offered *no* advantage over cfDNA for colon cancer mutation detection (61.9% vs 66.7% sensitivity) [189], against studies reporting additive value in other contexts [174, 182]. EV-DNA biology itself remains unsettled a decade after its discovery [186, 188].

**Implication:** whether and when multi-analyte fusion helps is unresolved and likely cancer-type- and stage-dependent; principled fusion methods with honest benchmarking are missing.

#### 3.2.5 EV PD-L1 and immunotherapy monitoring

EV-borne immune checkpoints are a clinically motivated monitoring axis. Read in full text: circulating PD1+ EVs drive anti-PD1 resistance in metastatic melanoma, and liquid-biopsy tracking of PD1/PD-L1+ EVs predicts response [130]; longitudinal exosomal PD-L1 tracking outperformed tissue PD-L1 for monitoring melanoma patients [133]. Around these: SPR-ExoPD-L1 [126] and HOLMES-ExoPD-L1 [135] assays; exosomal PD-L1 predicting immunotherapy efficacy in NSCLC [132]; PD-L1 on *large* EVs predicting response specifically in tissue-PD-L1-low/negative NSCLC [138]; a validated multiomic model combining plasma EV PD-L1 dynamics with radiomics to predict checkpoint-inhibitor non-response [137]; and reviews of exosomal PD-L1 biology, mechanisms of checkpoint-blockade resistance, and clinical implications [127, 128, 129, 131, 134, 136, 139, 140]. The therapeutic corollary — removing exosomal PD-L1 to sensitize tumors to anti-PD-1 — is explicitly proposed [128] and connects this facet to Section 3.3.1.

#### 3.2.6 MRD and longitudinal monitoring

EVs are proposed as MRD and recurrence-monitoring substrates across the cancer journey [243, 247, 248, 250, 251]: myeloma MRD biomarkers in blood EVs [245]; longitudinal multi-parametric liquid biopsy in metastatic breast cancer [246]; small EVs from surviving cancer cells as multiparametric MRD monitors [249]; EV MGMT methylation in IDH-wildtype glioblastoma at 85.7% sensitivity [256]; miR-203a-3p in tumor-draining-vein EVs as a lung cancer relapse biomarker [257]; and four targetable myeloma EV surface subpopulations (BCMA, CD38, GPRC5D, CD319) tracking CAR-T response depth [258]. Methodologically, tumor-EV digital scoring assays provide absolute RNA-based quantification suitable for longitudinal monitoring [253], and a 2024 consensus statement calls for standardization of EV liquid biopsy to advance the clinic [252]. General liquid-biopsy challenges (cost, standardization, accessibility) apply [242, 244].

**Gap:** longitudinal EV data exist, but quantitative frameworks for thresholds, dynamics, and lead-time estimation (analogous to ctDNA MRD modeling) do not.

#### 3.2.7 Databases, standards, and computational resources

The public data substrate is substantial: **Vesiclepedia 2024** catalogs 62,822 quantitative entries across 47 studies (50,550 RNA; 566,911 protein; 3,839 lipid; 192 metabolite; 167 DNA entries) [119], continuing the Vesiclepedia lineage [118, 120]; **exoRBase** has reached version 3.0, an atlas of extracellular RNAs from human biofluids for context-specific biomarker identification [111, 107]; **SVAtlas** integrates single-EV omics [109]. Standards and reporting frameworks include MISEV2018 [125] and its update process [124], and the EV-TRACK knowledgebase for transparent experimental metadata [122, 123]. The EV RNA-seq community has published explicit bioinformatics challenges and recommendations, flagging the lack of standardized pipelines for delineating EV-associated small RNAs [110]. Multi-omics integration of EV data is repeatedly called for [108, 113, 115, 116], and EV proteo-transcriptomic signatures are known to diverge from their cell of origin [112].

**Gap:** these resources are siloed; isolation-method bias (documented experimentally, e.g., ACE vs ultracentrifugation protein contamination differences [47]) propagates into every atlas, and no cross-method harmonization or batch-correction framework exists for EV omics.

### 3.3 Treatment

#### 3.3.1 Therapeutic EV removal and scavenging

Removing tumor-derived EVs is a mechanistically motivated therapeutic strategy: tumor exosomes drive immune suppression [161, 171] and were proposed as an immune-escape mechanism whose extracorporeal removal could restore immunity as early as 2008 [169]. Read in full text: **Marleau et al. 2012** introduced the ADAPT hemofiltration concept — affinity agents (antibodies against surface proteins, lectins against glycoproteins) capturing tumor exosomes from the entire circulation as a therapeutic adjuvant [157]; **Sun et al. 2023** reported a Pro-Metastatic Derivatives Eliminator (PMDE) achieving in vivo dual removal of circulating tumor cells and tumor-derived exosomes, reducing their biodistribution to distant organs [162]. The translational frontier: the Aethlon Hemopurifier clears exosomes with malignancy/immunosuppressive signatures from cancer patient plasma [158] and is now in a safety, feasibility, and dose-finding trial (3+3 design; one-to-three treatment cohorts) combined with pembrolizumab/nivolumab in solid tumors [167, 168]. Alternative scavenging modalities: a curvature-sensing peptide that disrupts tumor-derived exosomes and enhances immunotherapy [164]; bioorthogonal metabolic engineering reducing circulating tumor-cell-derived exosomes after chemotherapy [165]; and therapeutic plasma exchange with EV filtration reducing tumor-derived EVs ~75% in vitro and potentially enhancing anti-PD-L1 effects [170]. Targeted inhibition of tumor-derived exosomes is reviewed as a therapeutic option [160].

**Gap (striking):** despite clinical translation, there is no quantitative model anywhere in the corpus for EV production/clearance kinetics, post-filtration rebound, dosing schedules, or combination-therapy timing — the pharmacometrics of EV scavenging does not exist.

#### 3.3.2 Engineered EVs for therapy

The engineering literature is large and almost entirely experimental: source selection and engineering for cancer-targeted therapy [88, 93]; surface-modification strategies for targeting, circulation time, and payload release [89, 90, 91, 92, 95, 96, 97, 99, 102, 104]; therapeutic RNA delivery [98]; immunotherapy applications [94, 100]; bench-to-clinic translation [101]; biodistribution and organotropism engineering [103]; EGFR-targeting monobody–C1C2 fusion display on EVs [105]; and phospholipid-anchored ligand conjugation modifying >90% of milk-derived EVs [106]. Immune-cell-derived EVs (including dendritic-cell exosomes) form an immunotherapy sub-thread [159, 166].

**Computational opening:** ligand selection, tether/fusion design, and biodistribution prediction are done empirically; protein-design and PK-modeling tools have barely touched EV engineering.

---

## 4. Gap map: where the evidence is strong, contested, or thin

**Strong (convergent evidence):**
1. EV cargo and surface markers carry cancer-diagnostic information across many cancer types and assay modalities [46, 47, 54, 1, 16, 18].
2. Immunoaffinity microfluidic capture can achieve high efficiency and multiplexed readout from small volumes [78, 84, 211, 218, 220, 87].
3. EV PD-L1/PD1 track and predict immunotherapy response [130, 133, 132, 138].
4. Single-EV heterogeneity is real, measurable, and clinically patterned [57, 63, 60].

**Contested (direct conflicts):**
1. **Multi-analyte fusion value** — EV-DNA adds nothing over cfDNA in colon cancer [189] vs. additive value in pancreatic [174] and HER2 [182] settings.
2. **Capture-efficiency claims** — incomparable across studies due to inconsistent metrics, against MISEV/EV-TRACK standardization calls [125, 122].
3. **Design intuition** — maximize-surface-area pillar packing is *worse* than sparse large pillars under flow [216], contradicting the implicit heuristic of many chip designs.

**Thin (near-empty, computationally tractable):**
1. Simulation-driven capture-device optimization (one paper: [216]).
2. Computational design of EV-capture aptamers (methods exist [143, 146, 148, 152, 156]; EV targets untouched).
3. Single-EV computational methods (data platforms [57, 63, 68, 109] outpace analysis methods [65]).
4. Cross-cancer EV glycan mining (single-cancer studies only [192, 193, 195, 198, 201, 202]).
5. Pharmacokinetic-style modeling of therapeutic EV removal (clinical trial running [167, 168]; zero models).
6. Harmonization/batch-correction of EV omics atlases across isolation methods [110, 111, 119, 122].
7. External-validation and benchmarking frameworks for ML-EV diagnostics [47, 49, 54].
8. Quantitative longitudinal/MRD modeling of EV signals [252, 253, 258].

---

## 5. Portfolio: ten computationally-tractable thesis topics

Each topic is scoped as a PhD/postdoc project executable primarily at a computer: public data, open tools, and at most optional collaboration for validation. Ordered roughly by novelty × feasibility.

### Topic 1 — ML-guided design of aptamer capture reagents for tumor-EV surface proteins
- **Rationale:** Capture aptamers are all SELEX-derived [20, 27, 28, 29], while computational aptamer design has matured independently: ML-guided refinement improved a 275 nM seed to 8 nM and enriched <128 nM binders 460-fold over particle display alone (187,499 aptamers tested) [156]; diffusion models (AptaDiff) [143], transformers (AptaTrans) [148], XGBoost interaction predictors (APIPred, 96.5% accuracy) [146], generative RBMs [149], and platforms like AptaBLE [152] are published and several are open-source.
- **Novelty:** first application of generative/predictive aptamer ML to EV surface targets (CD63, EpCAM, PD-L1, cancer-associated glycoprotein ectodomains); includes multivalency-aware design motivated by [37].
- **Data/tools:** AptaDB/SELEX-derived datasets used by [146, 148, 149]; target structures from AlphaFold DB/PDB; AptaTrans/APIPred/AptaDiff code; docking + MD for structure-based filtering [151, 155].
- **First experiments (3–6 mo):** (i) reproduce the [156] ML-refinement pipeline on public particle-display data; (ii) benchmark AptaTrans vs APIPred vs AptaDiff on held-out aptamer–protein pairs; (iii) generate and rank candidate aptamers against one EV target (e.g., CD63 extracellular loop) with docking/MD triage.
- **Risks:** EV-target training data are scarce → mitigate with transfer learning from protein–aptamer pairs; experimental validation requires a collaborator (plan as optional).

### Topic 2 — Simulation-driven optimization of EV capture devices with ML surrogates
- **Rationale:** the only computational device-optimization study showed counterintuitive optima and +102% experimental improvement [216]; the rest of the device literature is empirical [78, 80, 82, 84, 211, 218, 220]. Adjacent CFD/transport methods are available [206, 207, 210, 215].
- **Novelty:** extend [216] from geometry-only to *multi-physics* optimization — flow + surface reaction kinetics (affinity capture) + particle margination — and replace COMSOL-in-the-loop with a trained surrogate (graph neural network or Gaussian-process) enabling multi-objective optimization (capture efficiency vs. purity vs. shear stress on EVs).
- **Data/tools:** OpenFOAM or COMSOL; Python optimization (nevergrad/ax); published device geometries as baselines [78, 211, 216].
- **First experiments:** (i) reproduce the triangular-pillar optimization of [216] in an open-source stack; (ii) add a Langmuir surface-binding reaction and show how optimum shifts with kon/koff; (iii) train a surrogate on 500–2,000 simulated designs and benchmark against pattern search.
- **Risks:** COMSOL licensing → use OpenFOAM; experimental validation optional (collaboration), in silico Pareto fronts are publishable alone.

### Topic 3 — Cross-method harmonization of public EV omics atlases
- **Rationale:** exoRBase 3.0 [111], Vesiclepedia 2024 [119], and SVAtlas [109] aggregate data generated under heterogeneous isolation methods; EV RNA-seq bioinformatics is explicitly non-standardized [110]; isolation choice measurably changes cargo (ACE vs UC protein contamination [47]); EV-TRACK metadata quantify reporting heterogeneity [122, 123].
- **Novelty:** first systematic quantification and correction of isolation-method/platform batch effects across EV atlases — an "EV omics harmonization" methods paper analogous to batch-correction milestones in scRNA-seq.
- **Data/tools:** exoRBase 3.0, Vesiclepedia 2024, SVAtlas downloads; EV-TRACK metadata; ComBat/harmony/scVI-style methods; MISEV annotations [125].
- **First experiments:** (i) assemble cancer-vs-control EV RNA matrices from exoRBase stratified by isolation method; (ii) quantify method-driven variance (PERMANOVA) and show its effect on biomarker signatures; (iii) benchmark 3–4 batch-correction methods by biomarker-signature stability.
- **Risks:** metadata incompleteness → restrict to well-annotated subsets; purely computational, low risk.

### Topic 4 — Honest benchmarking and external-validation framework for ML-EV diagnostics
- **Rationale:** headline AUCs of 0.95–0.99 [47, 49, 54] coexist with author-acknowledged absence of external validation [47, 54]; [49] shows external validation is feasible but rare; reviews call for roadmaps [43, 48].
- **Novelty:** a reproducible benchmarking harness for EV-diagnostic ML — nested CV, leave-one-study-out evaluation, leakage audits (spectral preprocessing, feature selection inside CV), calibration analysis — applied across public EV datasets; quantifies the retrospective-to-external performance drop.
- **Data/tools:** public EV proteomic/RNA datasets (exoRBase [111], Vesiclepedia [119], study supplements from [47, 51, 52]); scikit-learn/pytorch; TRIPOD/PROBAST reporting checklists.
- **First experiments:** (i) reimplement two published classifiers (e.g., [47]-style protein panel, [52]-style RNA model) from public data; (ii) demonstrate performance inflation under naive vs. nested evaluation; (iii) release the harness as a community resource.
- **Risks:** raw-data availability varies → pre-register dataset list; high community value even with 3–4 datasets.

### Topic 5 — Computational deconvolution of single-EV subpopulations
- **Rationale:** single-EV platforms (SEVEN [57], SEVEN-UP [68], NOBEL-SPA [63]) and the SVAtlas resource [109] generate per-vesicle multi-marker data; bulk analysis misses disease-relevant rare subpopulations [63]; the field lacks standardized heterogeneity quantification [65].
- **Novelty:** mixture models / deep clustering tailored to single-EV count and intensity data; multi-view integration of protein+miRNA colocalization [63]; statistical tests for subpopulation frequency differences between clinical groups.
- **Data/tools:** SVAtlas [109]; published single-EV datasets [57, 63, 68]; scvi-tools-style probabilistic modeling; flow-cytometry-inspired gating benchmarks.
- **First experiments:** (i) benchmark clustering methods on SVAtlas; (ii) develop a noise-aware mixture model for per-vesicle marker counts; (iii) re-analyze a published cancer cohort to test subpopulation–stage associations.
- **Risks:** data sparsity per vesicle → simulate ground truth from platform noise models first.

### Topic 6 — Cross-cancer computational mining of EV glycan signatures
- **Rationale:** glycan/lectin capture shows strong single-cancer results (TNBC 91%/96% [192]; gastric staging [193]; CRC subtypes [195]; pancreatic O-glycans [198]; HCC fucosylation [201]; ovarian glycosignatures [202]; 344 differential N-glycopeptides in HCC urinary EVs [205]) but no cross-cancer synthesis.
- **Novelty:** first pan-cancer EV glycosignature atlas: harmonize lectin-array/glycoproteomic datasets, link EV glycan patterns to tumor glycosyltransferase expression (TCGA), and nominate lectin-capture panels per cancer type.
- **Data/tools:** published lectin microarray data [193, 195, 198]; glycoproteomics [205]; TCGA glycogene expression; glycan databases (GlyTouCan); ML classifiers.
- **First experiments:** (i) curate 4–6 public EV glycan datasets into a common feature space; (ii) pan-cancer vs. cancer-specific glycan-signature classification; (iii) correlate EV glycan signatures with matched-tissue glycosyltransferase expression.
- **Risks:** heterogeneous glycan feature spaces → use motif-level representations; moderate curation burden.

### Topic 7 — Principled multi-analyte fusion: when do EVs add value over cfDNA?
- **Rationale:** direct conflict in the literature — no EV-DNA advantage in colon cancer [189] vs. additive value in pancreatic panels [174] and HER2 detection [182]; EV-DNA biology unsettled [186, 188]; fusion reviews are qualitative [178, 172, 173].
- **Novelty:** a decision-theoretic fusion framework: early/late/joint fusion benchmarks on paired EV+cfDNA datasets, with value-of-information analysis identifying cancer types/stages where EV measurements add discriminative signal beyond cfDNA.
- **Data/tools:** public paired EV-cfDNA datasets (e.g., [182] ONCE data, [174] panel data); cfDNA fragmentomics/methylomics resources [185]; scikit-learn/pytorch.
- **First experiments:** (i) reproduce [189] and [182] analyses from public data; (ii) implement 3 fusion strategies and compare under matched evaluation; (iii) stratify added-value by stage and cancer type.
- **Risks:** few truly paired public datasets → augment with simulation calibrated to published effect sizes.

### Topic 8 — In silico pharmacokinetics of therapeutic EV scavenging
- **Rationale:** EV removal is clinically translated (Hemopurifier + anti-PD-1 trial [167, 168]; ADAPT concept [157]; PMDE in vivo [162]; plasma exchange ~75% EV reduction [170]; curvature-peptide disruption [164]) with **zero** quantitative models of EV production, clearance, rebound, or dosing.
- **Novelty:** first compartmental PK/PD-style model of circulating tumor-EV dynamics under extracorporeal scavenging: production rate from tumor burden, first-order clearance, device extraction ratio, rebound between sessions; simulate trial dosing cohorts from [168] and combination timing with checkpoint blockade [128, 130].
- **Data/tools:** published EV half-life/clearance estimates and trial protocols [157, 162, 167, 168, 170]; Python ODE solvers (scipy) or nlmixr2; parameter identifiability and sensitivity analysis.
- **First experiments:** (i) build the baseline two-compartment model parameterized from literature; (ii) fit extraction ratios to [170]'s ~75% in vitro reduction and [158]'s clearance data; (iii) simulate the [168] 3+3 cohorts and derive model-informed dosing hypotheses.
- **Risks:** sparse kinetic parameters → global sensitivity analysis and explicit uncertainty quantification; high novelty offsets parameter uncertainty.

### Topic 9 — Quantitative longitudinal models for EV-based MRD and response monitoring
- **Rationale:** EV MRD markers are proliferating (myeloma subpopulations tracking CAR-T response [258]; MGMT methylation [256]; digital scoring assays [253]; relapse miRNAs [257]) and consensus calls for standardization [252], but no quantitative frameworks define thresholds, dynamics, or lead time.
- **Novelty:** state-space / mixed-effects longitudinal models of EV biomarker trajectories; data-driven definition of MRD positivity thresholds and recurrence lead-time estimation; digital-scoring-aware measurement-error models building on [253].
- **Data/tools:** published longitudinal EV cohorts [133, 246, 258]; ctDNA-MRD modeling literature as methodological template; Stan/pyMC or lme4.
- **First experiments:** (i) re-analyze a public longitudinal EV PD-L1 cohort [133] with a state-space model; (ii) derive threshold rules and compare against static cutoffs; (iii) simulate monitoring schedules to optimize sampling frequency.
- **Risks:** public longitudinal EV datasets are few → start from [133]/[246] supplements; simulation-based validation.

### Topic 10 — Data-driven prioritization of EV surface targets for capture and therapy
- **Rationale:** capture reagents target a small, historically-fixed marker set (CD9/CD63/CD81, EpCAM, PD-L1) [20, 29, 130], while Vesiclepedia 2024 holds 566,911 protein entries [119] and engineered-EV targeting uses empirical ligand choices [105, 106]; EV proteomes diverge from cells of origin [112].
- **Novelty:** a reproducible target-prioritization pipeline integrating EV proteomics (Vesiclepedia), surfaceome annotations, tumor vs. healthy tissue expression (TCGA/GTEx), and existing capture-reagent performance — ranking pan-cancer and cancer-specific EV surface targets for capture (diagnostics) and for engineered-EV homing (therapy).
- **Data/tools:** Vesiclepedia 2024 [119]; TCGA/GTEx; Human Protein Atlas; UniProt surfaceome; network propagation.
- **First experiments:** (i) build the integrated target table; (ii) benchmark the ranking against markers used in published devices [78, 211, 218, 220]; (iii) nominate 5 novel capture targets with supporting evidence profiles.
- **Risks:** EV protein detection bias toward abundant proteins → incorporate detectability modeling; purely computational.

**Portfolio cross-check against constraints:** all ten topics use public data and open tools; none requires wet-lab access for the core contribution; each has a 3–6-month first experiment yielding a concrete artifact (benchmark, model, or ranked list); each cites a specific, evidenced gap rather than a generic "apply ML" claim.

---

## 6. Limitations of this review

1. **Evidence depth is asymmetric:** 13 papers were read in full text; the remaining 245 are represented by abstracts/metadata, so quantitative claims beyond the full-text set should be verified against source before being built upon.
2. **Retrieval boundaries:** the 15 queries define the corpus; adjacent areas (e.g., EV-based vaccines, non-cancer EV diagnostics) were out of scope. One pivotal therapeutic-removal paper [164] had no accessible full text and is cited from abstract.
3. **Publication and language bias:** the corpus over-represents positive device/diagnostic performance; negative or null capture/diagnostic studies are rarely published.
4. **Rapid field movement:** several 2025–2026 records are early-access; performance claims (e.g., [49, 87]) await independent replication.
5. **The topic portfolio is a hypothesis-generating artifact:** feasibility judgments reflect public-data availability as of 2026-08; dataset licensing and access should be confirmed at project start.
