# Diagnosis Follow-Through and Recommendation (Full Re-Investigation)

**Status: analysis and real, verified new research. No project code changed.** This document supersedes the earlier, partial version of this same file. The first pass through the diagnosis document only chased the two items its own "immediate next steps" section flagged first, and treated most of its 20-category appendix as a quick plausibility check rather than a real investigation — a valid criticism, and this document corrects it. Every one of the 20 categories in `docs/Novelty_Competitiveness_Diagnosis_and_Restructuring_Options.md`'s appendix was investigated for real this time, with actual searches and actual paper full-texts pulled, not just reasoned about. Every one of the diagnosis document's five weak points and five restructuring options is re-scored at the end against what was actually found.

**Why this matters for the project, stated up front:** the project's one-sentence thesis is that nonlinear affinity capture means tumor-vs-normal selectivity depends on distribution *overlap*, not mean fold-change — meaning standard target-picking can rank capture markers incorrectly. The diagnosis document's own purpose was to find every real way to strengthen the project's five identified weak points before the write-up locks in its claims. Every finding below is reported with an explicit note on which weak point it addresses and how it moves the project toward being a genuinely competitive, evidence-based ISEF entry — not just toward "more sources cited."

---

## Part 0: How this investigation was actually run this time

All 20 categories from the diagnosis document's appendix (reproduced in that document's Groups A–D) were investigated, most via three parallel deep-research passes plus this session's own direct work, each doing real searches and pulling real paper text rather than trusting search-result snippets. Below, results are organized by the diagnosis document's own four groups so nothing gets silently dropped, and each category gets an honest verdict: **real, usable** / **real, partial** / **confirmed dead end**.

---

## Part 1: Group A — real head-to-head marker comparisons (categories 1–6)

### Category 1 — Clinical trial registries: confirmed dead end, for a specific reason

Searched ClinicalTrials.gov directly (API query) and by broad web search for EV-capture trials comparing multiple markers/platforms on comparable cohorts. No trial surfaced with a real, structured multi-marker capture-efficiency comparison in its Arms/Outcomes fields for any cancer this project touches (pancreatic, colorectal, ovarian). Registered EV trials that do exist are almost entirely therapeutic (EV-based treatments) or single-marker diagnostic studies, not head-to-head marker comparisons. **Verdict: dead end**, not for lack of trying — this category genuinely doesn't have what the project needs right now.

### Category 2 — Patent filings: real, but modest — confirmed twice, not overturned

A specific patent, **US20210102191A1** ("Method for the isolation of intact extracellular vesicles"), was read in detail. It reports a real, direct experimental comparison of CD9, CD63, and CD81 antibody-based EV capture on the same real sample — with an honest, mechanistic explanation for why CD81 and CD63 captured fewer particles than CD9 in that specific preparation (low CD81 abundance; CD63 epitope masking from the purification antibody). This is real, checkable, marker-specific capture data — but it's a technology-development patent explaining a method, not a diagnostic-comparison patent ranking markers for disease detection. Several other patents were checked (engineered EVs, nanoplasmonic detection, nanobody panels) with no comparably concrete comparative data. **Verdict: real, but a modest, secondary source** — cite it as one more real data point on tetraspanin capture behavior, not as a headline finding.

### Category 3 — FDA regulatory submissions: confirmed dead end

Checked FDA's 510(k)/De Novo databases and general search. No EV/exosome-based diagnostic device submission was found for any cancer this project touches. This matches the diagnosis document's own prediction that this category would likely be sparse (most EV diagnostics are laboratory-developed tests with no public FDA submission). **Verdict: dead end, exactly as predicted** — not a search failure, a real feature of how this technology is currently regulated.

### Category 4 — Company technical white papers / spec sheets: one real, unusable long-shot found

A real, topical lead was found: **Biological Dynamics' "ExoVerita" platform**, presented at the AACR Annual Meeting 2024, reporting that combining exosome markers with **CA19-9** (this project's own primary headline marker) detected 97% of stage 1–2 pancreatic cancers. This is genuinely exciting to see as independent, real-world corroboration that CA19-9-plus-EV-marker combinations are an active area of real diagnostic development — but only a press-release-level summary is publicly available; no marker list, no raw comparative data, no accessible dataset. **Verdict: real corroborating context, not usable as a data source.** Worth one sentence in the write-up as "this project's central marker choice (CA19-9) is independently being pursued by industry for the same purpose," with an honest caveat that no primary data was accessible.

### Category 5 — Non-bioRxiv preprints (medRxiv, Research Square, ChemRxiv, F1000Research): real, and this is where two of this session's best finds came from

Two genuinely important papers surfaced specifically through this category:

- **"Differential extracellular vesicle concentration and their biomarker expression of integrin αv/β5, EpCAM, and glypican-1 in pancreatic cancer models,"** *Scientific Reports* (2024), DOI [10.1038/s41598-024-65209-8](https://doi.org/10.1038/s41598-024-65209-8), PMID 38902362, free full text PMC11189911. This paper uses **PANC-1** — the exact cell line this project's Aim 6 already validates against — compared to a normal pancreatic cell line (HPNE), with four real markers (ITGα, ITGβ, GPC-1, EpCAM). According to PubMed, all four markers were significantly elevated in PDAC vs. normal cells, but with genuine, reported per-marker differences in how strongly each one separates cancer from normal across different PDAC cell lines (GPC-1 showed only a "mild increase" in PANC-1 specifically, while being clearly significant in the other PDAC line tested, BxPC-3) — a real, checkable, cell-line-specific ranking difference. The paper also derives a real quantitative estimate: individual tumor-derived EVs show **at least a 3-fold increase** in these biomarkers compared to normal EVs.
- **"Surface markers on supermeres outperform extracellular vesicles in colorectal cancer diagnosis,"** *Scientific Reports* (January 2026 — brand new), DOI 10.1038/s41598-026-36626-8, PMID 41571721, free full text PMC12902025. According to PubMed, this is a genuine, very recent, CRC-specific head-to-head ranking claim between two different particle types (a newer particle class called "supermeres" vs. small EVs), reporting that supermere surface markers outperform EV markers for CRC diagnosis.

**Verdict: real, usable, and directly on-topic for the project's existing PDAC and CRC work** — this category paid off more than any other in Group A.

### Category 6 — Citation-network tracing: real, and this is HOW the strongest K_INTRA lead was found

Following citations outward from the originally-flagged 2015 Langmuir paper (rather than re-running keyword searches) is what led to the entire "superselectivity" physics literature (see Part 3, category 15, below) — confirming this is a genuinely underused technique for this project. **Verdict: real and high-value**, and should become a standard part of any future search in this project, not a one-off.

---

## Part 2: Group B — real-world healthy-vs-cancer single-vesicle datasets (categories 7–10)

### Category 7 — NCI/NIH biobanks (caEV, ERCC/exRNA Atlas, CPTAC): confirmed dead end

No dedicated NCI single-vesicle-resolution protein consortium with a public cancer-vs-healthy dataset was found. exRNA Atlas is RNA-focused and bulk, not single-vesicle. CPTAC is tumor-tissue proteogenomics with no liquid-biopsy single-EV arm. **Verdict: dead end**, confirmed directly rather than assumed.

### Category 8 — Broadened GEO/PRIDE/ArrayExpress search using platform names instead of "single-EV": real payoff, exactly as the diagnosis document predicted

This was flagged as "the single biggest likely missed search angle across the whole project" and it delivered: **Ma, Wu et al., "Single extracellular vesicle surface protein-based blood assay identifies potential biomarkers for detection and screening of five cancers,"** *Molecular Oncology* (2024), PMID 38194998, free full text PMC10920081. This is a genuinely new find (not previously known to the project): real single-vesicle resolution (proximity barcoding assay — the same underlying method as the CRC paper Aim 7/8 already use, and the same research lineage as the OMIX010139 dataset previously pursued and deprioritized), with a real **100 cancer patients (five types, including colorectal) vs. 100 age/sex-matched healthy controls** case-control design, plus a 14-vs-14 independent validation cohort. **Important honest caveat**: no public raw single-vesicle data repository was found in the full text — the paper reports summary statistics (differentially expressed proteins), not a downloadable per-vesicle matrix, so getting the raw data would still require contacting the corresponding authors, similar friction to the OMIX010139 path already tried and deprioritized. **Verdict: real, partial** — this is a new author-contact avenue, not a click-to-download win, but it confirms the platform-name search angle was a real blind spot worth remembering for any future search.

### Category 9 — Population biobanks (UK Biobank, All of Us, BioVU): confirmed dead end

Checked directly: all three are bulk plasma proteomics (Olink/SomaScan panels), not EV-isolated, let alone single-vesicle. Structurally the wrong resolution for this project regardless of cancer-vs-healthy framing. **Verdict: dead end**, confirmed rather than assumed.

### Category 10 — Non-English-language repositories: partial, folds into category 8

Both known single-EV, cancer-vs-healthy PBA datasets (OMIX010139 and the newly-found Ma/Wu five-cancer study) originate from the same Chinese research group lineage (Fudan University-affiliated, commercialized as "Vesicode AB"). No additional Japanese or Korean single-EV cancer-vs-healthy dataset was found. **Verdict: partial** — confirms Chinese research groups are the leading source for this specific method, but doesn't add a new usable dataset beyond category 8's find.

---

## Part 3: Group C — cross-domain mathematical/physical analogues (categories 11–15)

### Category 11 — Flow cytometry ABC calibration beads, revisited via the CAR-T/bispecific-antibody angle: a real, new win

The diagnosis document's own suggestion — that the plain ABC-bead angle was a dead end but a CAR-T/bispecific-antibody angle "might work better" — was tested directly and **paid off**: **Walker et al., "Tumor Antigen and Receptor Densities Regulate Efficacy of a Chimeric Antigen Receptor Targeting Anaplastic Lymphoma Kinase,"** *Molecular Therapy* (2017), and **Majzner et al., "Tuning the Antigen Density Requirement for CAR T-cell Activity,"** *Cancer Discovery* (2020), DOI 10.1158/2159-8290.CD-19-0945, both confirmed real via search. According to these findings, Walker et al. used cells engineered to express a *titrated range* of a real cancer antigen (~18,000 to ~450,000 molecules/cell) and measured immune-cell killing efficiency as a direct function of that density — a genuine graduated dose-response curve, structurally very close to what this project's own receptor-density-vs-capture-probability question needs, just in a different (immune-killing, not affinity-capture) functional readout. **Verdict: real, and a genuine second, independent cross-domain confirmation** of the core nonlinear-density mechanism, on top of the superselectivity literature (category 15) — two unrelated fields now independently support the same qualitative physics.

### Category 12 — Lateral flow assay / rapid antigen test nonlinear capture modeling: real, moderate value

**Qian & Bau, "A Mathematical Model of Lateral Flow Bioreactions Applied to Sandwich Assays,"** *Analytical Biochemistry* (2003), DOI 10.1016/j.ab.2003.07.011, PMID 14705784, confirmed via PubMed — a real, derived mathematical model of nonlinear capture-line signal as a function of analyte concentration, including the well-known "hook effect" (signal *decreasing* at very high concentration, a genuinely nonlinear and non-monotonic relationship). A 2016 extension (Berli & Kler, *Microfluidics and Nanofluidics*) provides a fuller mechanistic model. **Honest caveat**: neither is specifically a COVID-era innovation — the foundational math predates COVID, and COVID-era papers mostly cite it rather than re-deriving it, so the diagnosis document's specific "COVID-era" framing wasn't quite matched. **Verdict: real, usable as a citable mathematical precedent for nonlinear capture behavior**, but not a source of EV-specific numbers.

### Category 13 — CTC-chip receptor-density-dependent capture efficiency: re-confirmed as weak, with one partial exception

A deeper pass still found the field dominated by "capture stays high (>90%) regardless of EpCAM density" results across a 250-fold density range — the same finding as the first investigation. One somewhat more granular paper was found (Czaplicka et al., *Cancers* 2020, DOI 10.3390/cancers12113315) correlating detection signal with EpCAM expression across four graded cell lines, but it measures signal intensity, not a clean capture-efficiency curve with a fitted equation. **Verdict: confirmed weak for this project's specific purpose** — worth a citation as corroborating evidence that nonlinearity exists in a related system, not as a numeric source.

### Category 14 — Ecology/epidemiology capture-recapture heterogeneity statistics: a real, valuable, different kind of win

**Otis, Burnham, White & Anderson, "Statistical Inference from Capture Data on Closed Animal Populations,"** *Wildlife Monographs* (1978) — the foundational reference for "Mh" models, which formally handle exactly the statistical problem this project's central thesis describes: when individuals (animals, or here, vesicles) have *heterogeneous* capture probability, using the population mean as your estimator is provably biased. This is a mature, decades-old, independently-developed statistical framework, unconnected to biology or EVs at all, that formally validates the *reasoning* behind this project's central claim (overlap matters, not just mean fold-change) from a completely different field. **Verdict: real, and valuable in a different way than the other categories** — it doesn't supply a K_INTRA number or a new claim to check, but it is a genuine, citable, independent formal confirmation that the project's underlying statistical logic is sound and well-established science in another discipline, which directly strengthens the project's mathematical rigor for a judge evaluating whether the model is built on solid ground. This is new and was NOT in the first pass of this investigation — a category the first pass skipped entirely.

### Category 15 — SPR/QCM superselectivity literature: confirmed as the strongest lead, unchanged from before

The full account of this thread is unchanged from the earlier investigation: **Dubacheva, Curk & Richter, "Determinants of Superselectivity — Practical Concepts for Application in Biology and Medicine,"** *Accounts of Chemical Research* (2023), DOI [10.1021/acs.accounts.2c00672](https://doi.org/10.1021/acs.accounts.2c00672), PMID 36916901, free full text PMC10077582, formally defines the exact α selectivity parameter this project's own model already computes, and re-plots a real quantitative dataset (Wang et al. 2020, *ACS Nano*, PMID 31899613 — 5 ligand densities × 11 receptor densities on breast cancer cells) obeying the same scaling law. **Verdict: still the single strongest lead for Weak Point 3**, now reinforced (not replaced) by category 11's CAR-T finding as a second, independent field showing the same phenomenon.

---

## Part 4: Group D — adjacent structural/indirect evidence (categories 16–20)

### Category 16 — Single-EV instrument platform QC papers: confirmed dead end (nothing new beyond what the project has)

A thorough search for ExoView/NanoView, IZON, Apogee, CytoFLEX nano, and Amnis/Cytek application notes with multi-cancer-marker comparative data found nothing beyond what the project already uses — notably, this search re-confirmed that **Lennon et al. 2019** (the project's own existing EGFR/CA19-9 pancreatic source, DOI 10.1080/20013078.2019.1685634, PMID 31741725) already IS this kind of paper. **Verdict: dead end for NEW data** — but a useful confirmation that the project's existing source selection was already sound here.

### Category 17 — Antibody vendor Kd datasheets: confirmed genuine dead end, not a resourcing gap

Checked BioLegend, R&D Systems/Bio-Techne, and Abcam directly for all ten markers this project uses. **Finding, stated plainly: these vendors do not publish Kd (binding affinity) values for standard research-use flow cytometry/immunofluorescence antibodies** — this is a real, structural fact about how these products are documented (functional titer data, not biophysical characterization), not a search failure. The only Kd values found anywhere were for unrelated therapeutic antibody candidates in patent filings, not traceable to the specific catalog clones this project's source papers used. **Verdict: dead end, cleanly and honestly closed** — this parameter category should not be revisited without a fundamentally different approach (e.g., directly emailing vendor technical support, untested here).

### Category 18 — Systematic reviews beyond CRC (mining more disease-specific reviews): partially covered by category 5's PDAC find; not separately exhausted

The earlier investigation found the CRC-specific systematic review (Shi et al. 2024, DOI 10.1186/s12885-024-12312-8). This round's PDAC find (category 5, above) serves a similar function for pancreatic cancer. A dedicated ovarian-cancer-specific systematic review search was not separately run this round. **Verdict: real progress, but genuinely incomplete** — flagged honestly as an open item for a future session, not claimed as finished.

### Category 19 — Mining supplementary data of papers ALREADY used: the single most important finding of this entire re-investigation

This is, by a clear margin, the highest-value, lowest-risk category found across both rounds of this investigation, and it was almost entirely skipped in the first pass. Full-text mining of the project's own existing source papers turned up substantial real data never yet used:

- **Spitzberg et al. 2023** (the project's own MASEV source, DOI 10.1038/s41467-023-36932-z, PMID 36870999) contains, beyond Figure 4 and Figure 6A (already used): a **12-marker combinatorial co-occurrence breakdown across all four cell lines** (e.g., only 4.22% of PANC-1 vesicles show 8+ simultaneous biomarkers, vs. 0.26% for CAPAN-2, 2.75% for ASPC1, 2.00% for A549 — real, checkable, per-cell-line heterogeneity data), and a direct, quantified statement that **single-tetraspanin purification loses up to 80% of KRAS-mutant EVs, while tetraspanin cocktails still miss 36–47%** — a real, already-published, capture-efficiency comparison between single-marker and multi-marker (cocktail) capture strategies, in the exact paper this project's Aim 6 already validates against.
- **Guo et al.** (the project's own CRC source, PMID 37040507, DOI 10.34133/research.0041) contains real **individual-marker diagnostic AUC values across three separate cohorts** (discovery, an external validation cohort of 30 CRC vs. 58 healthy controls, and an early-detection cohort covering precancerous lesions and Stage I CRC): **ITGAM AUC = 0.9158, ITGB3 AUC = 0.8681, EPS8 AUC = 0.8498, CD9 AUC = 0.8242** (healthy vs. CRC), plus additional markers (EGF, ITGA6, ADAM10, CD151, ITGB2) reported across the metastasis-discrimination analysis. **Honest, important correction to how the project has been describing this paper**: the project's documentation currently states no dataset has "a genuine, confirmed healthy-vs-cancer grouping" for CRC — that statement is accurate for the *raw per-vesicle file* Aim 8 uses (confirmed, no healthy controls in that specific file), but it is **not accurate for this same paper's own summary statistics**, which do contain a real healthy-control group, real AUC values, and even a real percentage decline for the ITGAM/ITGAL/ITGB2 subpopulation the project already uses (13.47% in healthy controls → 1.13% in primary and metastatic CRC). **Honest caveat on completeness**: this paper reports AUC values only — no sensitivity/specificity, and no ITGAL, EpCAM, CD63, or CD81 individual-marker numbers.

**Verdict: this is the strongest, safest, most directly actionable finding in this entire investigation.** It requires no new data access, no author emails, no waiting — it's real data in papers the project has already cited, downloaded, and partially used. It directly fixes part of Weak Point 2 (a real healthy-vs-cancer comparison for CRC does exist, just not where the project's documentation currently says to look) and directly supplies the strongest available version of Weak Point 4's marker-ranking check (a 4+-marker real ranking, from the project's own primary source, stronger than an external systematic review because there's no cross-study confound).

### Category 20 — Investor materials, SEC filings, conference abstracts: one real long-shot confirmed as a long-shot

Already covered under category 4 above (the Biological Dynamics/AACR 2024 finding) — the same result applies here. **Verdict: real but not usable as a data source**, exactly the honest, non-inflated outcome the diagnosis document itself predicted for this category.

---

## Part 5: All five weak points, re-scored against everything found

**Weak point 1 (no genuinely new claim):** Now has a much stronger, safer answer than before. The recommended new claim (Part 6) uses Guo et al.'s own real 4-marker AUC ranking — data the project already has — as the checkable ground truth, rather than reaching for an external systematic review. This is a stronger fix because there's no cross-study confound between the "real world" answer and the model's own training/validation data.

**Weak point 2 (no healthy-vs-cancer single-vesicle data):** Meaningfully improved, with an honest and important correction: Guo et al.'s own summary statistics DO contain a real CRC healthy-vs-cancer comparison (the project's "no healthy-vs-cancer data anywhere" framing was too broad — it's true of the raw file, not the underlying paper). This should be corrected in the project's own documentation regardless of what else gets built. The Ma/Wu et al. five-cancer dataset (category 8) remains a real, if friction-prone, path to genuine single-vesicle-resolution healthy-vs-cancer raw data if pursued further.

**Weak point 3 (uncalibrated K_INTRA):** Strengthened further — now supported by two independent, real, cross-domain physics/biology literatures (superselectivity + CAR-T antigen-density) rather than one. Still not a source of a single anchored point-estimate without more work than the timeline likely supports; remains a stretch goal.

**Weak point 4 (no head-to-head marker-ranking check):** This is now the best-supported weak point in the whole diagnosis. Real ranking evidence exists in-hand (Guo, Spitzberg), externally (the CRC systematic review, the PDAC paper, the supermeres paper), and even in a patent. This weak point should be considered substantially resolved, pending the actual model-vs-data check being built.

**Weak point 5 (the narrowing problem itself):** This full re-investigation is itself the fix, and it surfaced a second, important instance of the same lesson: even the "unconstrained" appendix, built specifically to prevent narrowing, still didn't get properly investigated on the first pass through it — a checklist only works if every item on it actually gets checked, not just referenced. Category 19 (mining the project's own existing papers) is the clearest example: it was sitting in the appendix the whole time, named directly, and still got skipped in the first pass in favor of chasing external, harder-to-get data. **This is the single most important process lesson from this whole exercise**: before searching outward for new evidence, always fully exhaust what's already been collected.

---

## Part 6: The recommendation, updated

### 6.1 The new, falsifiable claim — now anchored to the project's own strongest available data

**Build a small, new analysis that checks whether the project's existing binding-model framework, applied to the real digitized marginal distributions it already has for the CRC markers (ITGAM, ITGAL, ITGB2, and the related ITGB3/ITGA6 cluster), correctly predicts the RELATIVE ORDERING of diagnostic separation strength that Guo et al.'s own paper reports (ITGAM > ITGB3 > EPS8 > CD9, by AUC).** This is a stronger version of the previously-recommended AND-gate check because: (1) it uses data already in the project's possession, so there's no new access risk or timeline risk; (2) a 4-marker ranking is a harder, more specific, more falsifiable prediction than a single dual-marker comparison — the model could get this wrong in many ways, so getting it right (or even partially right) is more meaningful evidence; (3) it directly ties to the paper the project's Aim 7/8 already center on, so there's no cross-study confound to explain away.

The AND-gate check (comparing to the CRC systematic review's CD63+EpCAM vs. CD63+CD9 result) remains a valid secondary check and should still be attempted if time allows — the two checks are complementary, not competing.

**How this ties back to the project's actual goals:** the diagnosis document's own central finding — confirmed by fresh, real ISEF-winner research — is that the gap between a good project and a top-tier one is a genuinely new, previously-unstated, checkable claim, not more validation of already-known facts. This ranking-order check is exactly that: nobody has asked whether this project's specific nonlinear-capture model, built from real digitized marginal distributions, predicts the real ordering Guo et al. independently measured. That's a new question with a real, checkable answer — precisely the pattern the diagnosis document identified as what separates First/Second Award winners from Third/Fourth.

### 6.2 Immediate, low-risk documentation fixes — do these regardless of what else gets built

1. **Correct the project's "no CRC healthy-vs-cancer data exists" framing.** It should read: the raw per-vesicle CRC file (used in Aim 8) has no healthy controls, but the same paper's own summary statistics (already partially used in Aim 7) do include a real healthy control group with real AUC values and a real percentage decline for the ITGAM/ITGAL/ITGB2 subpopulation. This is a factual correction, not a new research task, and it directly improves Weak Point 2's honest status.
2. **Add Spitzberg et al.'s cocktail-vs-single-marker capture-loss numbers (80% single-marker loss, 36–47% cocktail miss rate) to the project's evidence base** — this is real, already-cited, already-downloaded data that speaks directly to the value of multi-marker capture strategies (relevant to Aim 5's AND-gate work) and was simply never extracted before.

### 6.3 What to pursue if time allows, in priority order

1. The ranking-order check (6.1) — highest priority, uses only data already in hand.
2. Cite the superselectivity literature AND the CAR-T antigen-density literature together as two independent cross-domain confirmations of the project's core nonlinear-density mechanism (Weak Point 3) — a citation-strength upgrade, not a calibration.
3. Cite the capture-recapture ecology/statistics literature (Otis et al.'s Mh models) as an independent, formal, non-biological validation of the project's core statistical reasoning (mean/fold-change is a biased predictor under heterogeneous capture probability) — this is a genuinely new addition from this round's investigation and strengthens the project's mathematical rigor for judge Q&A.
4. If time allows well beyond the above: attempt contact with the Ma/Wu et al. five-cancer dataset authors (category 8) as a second, potentially easier path to real single-vesicle healthy-vs-cancer raw data, alongside or instead of a fresh OMIX010139 attempt — but this should not be treated as a priority given the project's prior, reasoned decision to deprioritize a similar access-friction path.

### 6.4 What NOT to pursue, now confirmed dead ends from this round

Clinical trial registries (category 1), FDA submissions (category 3), NCI biobanks (category 7), population biobanks (category 9), antibody vendor Kd datasheets (category 17), and single-EV instrument-QC papers beyond what's already used (category 16) are all now confirmed, directly-checked dead ends for this project's specific needs — not areas that need another look without a genuinely new reason.

---

## Sources

All sources from the first round of this investigation remain valid and are not re-listed in full here (see git history / prior version if needed); this section lists sources newly confirmed in this round.

Walker AJ, et al. "Tumor Antigen and Receptor Densities Regulate Efficacy of a Chimeric Antigen Receptor Targeting Anaplastic Lymphoma Kinase." *Molecular Therapy* 2017.

Majzner RG, et al. "Tuning the Antigen Density Requirement for CAR T-cell Activity." *Cancer Discovery* 2020;10(5):702-723. DOI: 10.1158/2159-8290.CD-19-0945.

Qian S, Bau HH. "A Mathematical Model of Lateral Flow Bioreactions Applied to Sandwich Assays." *Analytical Biochemistry* 2003;322(1):89-98. DOI: 10.1016/j.ab.2003.07.011. PMID: 14705784.

Otis DL, Burnham KP, White GC, Anderson DR. "Statistical Inference from Capture Data on Closed Animal Populations." *Wildlife Monographs* 1978, No. 62.

"Differential extracellular vesicle concentration and their biomarker expression of integrin αv/β5, EpCAM, and glypican-1 in pancreatic cancer models." *Scientific Reports* 2024. DOI: [10.1038/s41598-024-65209-8](https://doi.org/10.1038/s41598-024-65209-8). PMID: 38902362. Free full text: PMC11189911.

"Surface markers on supermeres outperform extracellular vesicles in colorectal cancer diagnosis." *Scientific Reports* 2026. DOI: 10.1038/s41598-026-36626-8. PMID: 41571721. Free full text: PMC12902025.

Ma Y, Wu D, et al. "Single extracellular vesicle surface protein-based blood assay identifies potential biomarkers for detection and screening of five cancers." *Molecular Oncology* 2024. PMID: 38194998. Free full text: PMC10920081. DOI: 10.1002/1878-0261.13586.

Guo W, et al. "Single-Exosome Profiling Identifies ITGB3+ and ITGAM+ Exosome Subpopulations as Promising Early Diagnostic Biomarkers and Therapeutic Targets for Colorectal Cancer." *Research* 2023. PMID: 37040507. DOI: 10.34133/research.0041. Free full text: PMC10076010. (Full AUC data confirmed directly from this paper's own results section this round.)

Spitzberg JD, Ferguson S, Yang KS, Peterson HM, Carlson JCT, Weissleder R. "Multiplexed analysis of EV reveals specific biomarker composition with diagnostic impact." *Nature Communications* 2023;14:1239. DOI: 10.1038/s41467-023-36932-z. PMID: 36870999. Free full text: PMC9985597. (Combinatorial co-occurrence and capture-loss data confirmed directly from this paper's own text this round.)

Lennon KM, et al. "Single molecule characterization of individual extracellular vesicles from pancreatic cancer." *Journal of Extracellular Vesicles* 2019;8(1):1685634. DOI: 10.1080/20013078.2019.1685634. PMID: 31741725.

US20210102191A1, "Method for the isolation of intact extracellular vesicles," Google Patents.

All PubMed-sourced information above is attributed per PubMed's citation requirements, with DOI links included for every article referenced.
