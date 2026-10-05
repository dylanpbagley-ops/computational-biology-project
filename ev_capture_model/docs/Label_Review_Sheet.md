# Superseded review sheet — use the October 2 score-free form

This historical sheet exposes the original labels to reviewers and contains known source errors. Use `strengthening_Oct2/reviewer_form_score_free.tsv` together with `Reviewer_Instructions.md` in that directory. Do not send the review key or model scores to the reviewer. The historical content below is preserved for provenance, not current clinical adjudication.

# Label review sheet — whole-body target-safety analysis

**Purpose.** Every result in section 23 of the results document depends on the outcome column
below. Those labels are literature judgements made by one person (the student), and the single
cheapest way to strengthen the work is for somebody else to check them. This sheet exists to
be checked.

**How to check it.** For each row, the question is only: *is the outcome label right?* Nothing
about gene expression, cell types or the model is involved. A reviewer needs no computational
background — only the clinical record.

**Rules the labels follow (fixed before any number was computed):**

- `safe` — the agent reached approval, or is documented as clinically tolerable, **and** no
  dose-limiting on-target toxicity in a normal tissue was reported.
- `toxic` — a **severe, dose-limiting or fatal** on-target toxicity in a normal tissue is
  documented. Note that a target can be `toxic` under this rule and still be an approved drug
  (FOLR1/mirvetuximab carries a boxed warning for on-target ocular toxicity).
- `efficacy_fail` — the programme failed, but for lack of efficacy rather than toxicity. These
  rows are **excluded from the toxicity test** and used only as a specificity control.
- `excluded_*` — stated exclusions, with the reason given in the row.

**Reviewer:** ______________________  **Date:** ____________
**Rows I disagree with (row number + what it should be):**

<br>

| # | gene | target | agent | outcome | evidence as recorded | expected normal site | agree? |
|---|---|---|---|---|---|---|---|
| 1 | `CD19` | CD19 | tisagenlecleucel / axi-cel | **safe** | Approved CAR-T; on-target B-cell aplasia is managed with immunoglobulin replacement | B cells | &#9744; |
| 2 | `MS4A1` | CD20 | rituximab | **safe** | Approved 1997; B-cell depletion clinically managed | B cells | &#9744; |
| 3 | `CD22` | CD22 | inotuzumab ozogamicin | **safe** | FDA approved 2017 (ALL); B-cell depletion is on-target and clinically managed | B cells | &#9744; |
| 4 | `CD79B` | CD79b | polatuzumab vedotin | **safe** | FDA approved 2019 (DLBCL); no dose-limiting normal-tissue on-target toxicity | B cells | &#9744; |
| 5 | `TNFRSF17` | BCMA | ide-cel / cilta-cel | **safe** | Approved 2021; plasma-cell aplasia managed | plasma cells | &#9744; |
| 6 | `SLAMF7` | CS1 | elotuzumab | **safe** | FDA approved 2015 (myeloma); well tolerated in combination regimens | plasma cells / NK | &#9744; |
| 7 | `CD38` | CD38 | daratumumab | **safe** | FDA approved 2015 (myeloma); on-target plasma-cell loss is clinically managed | plasma cells | &#9744; |
| 8 | `TNFRSF8` | CD30 | brentuximab vedotin | **safe** | FDA approved 2011 (Hodgkin, ALCL); neuropathy is payload-driven, not on-target | activated lymphocytes | &#9744; |
| 9 | `CD33` | CD33 | gemtuzumab ozogamicin | **safe** | Approved 2000/2017; myelosuppression is on-target but marrow regenerates | myeloid | &#9744; |
| 10 | `GPRC5D` | GPRC5D | talquetamab | **safe** | Approved 2023; on-target skin/nail/taste effects are common but not dose-limiting | plasma cells, keratinised tissue | &#9744; |
| 11 | `FOLH1` | PSMA | 177Lu-PSMA-617 | **safe** | FDA approved 2022 (mCRPC radioligand); salivary uptake is on-target but not dose-limiting | prostate epithelium | &#9744; |
| 12 | `TACSTD2` | TROP2 | sacituzumab govitecan | **safe** | FDA approved 2020 (triple-negative breast); toxicity is payload-driven neutropenia | epithelium | &#9744; |
| 13 | `NECTIN4` | Nectin-4 | enfortumab vedotin | **safe** | Approved 2019 (skin effects on-target but not dose-limiting) | urothelium/skin | &#9744; |
| 14 | `CLDN18` | claudin-18.2 | zolbetuximab | **safe** | Approved 2024; nausea/vomiting on-target on gastric mucosa, manageable | gastric mucosa | &#9744; |
| 15 | `DLL3` | DLL3 | tarlatamab | **safe** | FDA approved 2024 (small-cell lung); DLL3 is essentially absent from adult normal tissue | neuroendocrine | &#9744; |
| 16 | `MET` | c-MET | telisotuzumab vedotin | **safe** | FDA approved 2025 (c-MET-overexpressing non-squamous NSCLC); no dose-limiting on-target normal toxicity | epithelium | &#9744; |
| 17 | `MSLN` | mesothelin | anti-MSLN CAR-T / anetumab | **safe** | Repeatedly tolerable in trials; serositis only at high intrapleural dose | mesothelium | &#9744; |
| 18 | `GPC3` | glypican-3 | GPC3 CAR-T | **safe** | Tolerable in hepatocellular trials | fetal liver (absent in adult) | &#9744; |
| 19 | `CD52` | CD52 | alemtuzumab | **safe** | Approved; profound but managed lymphodepletion | lymphocytes | &#9744; |
| 20 | `ERBB2` | HER2 | HER2 CAR-T (Morgan 2010) | **toxic** | Fatal pulmonary event within hours of infusion; trastuzumab also causes on-target cardiotoxicity | lung epithelium, cardiomyocytes | &#9744; |
| 21 | `EGFR` | EGFR | anti-EGFR CAR-T / cetuximab | **toxic** | Dose-limiting on-target skin and pulmonary toxicity | keratinocytes, lung | &#9744; |
| 22 | `CEACAM5` | CEA | CEA CAR-T (Parkhurst 2011) | **toxic** | Severe transient colitis, dose-limiting | colonocytes | &#9744; |
| 23 | `CA9` | CAIX | CAIX CAR-T (Lamers 2006) | **toxic** | Grade 2-4 on-target liver enzyme toxicity, trial halted | bile duct epithelium | &#9744; |
| 24 | `CD47` | CD47 | magrolimab | **toxic** | On-target anaemia; programme discontinued 2024 | erythrocytes, ubiquitous | &#9744; |
| 25 | `CD44` | CD44v6 | CD44v6 CAR-T / bivatuzumab | **toxic** | Fatal on-target skin toxicity with bivatuzumab mertansine | keratinocytes | &#9744; |
| 26 | `NCAM1` | CD56 | lorvotuzumab mertansine | **toxic** | Dose-limiting on-target peripheral neuropathy; programme halted | neurons | &#9744; |
| 27 | `EPCAM` | EpCAM | catumaxomab / adecatumumab | **toxic** | On-target hepatotoxicity; catumaxomab withdrawn | hepatocytes/biliary | &#9744; |
| 28 | `FOLR1` | FRalpha | mirvetuximab soravtansine | **toxic** | Dose-limiting on-target ocular keratopathy (boxed warning) despite approval | corneal epithelium | &#9744; |
| 29 | `FAP` | FAP | FAP CAR-T | **toxic** | On-target toxicity to FAP+ bone-marrow stromal cells; cachexia and marrow failure | marrow stroma | &#9744; |
| 29b | `F3` | tissue factor | tisotumab vedotin (Tivdak) | **toxic** | FDA approved 2021/2024 (cervical cancer); boxed warning for ocular toxicity, up to 60% of patients (conjunctivitis, dry eye, corneal damage, blepharitis), dose-limiting, requires dedicated eye-care REMS | corneal/conjunctival epithelium | &#9744; |
| 30 | `FCER2` | CD23 | lumiliximab | **efficacy_fail** | Phase 3 stopped for futility, not toxicity | B cells | &#9744; |
| 30b | `GPNMB` | GPNMB | glembatumumab vedotin (CDX-011) | **efficacy_fail** | METRIC phase 2b randomized trial (gpNMB-overexpressing triple-negative breast cancer) did not meet primary endpoint (PFS); no on-target toxicity signal reported | melanocytes / osteoclasts | &#9744; |
| 31 | `CD40` | CD40 | dacetuzumab | **efficacy_fail** | Stopped for lack of efficacy | B cells / APC | &#9744; |
| 32 | `CD80` | CD80 | galiximab | **efficacy_fail** | Stopped for lack of efficacy | B cells / APC | &#9744; |
| 33 | `HLA-DRA` | HLA-DR | apolizumab | **efficacy_fail** | Stopped for lack of efficacy | all APC | &#9744; |
| 34 | `SDC1` | CD138 | indatuximab ravtansine | **efficacy_fail** | Development discontinued, efficacy | plasma cells / epithelium | &#9744; |
| 35 | `MUC1` | MUC1 | several | **efficacy_fail** | No agent has reached approval; failures attributed to efficacy/shedding | glandular epithelium | &#9744; |
| 36 | `MUC16` | CA125 | oregovomab | **efficacy_fail** | Phase 3 negative for efficacy; no on-target normal-tissue toxicity signal reported | coelomic epithelium | &#9744; |
| 37 | `SLC34A2` | NaPi2b | lifastuzumab vedotin | **efficacy_fail** | Did not beat comparator. **Documentation note added September 10, 2026:** a related ADC on the same target, upifitamab rilsodotin, showed a bleeding-event signal in trials; it was never mechanistically confirmed as on-target, and the sponsor's own stated reason for discontinuation was efficacy, not safety — so the label is unchanged, but the wrinkle is recorded here rather than silently left out. | lung/thyroid | &#9744; |
| 38 | `CD74` | CD74 | milatuzumab | **efficacy_fail** | Development discontinued | B cells / APC | &#9744; |
| 38b | `ERBB3` | HER3 | patritumab deruxtecan | **efficacy_fail** | FDA BLA voluntarily withdrawn 2026 after the confirmatory phase 3 HERTHENA-Lung02 trial (EGFR-mutated NSCLC) missed its overall-survival endpoint despite an earlier positive PFS signal; sponsors (Daiichi Sankyo/Merck) attributed the withdrawal to the OS miss, not a safety finding. Two ILD-related deaths occurred in the treatment arm but were not cited as the reason for withdrawal — a clean efficacy failure by the sponsors' own account, same pattern as GPNMB. | epithelium (broadly expressed, low level) | &#9744; |
| 38c | `AXL` | AXL | enapotamab vedotin | **efficacy_fail** | Genmab discontinued development in 2020 after early-phase expansion cohorts; company statement cited not meeting the "stringent criteria for proof-of-concept," with "some evidence of clinical activity" that could not be improved by alternate dosing or biomarker selection. Explicitly an efficacy-driven decision, no toxicity concern cited. | epithelium / mesenchymal (broadly expressed) | &#9744; |
| 39 | `PTPRC` | CD45 | (none) | **excluded_control** | No cytotoxic agent ever reached patients; retained as a restriction negative control | all leukocytes | &#9744; |
| 40 | `IL3RA` | CD123 | tagraxofusp vs talacotuzumab | **excluded_ambiguous** | One approval and one toxicity failure on the same gene | myeloid/pDC | &#9744; |

## Counts

| outcome | n | used in the toxicity test? |
|---|---|---|
| safe | 19 | yes |
| toxic | 11 | yes |
| efficacy_fail | 12 | no — specificity control only |
| excluded | 2 | no |

**Added September 9, 2026 (rows 29b, 30b):** two targets found by re-checking the ADC trial
registry (De Schepper et al. 2025, Supplementary Data 1/2 — 72 genes) against this project's
40-target panel for genes with a settled clinical verdict not yet recorded. Both were verified
against primary/near-primary sources (FDA label for F3/tisotumab vedotin; the published METRIC
trial result for GPNMB/glembatumumab vedotin) before being added — this sheet is not the place
to record an unverified or preliminary-looking result. Every other gene checked in that pass
(roughly 50 candidates) had no settled outcome yet — still in active trials, which is exactly
why they appear in `prospective_predictions.json` as forward predictions rather than as labels
here. Growing this panel further means periodically re-checking those same genes for a result,
not inventing new ones.

**Added September 10, 2026 (rows 38b, 38c):** two more targets found the same way, re-checking
the remaining ~30 unresolved registry genes for a verdict that has since settled. ERBB3/patritumab
deruxtecan and AXL/enapotamab vedotin were both verified against primary/near-primary sources
(the FDA withdrawal announcement and HERTHENA-Lung02 trial reporting for ERBB3; Genmab's own 2020
discontinuation statement for AXL) before being added. This growth moved the specificity control
(section 23.4) from borderline non-significance (p = 0.043, 19 vs 10) to a small but statistically
real separation (p = 0.011, 19 vs 12) — recorded and discussed honestly in the results document
rather than left unmentioned. It does not change the safe-vs-toxic result this sheet exists to
check (sections 23.1–23.3, 23.5–23.7 are numerically unchanged). Also documented this round: a
wrinkle on the existing SLC34A2 row (see row 37) that was investigated and found not to warrant a
relabel, but is worth a reviewer's independent look.

## What a reviewer's disagreement would do

The leave-one-out analysis (section 23.3) drops each of the outcome-bearing targets in turn and
re-runs: as of the 44-target panel (30 safe+toxic targets), separation stays between 0.905 and
0.947 and remains significant in **30 of 30** drops. So a single relabelled row cannot overturn
the result — but a systematic pattern of relabelling could, and that is exactly what an outside
check is for. Rows 38b and 38c (ERBB3, AXL) are new since the last review pass and have not yet
been checked by an outside reviewer; the SLC34A2 wrinkle note on row 37 is likewise new.
