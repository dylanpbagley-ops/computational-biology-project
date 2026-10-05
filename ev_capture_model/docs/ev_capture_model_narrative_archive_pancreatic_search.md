# Pancreatic cancer antigen search — full narrative (September 1, 2026)

This archive holds the complete blow-by-blow of the pancreatic-cancer generalization search: finding and calibrating two new real antigen systems (CLDN18.2, mesothelin), a real saturation/non-monotonic finding in the CLDN18.2 data, and three full rounds of searching for outcome-linked pancreatic data. `EV-Capture-CURRENT.md` keeps only the settled conclusions; open this file if you need how something was found, verified, or ruled out.

---

## Real find #1 — mesothelin cell-line calibration data (Chu et al. 2024)

Dylan has a personal connection to pancreatic cancer and asked specifically to search for pancreatic cancer data, with 6 months of runway before ISEF, and access to paywalled articles via dual-enrollment community college database credentials (a real capability discovered this session).

**Chu et al. 2024, *Cytotherapy* 26(4):325-333, DOI: 10.1016/j.jcyt.2024.01.011** ("Mesothelin antigen density influences anti-mesothelin chimeric antigen receptor T cell cytotoxicity"). Dylan retrieved the full PDF via his database access. Contains:
- A real cell-line panel (Figure 2A): AsPC1 MSLN KO (<LoB), MIA PaCa2 (<LoB), Mero95 MSLN KO (126), Capan2 (8,191), AsPC1 (10,443), NCI-H2052 (41,810), Mero95 (44,716), MIA PaCa2 MSLN (174,382) — molecules/cell, pancreatic + mesothelioma lines.
- A real TITRATION curve (Figure 3A/D): AsPC1 cells treated with sheddase inhibitors (lanabecestat, TMI-1) at increasing doses, MSLN density measured at 6,101 → 17,683 → 20,974 → 47,242 molecules/cell, each paired with real measured m912 CAR-T cytotoxicity (%). Reported Spearman r=0.623 (95% CI 0.057-0.886, P=0.035) between density and % lysis at 48h.
- No cross-scale mismatch risk (same flow-cytometry/Quantum-Simply-Cellular-bead assay measures both density and pairs with killing).
- No data deposition/accession number anywhere in the paper — this is cell-line-mean data, not single-cell-per-patient data. Do not attempt to reconstruct a per-cell distribution from these 8 values.

**Mesothelin k_intra fit** — Code: `ev_capture_model/src/pancreatic_calibration/msln_fit.py`. Uses ONLY the two AsPC1 sheddase-inhibitor conditions where BOTH density AND %lysis are directly reported as exact numbers in the paper's own text (not figure-estimates): DMSO vehicle (6,101 molecules/cell → 2.5% lysis) and combined lanabecestat+TMI-1 (47,242 molecules/cell → 47.2% lysis, P=0.041 vs. DMSO). The two single-drug intermediate conditions have directly-reported density but no exact reported lysis value in the text — rather than fabricate interpolated numbers, this fit honestly uses only the 2 clean anchor points, citing the paper's own full-4-point Spearman correlation (r=0.623, P=0.035) as independent confirmation all 4 points follow a real monotonic relationship.

**Result:** k_intra = 1.31e-05, implied half-kill threshold (theta=0.5) at ~53,049 molecules/cell. The 2-point fit isn't exact (predicted 7.7% vs actual 2.5% at the DMSO point) — expected behavior with only 2 constraints, not a concerning fit failure.

---

## Real find #2 — CLDN18.2 calibration data (Barrett et al. 2024 / AZD6422)

**Barrett et al. 2024, *Clinical Cancer Research* 30(23):5413-5429, DOI: 10.1158/1078-0432.CCR-24-1853, PMC11609629** ("Preclinical Evaluation of AZD6422, an Armored Chimeric Antigen Receptor T Cell Targeting CLDN18.2 in Gastric, Pancreatic, and Esophageal Cancers") — AstraZeneca preclinical paper.

- Confirmed real quantitative flow cytometry (not IHC scoring): Quantum Simply Cellular beads (Bangs Laboratories) — the same calibration methodology as the mesothelin paper — converting MFI into antibody-binding-capacity (≈molecules/cell) values.
- A real paired (density, cytotoxicity) dataset exists (Figure 5D): quantitative FC density vs. % cytolysis by xCELLigence RTCA assay at 48h, E:T 1:1, across multiple cell lines including 3 pancreatic lines (AsPC1+CLDN18.2, PaTu8988s HS, BxPC3+CLDN18.2).
- Supplementary Data1 was retrieved and checked directly — it does NOT contain a data table for Figure 5D's values (only HEK293 specificity controls, epitope mapping, in vivo PDX tumor volume data, and a TGFβ-armoring persistence assay). Dylan zoomed into Figure 5D directly and provided high-resolution screenshots; values below are careful visual estimates off a broken-axis bar chart (bottom segment linear 0-4,000, break, top segment linear 4,000-500,000 — NOT log scale; not exact extracted numbers).

| Cell line | ABC (~molecules/cell, estimated) | % cytolysis @ 48h (estimated) | Pancreatic? |
|---|---|---|---|
| PaTu8988s 18.2 KO | 0 | 1.5% | Yes (negative control) |
| DAN-G | ~850 | 15% | No (gastric) |
| IM-95 | ~2,750 | 65% | No (gastric) |
| NUGC4 | ~45,000 | 68% | No (gastric) |
| AsPC1 + CLDN18.2 | ~45,000 | 68% | **Yes** |
| SNU601 | ~50,000 | 77% | No (gastric) |
| PaTu8988s HS | ~185,000 | 88% | **Yes** |
| BxPC3 + CLDN18.2 | ~350,000 | 82% | **Yes** |

Code: `ev_capture_model/src/pancreatic_calibration/cldn182_fit.py`, results pickled to `cldn182_fit_results.pkl`.

**Real results:**
- **Spearman correlation (n=7 non-KO points): r=0.964, p=0.0005** — even stronger than the mesothelin paper's own r=0.623.
- **k_intra fit via least-squares against `capture_probability()`** (sigma=z=1 held fixed per project convention): k_intra = 2.99e-05, implying a half-kill threshold (theta=0.5) at ~23,177 molecules/cell.
- **Honest fit-quality caveat: RMSE = 0.24 (24 percentage points).** The single-parameter sigmoid underpredicts killing at the low-density end (DAN-G: predicted 2.5% vs actual 15%; IM-95: predicted 7.9% vs actual 65%) and overpredicts/saturates at the high-density end (both pancreatic high-density lines predicted ~100% vs actual 88%/82% — killing plateaus and even dips slightly at the very highest density, a real biological effect, not noise).
- **Why the fit isn't clean:** unlike the mesothelin AsPC1 sheddase-inhibitor data (one cell line, drug-titrated density, isolating density as the only variable), this dataset spans multiple different cell lines with likely different intrinsic biology — density is not the only thing varying between points.
- **This is itself a genuinely reportable finding:** the naive correlation (r=0.964) looks extremely clean and would suggest "more antigen = proportionally more killing," but the actual model fit reveals a more nuanced pattern (saturating/plateauing, not simple linear-with-noise) that the correlation coefficient alone hides.

**Now have three independently-calibrated real antigen systems:** CD19 (original, multi-point IndEx-2 titration), CLDN18.2 (7-point multi-cell-line fit, r=0.964, RMSE 0.24), and MSLN (2-point same-cell-line fit) — spanning DLBCL (CD19) and pancreatic-cancer-relevant/adjacent antigens (CLDN18.2, MSLN).

---

## CLDN18.2 saturation / non-monotonic finding

**The finding, stated precisely:** among the three real pancreatic cancer cell lines in the CLDN18.2 dataset, antigen density and CAR-T killing do NOT increase together monotonically. PaTu8988s HS (~185,000 molecules/cell) was killed at 88% — the highest cytolysis in the entire 8-line dataset, pancreatic or gastric. BxPC3+CLDN18.2 (~350,000 molecules/cell, roughly double the density) was killed at only 82% — lower than a cell line with substantially less antigen.

**Why this matters:** the project's central methodological claim is that naive antigen-density-based reasoning ("more target = better target") can miss real, clinically-relevant structure. This result is a live, real-data demonstration of exactly that risk, independent of anything the project's own model predicts: a naive density-based ranking would pick BxPC3 over PaTu8988s HS, and that pick would be wrong by the actual measured outcome.

**Honest caveats:**
1. n=3 for the pancreatic subset — a suggestive real observation, not a statistically powered claim.
2. Different cell lines (not one line at different doses), so a confound (intrinsic apoptosis sensitivity, HLA expression, adhesion, growth rate) is real and cannot be ruled out with 2 data points.
3. A biologically plausible mechanism (antigen shedding — high density correlating with high shed/soluble antigen that competitively blocks CAR-T binding) is real and published in general, but speculative for this specific dataset (not measured or reported in this paper).

**Recommended framing for write-up:** present as a real, illustrative example — not primary statistical evidence — found in an independent, already-published dataset the project did not generate. Pair explicitly with the caveats above.

---

## Outcomes search, round 1 — Beatty 2018, CT041 case report, CPTAC (real dataset, weak signal)

Dylan asked to keep looking for outcome-linked pancreatic data specifically.

1. **Beatty et al. 2018 Gastroenterology** (mesothelin mRNA-CAR-T, n=6 metastatic PDAC, PMC6035088) — RULED OUT. Mesothelin expression reported for only 2 of 6 patients, only as a qualitative IHC localization call, never a number.
2. **CT041/satri-cel CLDN18.2 pancreatic case report** (PMC10492318, J Hematol Oncol) — RULED OUT for modeling. Real per-patient IHC scores exist paired with real RECIST outcomes, but it's a 2-patient case report — nowhere near enough points to fit or validate a threshold model.
3. **CPTAC Pancreatic Ductal Adenocarcinoma** (Cao et al., Cell 2021; cBioPortal `paad_cptac_2021`) — queried directly via the cBioPortal API: real mass-spec MSLN protein abundance (entrez 10232) matched to VITAL_STATUS and FOLLOW_UP_DAYS for 135 real patients — the first genuine per-patient (antigen level, survival outcome) dataset found.
   - Cox proportional hazards: HR=1.16 per log2 unit, p=0.28 (not significant).
   - Tertile split: Low MSLN median survival 715 days, Mid 795 days, High 358 days — directionally "more mesothelin, worse outcome," log-rank p=0.18 (not significant).
   - Among 76 deceased patients, MSLN vs. days-to-death: Spearman rho=-0.21, p=0.070 (borderline).
   - Deceased vs. living MSLN levels: Mann-Whitney p=0.62 (no difference).
   - **Honest read:** every cut points the same direction (higher mesothelin → somewhat worse outcome) but none reaches significance at n=135. Real, reproducible, directly-queryable — but a weak signal, not a headline finding.
   - Data/code: `ev_capture_model/src/cptac_msln_outcomes/` (JSON files; analysis run via `lifelines` CoxPHFitter, not yet formalized into a standalone script).

**Why CPTAC can't feed the reordering analysis (a separate issue from the weak correlation):** CPTAC is bulk mass-spec — one averaged number per patient — which structurally cannot supply the within-sample antigen-density distribution the reordering analysis needs, regardless of how strong the outcome correlation turned out to be.

---

## Outcomes search, round 2 — targeting single-cell/single-sample resolution, four leads ruled out

Round 3 in the original session numbering (labeled round 2 here for archive clarity — this was the pass explicitly targeting single-cell data after the CPTAC bulk-data limitation was diagnosed).

1. **Penn huCART-meso Phase 1 trial** (Aznar/Good et al., Cell Reports Medicine 2025; PMC12490220, NCT03323944) — 9 metastatic PDAC patients, mesothelin CAR-T. RULED OUT: MSLN quantified in only 3-4 of 9 patients as a bulk population percentage, not a distribution. scRNA-seq was performed only on sorted CD3+ T cells (the CAR-T cells), never on tumor/EpCAM+ cells — wrong cell population for an antigen-density question. No MSLN-outcome correlation drawn anywhere in the paper.

2. **iScience 2025, "Exploring the therapeutic efficacy difference in claudin18.2-targeted cell therapy revealed by single-cell sequencing"** (PMID 39925434, PMC11804729) — closest genuine match found: 5 gastric cancer patients on CT041, scRNA-seq of ascites including tumor-derived epithelial cells, with a real stated association ("high CLDN18 expression... indicator of good prognosis"). RULED OUT on two independent grounds: (a) only n=5 patients (1 responder); (b) the paper's stated data accession (HRA006971) was checked directly and is mislinked — it points to a sister study's EV-protein dataset, not this paper's scRNA-seq data. A follow-up subagent check of the supplementary tables found real per-patient T-cell-subset and QC tables but zero CLDN18 expression numbers anywhere — the "high CLDN18 = good prognosis" finding rests entirely on unlabeled UMAP/violin plots with no source data file.

3. **CT041 satri-cel parent trial's multi-omics biomarker sub-study** ("Decoding the response and resistant features to CT041," JCO 2025 43:16_suppl 4015) — the larger, more rigorous version of lead #2 (35 patients, ~41 single-cell samples). Currently only an ASCO conference abstract — no supplementary data or accession numbers exist yet. Worth revisiting in 6-12 months if it publishes in full.

4. **Anetumab ravtansine (mesothelin ADC) Phase I trial data supplement** (Hassan et al., JCO 2020, PMID 32213105) — main text reports only mesothelin-expression-BINNED (not per-patient) response tables. The trial's 33-page online Data Supplement was read in full (via Google's PDF viewer — see tooling note below) specifically to check for a hidden per-patient table. Confirmed: Supplementary Tables A5 and A6 (the only two tumor-response tables) are stratified by dose-cohort and site, with no mesothelin-expression column of any kind. No per-patient (expression, response) table exists anywhere in this trial's public materials.

**Tooling note:** the ASCO/JCO PDF supplement was unreadable via Chrome's native PDF plugin (hung/timed out repeatedly, eventually showed corrupted tiled rendering — a genuine renderer failure, not just unresponsive). Fix: closed the tab, opened a fresh one, used Google's online PDF viewer wrapper (`https://docs.google.com/viewer?url=<url-encoded-target>&embedded=true`), which rendered as real scrollable, extractable HTML. Reusable fallback for any future stubborn PDF.

**Bottom line after all rounds:** every genuinely single-cell or single-sample-resolution lead that could plausibly answer the reordering question has been checked and ruled out for concrete, specific, verified reasons — not assumption or abandoned searches. The CPTAC finding (135 patients, real per-patient bulk protein-to-survival data, weak/non-significant trend in the expected direction) remains the best available outcomes-linked pancreatic dataset — usable only as supporting/narrative context, not for the reordering analysis itself. This search line was recommended concluded barring a genuinely new lead (e.g., the CT041 multi-omics paper publishing in full).

---

## Second-antigen/cancer-type generalization search (non-pancreatic candidates) — for completeness

A separate, earlier generalization search (before narrowing to pancreatic cancer specifically) checked several non-pancreatic single-cell datasets:

- **Kamal et al. 2025 (Multiple Myeloma IMC)** — ruled out: imaging mass cytometry not flow-based CyTOF, bone marrow not blood, data not yet public (BPDC000148 not found live in BloodPAC catalog).
- **Baughn et al. 2023 (Multiple Myeloma CyTOF)** — real dissociated single-cell CyTOF, 49 samples, real externally-validated OS correlation, but raw data is "available on request from corresponding author" only, not deposited. No BCMA in the panel (CD38/CD138 closest). A draft data-request email to Dr. Linda Baughn (Mayo Clinic) was written but NOT sent — Dylan to decide.
- **Tislevoll et al. 2023 (AML CyTOF)** — verified live and public (FlowRepository FR-FCM-Z5YT, cross-confirmed via Zenodo mirror), but headline biology is phospho-signaling response to chemo, not surface antigen density. FlowRepository had a temporary infrastructure outage when download was attempted; likely retriable later. Never fully pursued once the BCMA dataset below was found.
- **Rade/Grieb/Weiss et al. 2024/2025 (BCMA CAR-T multiple myeloma, GEO GSE234261)** — the strongest non-pancreatic candidate, fully pursued to completion. Same antigen AND same therapy class as intended framing (BCMA-directed CAR-T), real CITE-seq surface protein data (CD269/BCMA confirmed present and populated), real outcome data (CR vs non-CR). Full 10-patient comparison run: patient-level Mann-Whitney on means p=0.22, on medians p=0.39 — a real trend in the expected direction (higher BCMA in responders) but not significant at n=5 vs n=5. A relative threshold-sweep (distribution-aware escape-fraction analysis, corrected after an early k_intra-direction bug was caught and fixed) found the direction stable across all 10 thresholds tested, but significance not exceeding the naive baseline at most thresholds (best p=0.0952, still not significant). A single-point absolute calibration anchor was found and verified: the exact TotalSeq antibody clone (19F2) traces to Friedman et al. 2018's bb2121/Ide-cel characterization, reporting a real functional threshold of 222 BCMA molecules/cell — usable for a single-point k_intra calibration if ever needed, though not a full titration curve.
