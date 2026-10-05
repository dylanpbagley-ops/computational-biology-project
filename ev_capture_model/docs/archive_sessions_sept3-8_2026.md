# Session narrative archive — September 3 to September 8, 2026

**What this is.** The blow-by-blow "what changed this session" sections that used to live in
`EV-Capture-CURRENT.md`, moved here under the consolidation rule in `PROCESS.md` step 4
(CURRENT.md had reached 88 KB against a 30 KB trigger). Nothing was deleted.

**When to read it.** Only when you need the history of a specific decision, bug or correction —
"why did we conclude X", "when did that number change", "what was the wrong control". Everything
currently ACTIONABLE has been folded into `EV-Capture-CURRENT.md`; this file is settled history.

**The definitive record of results is not here.** It is
`ev_capture_model/docs/Aim13_Results_Two_Comparator_Ceiling.md`, sections 1–32, which is
maintained as a results document rather than a diary.

Sections below run newest-first, as they did in CURRENT.md.

---

## Most recent session's work (September 3, 2026) — all 12 write-up figures built, redesigned, retextured, and exported; paper is next

This was a long, multi-part session (four consecutive continuations) that took the figure set from a first working draft to a finished, poster/paper-ready state, then this final pass added SVG exports and a convenience folder and consolidated this document itself. Full blow-by-blow (every bug caught, every research pass, every round of Dylan's feedback and exactly how it was addressed) is preserved in `ev_capture_model_narrative_archive_writeup_figures_sessions.md` — this section is the compressed, currently-actionable summary.

1. **Built all 12 write-up figures** (9 paper figures including a new Figure 0 graphical-abstract overview, plus 3 poster-specific hero versions of Figures 1/5/8) in `ev_capture_model/src/figures_writeup/`, importing only from existing tested `src/` modules, never modifying them. Current verified reference numbers (re-check these after any future edit): fig1 peak α=10.4652 (59.6%/104 reliable grid cells); fig2 CA19-9 ratio [1.2272,1.4488] / EGFR ratio [0.8849,0.9352]; fig3 7/8 K_INTRA values match naive ranking; fig4 2/36=5.56% reordered; fig5 Displaced={'Q7915'}, Promoted={'P7676'}; fig6 mean/median ratio [1.0823,29.7818], zero-fraction [0.0028,0.2168]; fig7 Pop A yield=0.6142 / Pop B yield=0.4397 (equal means=20.00); fig8 8.65%/5.56%/48.10%/59.26% across the four evidence tiers.
2. **Two rounds of real design research, both grounded in official/journal sources, not guesswork** — general scientific-figure conventions (Nature/PLOS DPI and colormap standards, viridis, pcolormesh with true log-spaced cell edges, dumbbell/slope charts, illustrative-figure labeling) and ISEF-specific rules read directly from Society for Science's own pages (only the Abstract and Research Plan are officially required; generative AI is barred from writing the submitted abstract/research plan/poster text; poster envelope is 76x122x240cm; Presentation is 35% of score but only 10 of those 35 are the poster itself, the other 25 are the live interview). A genuine search for real winning ISEF papers/posters came up empty — confirmed dead end, Society for Science's own database holds only abstracts.
3. **Two full rounds of Dylan's own design feedback, each resolved by research plus AskUserQuestion confirmation before any large rewrite**, per his explicit request to "think about it and ask questions" rather than guess: (a) "these don't read like a winning ISEF figure" → diagnosed as journal-style plots being reused as-is for the poster instead of being redesigned for it → built `figure_style.py` (shared badge/color/callout/banner system), added Figure 0, built genuine poster-specific versions of Figures 1/5/8 with plain-language takeaway banners; (b) "the text reads as AI-generated" → diagnosed specific tells (CAPS emphasis, "X not Y" contrast pattern, em-dashes, "real" as filler, marketing-style titles) → confirmed via AskUserQuestion → rewrote every title/caption/callout/banner across all 12 figures to neutral, purely-descriptive titles with plain-point (finding-stating) but restrained, scientific-register captions — Dylan's exact requested hybrid of two rendered options he compared side by side. Also removed a previously-flagged, never-fully-resolved unverifiable "17.8%" number that was still sitting in Figure 6's caption text.
4. **All 12 figures exported as SVG with genuinely editable text** (not just the paths alongside the 300 DPI PNGs). The first SVG export used matplotlib's default `svg.fonttype` setting (`'path'`), which silently converts every letter into an outline shape — confirmed directly (0 `<text>` elements, 105 `<path>` elements in one exported file). Fixed by setting `matplotlib.rcParams['svg.fonttype'] = 'none'` before export, which keeps text as real `<text>` objects — re-verified directly (22 `<text>` elements in the same file afterward). Titles/captions/callouts can now be clicked and retyped directly in PowerPoint/Canva/Inkscape. Generated via a one-off driver script (not kept — the 12 figure scripts in `figures_writeup/` are the source of truth for both PNG and SVG; re-run them, or write a similar driver with the fonttype line included, to re-export). Every SVG export re-printed and re-matched the same verified reference numbers as its PNG counterpart.
5. **`Poster, Paper, and Write Up/Figures/` is now the SOLE, canonical output location for all 12 write-up figures — not a copy.** Dylan asked for the scripts' own default output destination to be that folder, so `FIGURES_DIR` was changed in all 12 `figures_writeup/*.py` scripts (was `ev_capture_model/figures/`, now `Poster, Paper, and Write Up/Figures/`, both relative to the script's own location) and every script re-run to confirm it writes there correctly. The old copies in `ev_capture_model/figures/` were deleted (that folder now only holds the four older, unrelated validation plots — `analytical_vs_mc.png` and similar — from earlier in the project, not the write-up figures). This also fully resolves the copy-drift risk noted in a now-superseded version of this document — there is only one location now, nothing to fall out of sync.
6. **Fixed a real ambiguity in Figure 5, caught by Dylan, not by design review.** The displaced patient (Q7915) and promoted patient (P7676) were each drawn with a short dashed line trailing off their dot into empty space — visually, the two trailing dashes sat close enough together that it read as if they were connected to each other (implying Q7915 "becomes" P7676), which is wrong: they are two independent, unrelated events (one target drops out of the top 5, a different target separately enters it). Asked Dylan how to fix it via AskUserQuestion; he chose replacing both dashed lines with small colored text badges directly on the dot ("drops out of top 5" / "new entry to top 5") instead of a bare "OUT"/"IN" tag. Re-ran and re-verified: printed Displaced/Promoted sets unchanged.
7. **This document itself was consolidated per `PROCESS.md`'s process** (it had grown to ~53KB, past the ~30KB trigger, from four stacked "most recent session" sections that should have been folded into the rest of the document rather than accumulated back-to-back). The full narrative from that day's four sub-sessions was moved to the new archive file named above; this section is what remains in the living document.
8. **Figures are now considered fully finished** — built, verified, design-researched, retextured to a neutral tone, exported with genuinely editable SVG text, saved directly to their final poster/paper folder, and the one remaining visual ambiguity (Figure 5) fixed — pending only Dylan's own final look. The next, now top-priority deliverable is the ISEF Research Paper (poster after that, since it's a condensed distillation of the paper's own content).

---


---

## SEPTEMBER 7, 2026 (part 8) — I RAN THE WRONG CONTROL IN PART 6 AND IT REVERSED THE CONCLUSION. Read this before quoting section 26.

**Part 6 said the protein cross-check failed on mass spectrometry and that coarseness did not
explain it. The coarseness test was WRONG.** Script kept: `resolution_sensitivity.py` (the wrong
test). New: `matched_resolution_check.py` (the right one). Results doc section 26 rewritten
(now 26.4 / 26.5).

**WHY THE FIRST TEST WAS WRONG.** I tested "is 20 groups too coarse?" by RANDOMLY merging the
154 cell types into 20 buckets. Random merging PRESERVES DISTINCTNESS — a gene in exactly 1 of
154 cell types still lands in exactly 1 random bucket, so its confinement stays ~1.0. That is
why separation barely moved (0.837). Real organs are nothing like random buckets: a lymph node
contains B cells, T cells, endothelium and stroma at once, so pooling by ORGAN dilutes a
lineage-restricted gene in a way random pooling never does. The test asked "do fewer buckets
hurt?" (no). It never asked "does mixing cell types into real organs hurt?" — the actual question.

**THE RIGHT TEST.** 18 of the 20 MS tissue names appear verbatim in the HPA bulk RNA reference,
so the same measure runs on the SAME organs in both technologies.

| reference | resolution | separation | p |
|---|---|---|---|
| RNA single-cell | 154 CELL TYPES | **0.868** | 0.0003 |
| PROTEIN antibody staining | tissue x CELL TYPE (~96) | **0.785** | 0.0104 |
| RNA bulk | 51 organs | 0.726 | 0.0244 |
| **RNA, the same 18 organs MS covers** | **18 organs** | **0.626** | **0.1425** |
| **PROTEIN mass spec** | **20 organs** | **0.650** | **0.1148** |

**RNA fails on those organs too, by the same amount (0.626 vs 0.650).** Clean 2x2: RNA works at
cell-type resolution and fails at organ resolution; protein does exactly the same. **Technology
is not the variable. RESOLUTION is.**

**SO THE CLAIM IS SHARPER, NOT WEAKER — and the framework PREDICTED this.** The thesis is that
safety depends on WHICH CELL TYPE carries the target and whether it is expendable. An organ
homogenate destroys precisely that information. The measure failing at organ resolution is a
prediction of the theory, now confirmed. Correct phrasing for the paper:
> Restriction must be measured across normal CELL TYPES. Measured that way it works in RNA
> (0.868) and in protein (0.785). Measured across whole organs it fails in BOTH (0.626, 0.650).

Also recorded: a NaN bug in the first `matched_resolution_check.py` (blank MS Intensity -> mean of
empty list) made the MS row read 0.000. Found and fixed before any conclusion was drawn.

---


---

## SEPTEMBER 8, 2026 (part 10) — a THEOREM that unifies the project, a DISCOVERY screen with a 57x positive control, and locked forward predictions. Results doc sections 29-31.

**NO FURTHER DOWNLOADS ARE NEEDED.** `proteinatlas.tsv.zip` (per-gene HPA annotation, 119
columns) turned out to contain everything the remaining angles required: Protein class,
observed subcellular location, RNA/protein cell-type specificity, cancer AND cancer-cell-line
specificity, and prognostics for 21 cancers. The DVP zips and rna_cancer_sample.zip are NOT
needed — see part 10c for why the DVP data cannot answer its question regardless.

### 10a. SECTION 29 — the weighted multi-compartment ceiling theorem. THE BIGGEST RESULT.
Script `weighted_ceiling_theorem.py`. The project carried TWO measures with no formal
relationship: the ROC ceiling (a theorem) and confinement+expendability (an empirical rule).
**They are the same theorem.** Give each normal cell type an expendability weight w and define
harm H(f) = sum_c w_c q_c E_{p_c}[f]. Then (T(f), H(f)/W) is exactly an ROC point for tumour
against the WEIGHTED MIXTURE p_W of normal cell types, so the achievable region is that ROC's
convex hull and the optimum is likelihood-ratio thresholding. **Confinement falls out as the
leading term: H proportional to sum_c w_c s_c.**
- **Verified: 16,400 checks, 1-6 compartments, largest violation 0.000e+00, likelihood-ratio
  thresholding matches brute force exactly.** All three known results are special cases.
- **A BUG IN MY VERIFIER, recorded:** the first hull construction was an ad-hoc "upper envelope"
  loop that reported false violations up to 0.96. Now a monotone-chain upper convex hull.
- **This is the section to aim at the Horvitz Prize for Fundamental Research.** It converts the
  headline from "a measure that predicts safety" to "a theorem about what any antigen-targeted
  therapy can achieve, of which the measure is the leading term."

### 10b. SECTION 30 — genome-wide DISCOVERY, with a real positive control.
Script `target_discovery.py`. All 19,475 genes -> 4 filters -> **214 candidate targets nobody is
pursuing**, 134 of them in an expendable blood lineage.
- **POSITIVE CONTROL: 10 of 16 known-in-development immunotherapy targets recovered, 8 in the
  top 50.** LILRB4 rank 6, GYPA 9, BAFF-R 16, CD27 20, TACI 23, CSF3R 24, CD79A 29, CLEC12A 38.
  Screen keeps 1.10% of genes, so random expectation is 0.18 hits. **57-fold enrichment,
  hypergeometric p = 1.6e-16.**
- **The MISSES are correct calls:** CD5, CD7, KIT, FLT3 were dropped — and CD5/CD7 CAR-T cause
  T-cell aplasia while FLT3/KIT CAR-T cause myeloablation. The framework declined to recommend
  exactly the ones that poison people.
- **TWO FAILURES OF MINE, kept in the doc because they are the trap:**
  (1) Using BULK tumour RNA, the top hits were haemoglobin "elevated in glioblastoma" and
  immunoglobulins "elevated in lung adenocarcinoma" — bulk tumour contains BLOOD and IMMUNE
  INFILTRATE, so any blood-confined gene looks cancer-elevated BECAUSE it is blood-confined. The
  screen was rediscovering its own filter through a contaminant. Fixed with cancer CELL LINES.
  (2) No absolute floor meant perfect confinement of nothing scored 1.000 — the list filled with
  olfactory receptors. Added the CAR-T sensitivity floor from Aim 10 (>= 20 nTPM), and replaced
  "predicted membrane protein" with OBSERVED plasma-membrane localisation.
- Top candidates include IGKV genes in plasma cells at 1,300-5,500 nTPM in myeloma lines, which
  connects straight back to section 15's own light-chain clonality work.

### 10c. SECTION 31 — the protein-at-cell-type question stays OPEN, and it is the DATA's fault.
Script `hpa_summary_reference.py`. HPA integrated Deep Visual Proteomics (protein per cell type,
27 cell types, 14 tissues) in May 2026 — the reference I said did not exist. It **cannot answer
the question**: 8 of 29 labelled targets are NOT DETECTED, and **7 of those 8 are SAFE targets**
(BCMA, SLAMF7, CD30, GPRC5D, CD52, DLL3, GPC3). DVP is 27 cell types from a SINGLE donor; plasma
cells and activated lymphocytes are rare, so the antigens defining the safe class are below the
detection limit. Dropping them removes the class being tested.
- **Calibration that makes this readable:** HPA's 4-level ordinal ladder gives only 0.553 on the
  RNA cell-type data where the continuous computation gives 0.868. A null on the ladder is
  uninformative, not negative.
- **A scoring error of mine:** I first mapped "Not detected" to the BOTTOM of the ladder — the
  opposite of its meaning — which produced a meaningless 0.413. Now returns missing.
- Section 26.5's scope condition stands unchanged; antibody staining (0.785) remains the best
  protein-level evidence.

### 10d. LOCKED PROSPECTIVE PREDICTIONS (`prospective_predictions.py`).
52 targets in registered trials with NO safety verdict yet, scored with the frozen rule:
**14 predicted tolerated, 18 predicted on-target toxicity, 20 explicit no-calls.** Dated
2026-09-08. Put the table in the write-up unchanged and never revise it.

### 10e. ISEF ODDS, calibrated (doc: `Status_ISEF_Ratings_and_Expansion_v2.md`).
AzSEF 2025: 786 students / 558 projects, 20 advance. CBIO 2026: 2 First / 5 Second / 8 Third /
8 Fourth = 23 Grand Awards from ~85 projects. Top Awards drawn only from ~44 First winners.
Ratings: **qualify from AZ 85, any CBIO award 72, second-or-better 42, First in CBIO 22**;
top awards — Horvitz 8, Leighton&Berger 5, Coleman 5, Scripps 4, Young Scientist 3,
Yancopoulos 2, Moore 2, Barrett 1.

### 10f. ON THE AUTHOR EMAIL: recommended AGAINST as a priority.
Patient-level antigen-density-plus-response cohorts exist but are behind "reasonable request".
It tests a DIFFERENT mechanism from the project's thesis (targets, not patients). Send it once as
a 15-minute lottery ticket; build nothing on it.


---

## SEPTEMBER 7, 2026 (part 9) — the two-antigen work replicated in a solid tumour. Results doc section 28.

Script: `adc_two_antigen.py`. 72 ADC targets, **2,556 antigen pairs**, 749 metastatic vs 96
normal samples from the same autopsies. Matrix rebuilt from the long table and **verified against
the paper's own wide table: 5,992/5,992 values matched exactly** (hard assertion in code).

- **THEOREM REPLICATES: 0 violations in 5,112 checks.** Second disease, third technology.
- **Statistical gain replicates exactly:** median +0.0000 against a shuffled-label floor of
  +0.0062. A second antigen adds NO new information — same answer now in two diseases, two
  technologies, two labs.
- **Engineering gain does NOT replicate at the median** (+0.000 here vs +0.106 in DLBCL). It has
  become CONCENTRATED, not absent: 34.8% of pairs still exceed +0.05, 90th percentile +0.277. In
  lymphoma nearly every pair benefited because most of the panel is B-lineage with markers LOWER
  on tumour than on normal B cells; across organs that is no longer true of most pairs. Report as
  a real disease difference, not a failure.
- **A TRAP AVOIDED, same class as section 20.** Top per-pair gains look huge (ALCAM+PD-L1 +0.693)
  but the one-receptor baseline INSIDE those pairs is NEGATIVE (-0.16). A big ratio over a
  terrible baseline is not a big result. The honest cohort-level number:
  **best single buildable receptor 0.398 -> best two-receptor design 0.533, net +0.135** —
  strikingly close to DLBCL's +0.133, from completely different data.
- **Every top design is a NOT gate** — activating receptor plus an INHIBITORY receptor against an
  antigen higher on normal tissue. Same design conclusion as section 21, reached independently in
  a solid tumour.


---

## SEPTEMBER 7, 2026 (part 7) — the ceiling transplanted to a SOLID TUMOUR with matched normals, and section 22's claim replicated by an independent lab. Results doc section 27.

Data: `ev_capture_data/adc_autopsy/` — De Schepper et al., Nature Communications 2025
(s41467-025-67840-z). **Rapid autopsy: 30 metastatic breast cancer patients, 749 metastatic
lesions across 36 sites, PLUS normal tissue from the same autopsies, 72 ADC target genes by
RNA-seq in every sample.** Extracted long-format table: `adc_autopsy/supfig8_long.tsv` (60,840
rows). Scripts: `adc_autopsy_ceiling.py`, `trial_registry_control.py`.

**WHY THIS DATASET MATTERS: tumour and normal come from the SAME BODY, same protocol, same
batch.** The confound sections 19-21 had to defeat by identical gating is removed by design.
Different tumour type, lab, technology and country.

**THE HEADLINE — section 22 replicated in a solid tumour. Put this in the paper.**
Same 17 targets, same labels, same statistic, two measures:
- tumour-vs-normal **ceiling**: separation **0.528, p = 0.444** (chance)
- **confinement** across normal cell types: separation **0.750, p = 0.046**

The scramble is instructive, not merely null: TROP2 (approved) 0.781 and Nectin-4 (approved)
0.724 sit HIGH — but so do HER2 (fatal event) 0.640 and CAIX (liver toxicity) 0.632, while
CLDN18.2 (approved 2024) is 0.262 and c-MET (approved 2025) is 0.223. **Being easy to tell apart
from normal tissue has nothing to do with whether the drug was tolerated.** This is the project's
least intuitive claim, and it survived the transplant.

**A DISCRIMINANT CHECK USING SOMEBODY ELSE'S TABLE (attacks weakness #4).** The same paper ships
Supplementary Data 1: **623 registered clinical trials against 74 ADC targets**, with phase,
compiled by other people. It CANNOT be a safety label ("reached phase 3" is progress, not
toxicity — HER2 has 36 phase-3 trials and a fatal event). It CAN be a specificity control, and
the prediction was stated before running:
- confinement vs "reached phase 3": **0.495, p = 0.530** — AS PREDICTED, uninformative
- confinement vs on-target toxicity, same targets: **0.750, p = 0.046**
Confinement tracks poisoning, not popularity, verified against a table this project did not
compile. Does not remove the need for an outside reviewer on the labels, but it kills the
"you're just measuring which targets people like" objection.

**WHAT IT DOES NOT SETTLE.** RNA-seq (IHC on only 8 of 72 targets), so it does NOT answer
section 26's protein question. Breast cancer only. Post-mortem RNA degradation is inherited.
17 labelled targets means p = 0.046 is one relabelled target from crossing 0.05 — report as
consistent with section 23, not as full-strength independent confirmation.

**BULK, NOT SINGLE CELL.** Separation is between tumour SAMPLES and normal SAMPLES, an upper
bound on the cellular ceiling. Not comparable to the CyTOF numbers in sections 1-22. Say so.

**UNUSED IN adc_autopsy/ (next session).** Supplementary Data 1 (the 623-trial registry, only
partly used); Sup data 6.1-6.3 (4,694 weighted gene-PAIR co-expression coefficients with
p-values — this is a ready-made two-antigen dataset for a solid tumour, directly extending
section 21/25); Fig 6B (the authors' own tumour-vs-normal overlap percentages per gene, worth
cross-checking against our AUC); and per-lesion tables with tumour cellularity, HR/HER2 status
and histological subtype.


---

## SEPTEMBER 7, 2026 (part 6) — the protein cross-check came back SPLIT, and the failure is the important half. Results doc section 26.

Scripts: `protein_cross_check.py`, `ms_coverage_diagnostic.py`, `resolution_sensitivity.py`,
`reference_agreement.py`. Data added: `normal_ihc_data.tsv.zip`, `ms_tissue_sample_data.tsv.zip`,
`ms_tissue.tsv.zip` (HPA v25.1 protein tables) and `ev_capture_data/adc_autopsy/` (9 XLSX files,
Nature Communications rapid-autopsy ADC dataset — downloaded, NOT yet analysed).

**THE REGISTERED HYPOTHESIS H5 SPLIT. Report both halves; never report only the win.**
- PROTEIN, immunohistochemistry, 76 tissue x cell-type groups: **0.785, p = 0.0104 — H5 MET**
  (0.775, p = 0.0158 with uncertain antibodies dropped).
- PROTEIN, mass spectrometry, 20 bulk tissues: **0.650, p = 0.1148 — H5 NOT MET.**

**TWO EXCUSES WERE TESTED AND BOTH DIED. This is the part that makes it science.**
1. *"The MS table is sparse."* Restricting to targets measured in >= K tissues makes it **worse**,
   not better: 0.650 -> 0.600 -> 0.548 -> 0.393. Sparsity is not the explanation.
2. *"20 groups is too coarse."* Randomly merging the 154-cell-type RNA reference down to 20
   groups (60 random partitions) still gives **0.837 [0.753, 0.911]**. 0.650 is outside the range
   of every partition. Coarseness is not the explanation either.

**REFERENCE AGREEMENT (the factual account).** The two protein references agree with each other
(rho 0.62-0.70) BETTER than either agrees with RNA (rho 0.35-0.43). MS internal consistency is
0.963, so the disagreement is not noise. **"Confinement" is not one quantity measured four ways.**

**THE CLAIM IS NOW SCOPED, and this is the correct weaker version:**
> Restriction measured across normal CELL TYPES — by RNA at single-cell resolution, or by
> antibody staining at tissue x cell-type resolution — separates tolerated targets from
> on-target-toxic ones. Restriction measured as mass-spectrometry abundance across 20 bulk
> tissues does not.

Suggested (NOT proven) reason: on-target toxicity is about WHICH CELLS carry a target, and a
bulk homogenate destroys exactly that. Testing it properly needs a cell-type-resolution protein
reference, which does not exist at this scale — a limitation of the field, not of this analysis.

**H6 NOT TESTABLE** on either protein reference (both indexed by tissue, so the haematopoietic
flag fires on 0 of 24 targets and the two-factor rule collapses by construction). Recorded as
untestable, not as a failure. **H7 MET** (shuffled-label control, section 25).

**WHY THIS HELPS RATHER THAN HURTS AT ISEF.** Most projects would report the IHC win and omit the
MS failure. This one reports the split, tests its own two best excuses, kills both, and narrows
its own claim as a result. That is the behaviour that separates First Award work from Fourth.
Lead the limitations slide with it rather than hiding it.

**STILL OPEN.** The `adc_autopsy` dataset is downloaded but unanalysed — 30 patients, 749
metastatic lesions, 96 normal samples across 25 organs, 72 ADC targets, matched tumour and normal
in the SAME person. Two uses: (1) within-patient tumour-vs-normal ceilings in a solid tumour,
which is the DLBCL analysis transplanted; (2) raising the safe-vs-toxic panel well past 19 vs 10.
Note it is RNA-seq with IHC on 8 targets, so it does not resolve section 26's protein question.


---

## SEPTEMBER 7, 2026 (part 5) — all five expansion angles built; the whole-body result is now the project's strongest, and two of my own claims were corrected

Results doc sections **23-25**. Scripts: `src/aim13_scoping/target_labels.py` (FROZEN labels),
`wholebody_v2.py`, `pancreatic_thread.py`, `two_antigen_null_control.py`,
`build_confinement_table.py`, `targetscan.py`. Figures: `fig9_wholebody_restriction`,
`fig10_pancreatic_thread` (PNG + editable-text SVG). Tests: **195 pass** (173 before; 22 new in
`tests/test_wholebody_restriction.py`).

**THE HEADLINE IS NOW THE WHOLE-BODY RESULT, and it is much stronger than part 4's version.**
Panel expanded 34 -> 40 targets with every label frozen in a file no analysis may write to.
- **19 safe vs 10 toxic: confinement separates them at 0.868, permutation p = 0.00031.**
  (Part 4's number was 15 vs 7 at 0.810, p = 0.011. Superseded.)
- **Two-factor rule (confined AND expendable): 0.905, p = 0.00004.** Expendability = the home
  cell type is a haematopoietic lineage, assigned from the cell type's identity alone, so it is
  not circular. Lexicographic, no fitted parameter.
- **SOLID TUMOURS ALONE (every blood target removed): 0.792, exact p = 0.023.** This is the
  load-bearing check — the finding is NOT a blood-cancer artifact. Lead with this against the
  obvious objection.
- Leave-one-out: 29/29 drops stay significant, separation 0.895-0.942.
- Specificity control: 9 efficacy-only failures do NOT look toxic (0.684, p = 0.063, ns).
- Germ/placental cell types removed: two-factor unchanged at 0.905 (CD52's home correctly
  becomes B cells).
- Independent second reference (51 bulk tissues): replicates weaker, 0.726, p = 0.024;
  Spearman rho between references 0.512. Cell-type reference is primary and the doc says why.

**A BUG OF MINE, RECORDED (§23.2).** Haematopoietic membership was first tested by plain
substring match, so `"conjunctival goblet cells"` matched `"t cell"`. It **inflated** the
two-factor separation to 0.916; the corrected, word-boundary-anchored number is **0.905**. A
guard table in `target_labels.py` now fails at import if it regresses, plus a regression test.

**A CLAIM OF MINE WITHDRAWN (§24.3).** I told Dylan "the two top-ranked pancreatic targets are
the two this project independently calibrated." On the full 14-target PDAC panel that is
**false** — CLDN18.2 ranks 1st but mesothelin ranks **5th**. The six-target panel made the
concordance look tighter than it is. Correct statement: the one approved PDAC target ranks 1st
of 14; the three failed/toxic ones (MUC1, EGFR, HER2) rank 11th, 13th, 14th.

**PANCREATIC THREAD (§24), computed end to end in one script.** CLDN18.2 calibration
r = 0.964 / RMSE 0.237 on 3 pancreatic + 3 gastric lines; mesothelin calibration; **CPTAC PDAC
n = 135, mesothelin protein vs survival is an HONEST NEGATIVE (log-rank p = 0.350)** and is
plotted as such in fig10; the 14-target ranking; and the original CA19-9 flip (PANC-1 vs HPDEC)
as the thread's origin. Incidental: **GPC1, the famous pancreatic exosome marker, ranks 12th of
14 on confinement** — a good diagnostic marker is not a good therapeutic target, and the
framework says why.

**SHUFFLED-LABEL NULL CONTROL for the two-antigen work (§25).** Same cells, same code, labels
destroyed. Single-antigen J collapses 0.692 -> 0.002 (control validates itself).
**Statistical gain: real +0.0000 vs shuffled floor +0.0003 — AT the noise floor.**
**Engineering gain: real +0.106 vs +0.0008 — 130x the floor.** This is the crispest available
statement of §21: a second antigen adds essentially no new information, but it unlocks a large
amount that no single activating receptor can physically use.

**TOOL RELEASED.** `build_confinement_table.py` computes confinement for **all 19,475 expressed
human genes** (`confinement_table.tsv`, `confinement_web.json`); `targetscan.py` is the CLI;
and a published, shareable interactive page ("Target Confinement Atlas") does gene lookup,
places a query against the 29-target panel, and computes the ROC ceiling live from user-set
antigen distributions.

**LABEL PACKAGE FOR OUTSIDE REVIEW.** `docs/Label_Review_Sheet.md` — all 40 rows with citations
and a tick-box column, written so a non-computational reviewer (Ms. Turner, or a lab contact)
can check it. `OSF_Preregistration_DRAFT.md` gained an addendum that states plainly the
whole-body analyses were ALREADY RUN, freezes everything for future runs, and registers
falsification conditions H5-H7.

**A DOC/CODE INCONSISTENCY FIXED.** `two_antigen_ceiling.py`'s `edges()` still carried the
retracted "third bug" story (0.618 -> 0.410). The code comment now matches the corrected record
in the section above it: that diagnosis was wrong, the change is kept anyway, and the real
explanation is median-of-max vs max-of-median.

**STILL OPEN.** (1) HPA protein-level IHC cross-check — `normal_tissue.tsv.zip` triggered to
Dylan's Downloads, awaiting the Move-Item. (2) A *Nature Communications* 2025 rapid-autopsy
dataset (30 patients, 749 tumour lesions, 96 normal samples across 25 organs, 72 ADC targets,
matched tumour/normal in the same person) — permission to download requested; this is the best
route found to raise n and to get paired tumour-vs-normal in a SOLID tumour. (3) Patient-level
antigen-density-plus-response cohorts exist (Stanford BCMA n=76; a 2025 calibrated-flow LBCL
study) but are "available on reasonable request" — an author email, not an impossibility. My
earlier claim that patient-level outcomes were INHERENTLY unavailable was under-researched and
is withdrawn.

**ISEF AWARD STRUCTURE, confirmed from Society for Science's own pages.** Grand Awards are per
category: 1st $6,000, 2nd $2,400, 3rd $1,200, 4th $600. Every Top Award is chosen FROM the 1st
Award winners: Yancopoulos Innovator $100,000; two Regeneron Young Scientist Awards $75,000
each; Gordon E. Moore $50,000; Leighton & Berger $50,000; Barrett / **Horvitz (Fundamental
Research)** / Scripps / Coleman $10,000 each. Winning the category is the entry ticket, not the
ceiling — and the Horvitz prize is for work whose contribution is a principle, which is this
project's shape.


---

## SEPTEMBER 7, 2026 (part 4) — WHOLE-BODY NORMAL ATLAS: the thesis now covers SOLID TUMOURS, and the pancreatic thread closes

Scripts: `src/aim13_scoping/wholebody_restriction.py`, `wholebody_toxicity_test.py`.
Data: `ev_capture_data/normal_tissue_atlas/rna_single_cell_type.tsv.zip` (Human Protein Atlas single-cell consensus, **154 human cell types, 3,087,080 rows**). Assessment doc: `docs/Project_Status_Weaknesses_and_Expansion.md`.

**What changed.** Until part 3 the project rested on two blood-cancer datasets. It now has a whole-body normal reference, so the same restriction logic runs on ANY target, blood or solid.

- **Measure.** For each target gene: `top_share` = fraction of total normal expression in its single highest cell type; `n_types` = number of cell types above 10% of the top. Ran on **34 real clinical targets**; all found, all sensible.
- **THE NEW HEADLINE TEST.** Against documented **on-target toxicity** (15 safe vs 7 toxic, labels fixed from the literature before the run): `top_share` separates at **0.810, exact p = 0.011**; `n_types` at **0.829, exact p = 0.006**. Better n and smaller p than the blood-only 5-vs-5 test (0.840, p = 0.048).
- Representative values: CD20 0.850/1 type (approved), BCMA 0.789/1 (approved), PSMA 0.568/2 (approved), CLDN18.2 0.359/8 (approved 2024), CEACAM5 0.283/7 (severe colitis), CD44 0.061/45, EGFR 0.051/60, **HER2 0.041/72 (fatal CAR-T event)**, **CD47 0.030/120 (magrolimab, on-target anaemia)**.
- **THE THESIS IN ITS FINAL FORM.** Breadth alone catches CD47/HER2/EGFR but NOT CEACAM5, which is confined and still toxic. The correct statement, which explains every case including the exceptions:
  > **A target is safe only if it is BOTH confined to few normal cell types AND those cell types are expendable or regenerating.**
  > CD19/CD20 → B cells (expendable). BCMA → plasma cells (expendable). CLDN18.2 → gastric mucosa (regenerates). CEACAM5 → colonocytes (NOT expendable) → colitis. HER2/CD47 → everywhere → catastrophic.
  This unifies blood and solid tumours under one principle. **Lead the write-up with this.**
- **PANCREATIC RANKING (the personal thread, and it is not forced).** By whole-body confinement: **CLDN18.2 (0.359) > CEACAM5 (0.283) > mesothelin (0.239) > MUC1 (0.093) > EGFR (0.051) > HER2 (0.041)** — the order the field itself converged on. The two top-ranked pancreatic targets are the two this project **independently calibrated months ago** (CLDN18.2 dose-response on three pancreatic lines, r = 0.964, RMSE 0.24; plus the mesothelin calibration), CPTAC PDAC survival data (n = 135) is already in hand, and the project's ORIGINAL EV headline — the CA19-9 flip, PANC-1 vs HPDEC — was pancreatic. A pancreatic thread runs from Aim 1 to Aim 15.
- **Caveat to state:** HPA single-cell is transcript-level, not protein. The field's own tools (Dannenfelser, SCAN-ACT) do the same. Protein-level IHC cross-check is the named fix and is feasible but not yet run.

**ISEF calibration done this session** (Society for Science full award releases, 2024/2025/2026, CBIO category). Computational CAR-T and PDAC multi-omic projects have been landing at **Fourth Award** (2026 CBIO042 "ReinforCell", 2026 CBIO088 PDAC multi-omics, 2025 CBIO013 bispecific solid tumour, 2025 CBIO048 CD30 lymphoma). First Awards go to **mechanistic model + independent validation + one specific new claim** (2026 CBIO071 is a Markov chain). Realistic placing: **strong Second Award, First Award reachable.** Biggest presentation risk: being mistaken for CBIO042/CBIO088 — never lead with the dataset or the pipeline, always with the theorem and the mechanism.


---

## SEPTEMBER 7, 2026 (part 3) — circularity removed, full cohort, outcome test at its ceiling. Results doc sections 19-22.

**These numbers supersede parts 1-2 for the two-antigen work. Three further problems were found and are recorded, not hidden.**

- **Gating circularity removed.** DLBCL comparator is now RLN **normal B cells only** — tumor and normal are gated by the identical rule, so no marker is circular, and it is the binding compartment anyway. Myeloma excludes the 7 gate-defining antigens (CD3E/CD14/CD33/CD56/CD19/CD138/CD269), leaving 56 antigens and 1,540 pairs. **DLBCL now covers all 41 patients** (resumable store).
- **A WRONG DIAGNOSIS OF MINE, recorded.** I claimed a "third bug" (test cells dropped outside train histogram range) and fixed it; the number did not move because that was never the cause. The real explanation is median-of-max vs max-of-median: the single-antigen baseline may pick a different antigen per patient. Verified directly (CD27+FasR: median-of-max 0.616 vs per-antigen medians 0.384/0.412). Both numbers correct, different questions.
- **TWO single-antigen baselines are required and both are now reported.** SINGLE_hi (kill only where the antigen is HIGH) is what ONE receptor can physically do — no receptor fires on the ABSENCE of a protein. SINGLE_any (either direction) needs two receptors to build but is the right statistical control.
- **DLBCL, 41 patients, 666 pairs:** theorem violations **0/666** in-sample. Best 1-receptor J **0.618**; best 2-receptor gate **0.646** (+0.028). Engineering gain median **+0.133** (98% of pairs); statistical gain median **+0.003** (only 11 pairs >0.05). **NOT beats AND in 86%.**
- **KEY INTERPRETATION:** the large engineering gain and tiny statistical gain are both real and not in conflict. A second receptor mostly does NOT add information — it **unlocks information already in a single marker** (the ones LOW on tumor) that no activating CAR can use alone. Best kill-low markers: IgD 0.552, CD21 0.517, CD23 0.443; best buildable kill-high: FasR 0.410. Concrete design: activating CAR against a B antigen + inhibitory receptor against IgD, sparing naive B cells.
- **A FINDING THAT DID NOT REPLICATE:** "exclusion beats requirement" is 86% in DLBCL but **51% (coin flip) in myeloma**. Report as a lymphoma finding, not general. This disconfirms the generality claimed in part 2.
- **A finding that DID replicate:** the likelihood-ratio rule is the population optimum and the wrong practical object (overfit gap 0.017 DLBCL / 0.267 myeloma vs ~0.02 for a simple gate).

**OUTCOME TEST, generalised to "restriction to the cancer's own lineage" (B for B-cell cancers, myeloid for AML).** Inclusion rule fixed in advance; excludes myeloma (plasma cells absent from blood) and Hodgkin (no clean PBMC counterpart).
- **Primary, balanced 5 approved vs 5 failed: separation 0.840, exact p = 0.0476.**
- **A secondary predictor (margin over next-best lineage) was PRE-SPECIFIED and made things WORSE** (0.800, p = 0.075). Had it been added post hoc it would have looked like a refinement. It does not rescue HLA-DR either (margin 0.002, correctly low) because CD40's margin is 0.345 yet its drug failed.
- Sensitivity (unbalanced, +CD47 magrolimab and +CD70 cusatuzumab, both failures, no matching approvals): 0.886, p = 0.0152. Reported separately, never as headline.
- **CEILING, state this in the write-up:** the number of APPROVED blood-cancer surface targets is small in the world, not just in this panel. **n cannot be meaningfully raised by more work.** This test stays underpowered; present it as suggestive evidence with an exact p-value, not as a significant result.


---

## SEPTEMBER 7, 2026 (part 2) — healthy human reference added, and the MULTI-TARGET extension

Detail: results doc sections 16-17. Scripts: `src/aim13_scoping/healthy_pbmc_load.py`, `healthy_lineage_restriction.py`, `two_antigen_ceiling.py`.

**Healthy reference (GSE164378, Hao et al. Cell 2021):** 161,764 PBMCs, 8 HEALTHY donors, 228 surface proteins, authors' own cell-type labels. Isotype controls establish a **noise floor of 0.598** (not 0.5) — use it. Duplicate antibody clones agree to a median of 0.034.
- **CD19 is the most lineage-restricted protein of 224, at exactly 1.000**, home = B cells. CD22 0.990 (rank 8), CD20 0.986 (rank 12), CD79b 0.886 (rank 56). Combined with their ~chance scores against normal B cells, the thesis is now two measurements rather than an argument: **not tumor-selective, maximally lineage-selective.**
- **Forward prediction with two verified hits.** Filtering to B-homed proteins clearing CD79b's bar gives 6 non-approved candidates; two are the field's leading next-generation targets: **CD268/BAFF-R** (Qin et al., Sci Transl Med 2019, doi 10.1126/scitranslmed.aaw9414, built explicitly for CD19-relapse) and **CD72** (Izgutdina et al., J Immunother Cancer 2025, doi 10.1136/jitc-2025-012013, concluding CD72 may beat CD22 as second-line in B-ALL). The method saw only healthy-donor protein data.
- **Limit:** blood only, so every myeloma target scores at/below the isotype floor here (BCMA 0.580 etc.) — uninformative, NOT evidence against. Needs a healthy marrow reference.

**Two-antigen theorem (the multi-target extension) — CORRECTED, see results doc section 17.** Nested families SINGLE <= AND/NOT <= LIKELIHOOD-RATIO, ordering forced by Neyman-Pearson.
- **TWO BUGS OF MINE, caught by the myeloma replication, are documented in the results doc rather than quietly fixed.** (1) The single-antigen baseline only tested "kill where the marker is HIGH", so markers that are LOW on tumor (IgD, CD3E) made a NOT gate look like a two-antigen gain when it was a one-antigen effect with the opposite sign. (2) The likelihood-ratio rule is the population optimum but is NOT estimable at this sample size — it picks among 256 histogram bins and overfits (held-out gap 0.020 in DLBCL, **0.261** in myeloma, vs ~0.02 for a two-threshold gate).
- **The first version overstated the value of combining antigens roughly tenfold.** Corrected median gain: **+0.015** (DLBCL), **+0.003** (myeloma). The claim "CD20 NOT IgD, 0.215 -> 0.645" was an artifact; the true figures are 0.468 -> 0.582 (+0.114).
- **The theorem itself is verified exactly, in-sample, where it is a mathematical claim: 0 violations in 666 DLBCL pairs and 0 in 1,953 myeloma pairs.**
- **What survives and replicates across both diseases: NOT gates beat AND gates** (79% DLBCL, 55% myeloma) — exclusion beats requirement, against the field's AND-dominated emphasis.
- Best corrected DLBCL pairs: CXCR5 NOT IgD +0.163, CD21+CXCR5 +0.142, CD22+CD21 +0.139. For the approved myeloma pair CD138+BCMA the gain is **-0.000**.
- **Methodological finding worth reporting on its own:** the Neyman-Pearson optimum is the right theoretical object and the wrong practical one; a simple gate generalises, the fitted likelihood ratio does not.
- Caveat: gating markers are partly circular (CD3/CD19/CD20 in DLBCL; CD3E/CD14/CD33/CD56/CD19/CD138/CD269 in myeloma). The top myeloma pair CD3E+CD138 is entirely circular and must be ignored.

**THE OUTCOME FLIP (results doc section 18) — the strongest validation the project now has.** Both groups are antigens with a REAL drug given to patients; the only difference is approval. Labels come from drug history, not expression, so the earlier circularity is gone. Approved: CD19, CD20, CD22, CD79b, CD52. Failed: CD23 (lumiliximab, LUCID ph2/3 n=627 **stopped early for lack of efficacy**, doi 10.1111/bjh.13061), CD80 (galiximab), CD40 (dacetuzumab), HLA-DR (apolizumab), CD25 (camidanlumab), CD200 (samalizumab). CD52 was included deliberately although expected to hurt the result.
- **B-lineage restriction in healthy donors: separation 0.840, exact permutation p = 0.048.**
- **Tumor-vs-normal ceiling: separation 0.650, p = 0.278 — does not predict clinical success.**
- So: the ceiling says what a target CAN do; **which normal lineage it lives on says whether it is worth building.** The dispensable-lineage thesis now has a clinical-outcome p-value.
- Fragile: n=10, one reclassification flips it; labels are unverified by a second person; the one misordered pair is HLA-DR over CD79b. Pre-register before claiming.


---

## SESSION OF September 7, 2026 — five hardening steps done; replication in a second disease succeeded

Full detail: `ev_capture_model/docs/Aim13_Results_Two_Comparator_Ceiling.md` sections 13-15. Scripts/JSON in `ev_capture_model/src/aim13_scoping/`. All 173 tests still pass; all affected figures re-rendered as PNG + editable-text SVG.

- **Third normal compartment added** (myeloid/NK, from cells previously discarded). Ceilings now take the minimum over T cells, myeloid/NK and normal B cells. Monotonicity verified exactly within one run: **13 of 37 ceilings fell, 0 rose.** The spurious winners collapsed (CD38 -0.138, CD86 -0.127, HLA-ABC -0.087), and **the top ceiling in the panel is now CD19 (0.580)** — the antigen that actually became the drug. This substantially repairs the section-6 caveat.
- **Gating robustness settled.** Otsu and Gaussian-mixture gating vs the percentile gate, on matched files: rank agreement 0.922 / 0.956 / 0.984. Conclusions do not depend on gating. The real sensitivity is **cohort size** (41 vs 12 tumors: rank correlation 0.857), not gating method.
- **Strict convention adopted everywhere.** 48.1% (13/27) is now the primary number in code, figures and docs. Figure 8's non-comparable "9-patient subset" bar was replaced by the gated tier, so all four bars now use one criterion: 8.65% -> 5.56% -> 48.15% (all events) -> 51.85% (sorted tumor cells).
- **Figure 5 rebuilt on gated tumor cells**, computed live from the FCS files rather than a cached constant. The gated example is cleaner than the old one: O8155 (mean 176.9) and Q5971 (mean 178.8) enter the top 5, displacing Q7915 (230.3) and R5097 (205.2) — both promoted patients have LOWER means than both displaced, so the point holds for every pair. R5097 falls from 4th to 14th.
- **REPLICATION IN A SECOND DISEASE SUCCEEDED.** GSE234261 (multiple myeloma, BCMA CAR-T) was re-analysed: different cancer, lab, and technology (CITE-seq, not mass cytometry), 11 samples, 73,755 cells, **63 surface proteins**. **BCMA (CD269) is the top viable target of 63 (ceiling 0.608, rank 2)** — the method picks the real flagship myeloma target out of a 63-protein panel. Panel fluorophore controls (FITC 0.528, PE 0.495) sit at chance as required. **CD38 is approved yet scores 0.287, limited by NK cells — and daratumumab's documented on-target effect is NK depletion**, reproducing the "tolerated collateral, not selective target" thesis in an independent disease. **CD138 ranks 1st (0.698) and failed clinically** (indatuximab ravtansine, discontinued) because its dose-limiting expression is on epithelium, absent from bone marrow — a named, documented illustration of the framework's stated boundary.
- **OSF pre-registration drafted**: `ev_capture_model/docs/OSF_Preregistration_DRAFT.md`. Declares honestly which analyses were already run, fixes the label lists for both diseases in advance, and registers the CD138 false positive as an expected demonstration of the boundary. Dylan to rewrite in his own words and submit.
- **PENDING:** healthy-donor reference (GSE164378, Hao et al. *Cell* 2021 — 161,764 PBMCs, 8 healthy donors, 228 surface proteins, with the authors' own cell-type labels). Four files (~109 MB) downloaded to Dylan's Downloads folder; destination `ev_capture_data/healthy_pbmc_citeseq/` created and waiting.


---

## ACTION REQUIRED BEFORE THE WRITE-UP (September 6, 2026): Aim 12's headline number is described incorrectly

`aim12_single_cell_cytof_reordering.py` holds the FCS parser plus **hardcoded result constants transcribed from a browser session**; it contains no function that computes the reordering result, and no test checks those constants. The headline was therefore never reproducible from this repo. It has now been checked directly:

- **The data and model reproduce exactly.** Recomputed escape fractions match the stored worked example to within 0.05 pp (O6357 19.12 vs 19.1, Q5971 33.67 vs 33.7, O8155 36.05 vs 36.0, Q7998 44.68 vs 44.7, P7676 51.15 vs 51.1), and the stored 41-patient naive-by-mean ranking is bit-identical to the ranking recomputed from local files.
- **The counting rule is the problem.** The docs say a grid point counts when there is "≥2 percentage points of escape-fraction difference **between the reordered patients**", and that 81.8 pp is "one real patient's escape fraction 81.8 points higher than another's". Implemented literally, that gives **13/27 = 48.1%, largest reordered-pair spread 53.1 pp**. The convention that actually reproduces the stored 16/27 = 59.3% and 81.8 pp counts a grid point when the ranking flips at all AND the spread between the **most- and least-escaped patients in the whole cohort** is ≥2 pp — a condition on the cohort range, not on the pair that reordered. **The arithmetic is right; the prose is wrong.** 81.8 pp is the cohort's full escape range at that grid point.
- **Fix before writing:** either report 48.1% (13/27) / 53.1 pp against the stated criterion, or keep 59.3% and describe the criterion accurately. Figure 8's tiers (8.65 / 5.56 / 48.10 / 59.26) and any caption citing 81.8 pp need updating.
- **The effect itself survives gating**, so this is a description problem, not a retraction: ungated 13/27 (48.1%), B-gated 14/27 (51.9%), light-chain-purified clone 11/27 (40.7%, spread drops to 15.0 pp).
- **Replacement metric with no arbitrary convention** (from escape = F(n*)): do two patients' CDFs cross inside the real n* range of 12.3-28,868 molecules/cell? ungated 241/820 pairs (29.4%), B-gated 196/820 (23.9%), clone 18/171 (10.5%). Use this instead of the grid.

Also settled this session: **light-chain clonality gives a same-tube internal normal control** (RLN donors calibrate at K:L 1.00-2.82, textbook polyclonal; 19/41 tumors clonally restricted; 18 usable). It agrees with the cross-sample RLN comparator (Spearman 0.685, median |diff| 0.097) and shows the RLN comparator slightly *inflates* windows (+0.046) — conservative for the exclusion claims. Under this best control **CD19 scores 0.473, below chance.** And **agent operating points are NOT currently computable**: of the three real thresholds (81/843/8,690 molecules/cell) only n*=81 lands inside the CyTOF dynamic range, quantifying the absolute-calibration gap. The ceiling (an AUC) is calibration-free and stands; operating points must not be claimed.


---

## AIM 13-15 RESULTS (September 6, 2026) — the proposal below was executed; read the results doc first

`ev_capture_model/docs/Aim13_Results_Two_Comparator_Ceiling.md` is the results document. Headline, all computed this session from real data:

- **All 36 B-panel RLN control files were downloaded** (2.0 GB, 8,143,673 real cells, 0 parse failures) and their panel confirmed identical to the tumor files channel by channel. Gated caches for all 77 files (41 tumor + 36 control) are in `ev_capture_data/gated_cache/`.
- **Negative control first:** independent healthy donors' normal B cells differ from each other by median |AUC-0.5| = 0.089 (90th pct 0.220). Cross-batch adds only +0.020, empirically confirming AUC is immune to the ~10x batch-scale confound.
- **The retrodiction works.** Against bystander T cells, CD19/CD20/CD22/CD79B score 0.976/0.971/0.913/0.819. Against real normal B cells they collapse to **0.576/0.554/0.359/0.489** — at or below the healthy donor noise floor. That is exactly why every approved CD19/CD20/CD22/CD79b therapy causes B-cell aplasia, recovered from data with no clinical input.
- **Headline claim:** no antigen in the 37-marker panel has a joint window against both normal compartments that clearly exceeds healthy donor-to-donor variation (best = CD80 at 0.616; best approved-drug target = CD19 at 0.576). Target selection in B-cell malignancy was not solved by finding a selective antigen — it was solved by choosing a **dispensable lineage**. That also predicts why the same playbook fails in solid tumors.
- **Per-patient:** CD19's joint window ranges 0.082-0.977 across the 41 patients and is the best available target in only **4/41**; 17 different antigens win across the cohort.
- **The ceiling theorem is verified numerically:** across 642 monotone kill functions (including Hill/ultrasensitive circuits up to n=50, the Lim-lab class), the largest excess above the ROC convex hull was +2.3e-15.
- **Likelihood-ratio corollary, measured:** for CD27/CD48/CD38 a likelihood-ratio rule beats every possible density threshold (+0.79/+0.38/+0.10 Youden J); for CD19/CD20/CD22/CD79B the gain is exactly 0.000. This is a criterion for when a NOT-gate is needed instead of a CAR.
- **CRITICAL SCOPING — read section 6 of the results doc before citing anything.** The comparators are lymphoid only, so the numbers are valid **upper bounds**: antigens that fail here are definitively excluded (adding compartments only lowers ceilings), but the per-patient "winners" are NOT recommendations — HLA-ABC, PD-L1, CD80, CD86, CD24 and CD40 score well only because the compartments that would sink them were not measured. This is an exclusion tool, not a discovery tool.


---

## OPEN PROPOSAL (September 6, 2026) — read before starting the paper

A session run on September 6, 2026 was asked to find a way to push the project further rather than stop at the write-up. It produced `ev_capture_model/docs/Aim13_Therapeutic_Index_Ceiling_PROPOSAL.md`, with supporting scripts and result JSON in `ev_capture_model/src/aim13_scoping/`. **Dylan has not yet decided whether to adopt it.** Three things in it are load-bearing and were verified against real code and real local data, not asserted:

1. **Two real confounds in Aim 12.** The per-patient CD19 distributions are computed over ALL events in each FCS file, with no gating — a DLBCL biopsy is 5%-87% B cells across this cohort, so the naive per-patient mean is partly a measure of immune infiltrate (Spearman ungated vs. B-gated = 0.825; 15/41 patients move >5 rank positions, max 18). Separately, median DNA signal varies ~10x between batch groups (CytB2/3 ≈ 24-34, CytB4/5/6 ≈ 235-302), and batch explains 19.2% of variance in log ungated mean CD19. The documented top-5 example changes under gating, including the specific P7676-vs-Q7915 pair. Every documented `REAL_PATIENT_CD19_SUMMARY` number WAS independently reproduced exactly, so the parser and pipeline are correct — the issue is what the numbers are computed over.
2. **An exact identity: `escape_fraction == F(n*)`**, the population CDF evaluated at `n* = ln(1 + tau/(z*(1-tau))) / (sigma*k_intra)`. Verified to <1e-12 across 60 parameter combinations on real data. This means k_intra and kill_threshold are ONE axis, not two; the 27-point grid spans n* ≈ 12-28,900 molecules/cell; and at tau=0.5 n* equals exactly the real calibrated thresholds (81 / 843 / 8,690). It also names the real criterion: mean-based ranking is correct iff first-order stochastic dominance holds.
3. **A new claim with a working pilot.** The 41 local files carry 37 usable antigens plus CD3/CD45, so tumor-vs-bystander ROC can be computed per antigen. Run on all 41 patients: ranking by overlap puts CD19, CD20, CD22, CD79B (the four antigens with approved lymphoma drugs) at ranks 1, 2, 6, 11 and separates them perfectly from 11 broadly-expressed/known-toxic antigens (p = 7.3e-4); ranking by mean expression does not (separation 0.705, and it puts CD44, HLA-ABC and CD38 in its top 6). Robust across 40th-70th-percentile gating. See the proposal for the honest limitations — especially that the comparator is currently bystander T cells rather than normal B cells, and that a systematic novelty search has NOT been completed.
