# Archive: Diagnosis Follow-Through, Aim 9 Build, and Open-Ended Novelty Research (Second/Third/Fourth Sessions, August 29-30, 2026)

Moved out of `EV-Capture-CURRENT.md` during a PROCESS.md-directed consolidation pass. This preserves the full narrative; `EV-Capture-CURRENT.md`'s "Current state" and "Known risks" sections carry forward every fact that's still actionable.

---

## Third session: documentation fixes, citations, and Aim 9

Instruction: fully read both diagnosis documents (`docs/Novelty_Competitiveness_Diagnosis_and_Restructuring_Options.md`, `docs/Diagnosis_Followthrough_and_Recommendation.md`) and implement their recommendations. (Later corrected in the fourth session below — this pass over-focused on the "immediate next steps" list rather than the documents' full content.)

**Two documentation fixes:** (a) corrected an overclaim that no CRC dataset has a healthy-vs-cancer grouping — true of the raw per-vesicle file, not true of Guo et al.'s own summary statistics, which the project already partially used. Added `GUO_ET_AL_MARKER_AUC_HEALTHY_VS_CRC` (ITGAM=0.9158, ITGB3=0.8681, EPS8=0.8498, CD9=0.8242) and `GUO_ET_AL_CLUSTER11_HEALTHY_TO_CANCER_DECLINE` (13.47%→1.13%) to `aim7_crc_plasma_exosome_cluster_replication.py`. (b) added Spitzberg et al.'s previously-unused capture-loss numbers (single-marker loses up to 80% of KRAS-mutant EVs; cocktails still miss 36-47%) to `aim6_masev_tetraspanin_validation.py`.

**Two citation additions to `binding_model.py`:** superselectivity literature (Dubacheva, Curk & Richter 2023) + CAR-T antigen-density dose-response literature (Majzner et al. 2020, Walker et al. 2017) as two independent cross-domain confirmations; capture-recapture ecology/statistics (Otis et al. 1978) as an independent formal-statistics confirmation of the project's core reasoning.

**Aim 9 built** (`src/aim9_marker_ranking_order_check.py`) — the project's first "new claim" check: does the model, fed real CRC marker data, predict Guo et al.'s own AUC ranking? Checked the real data situation first: Guo's AUCs exist for 4 markers (ITGAM, ITGB3, EPS8, CD9); real raw positivity data exists for 4 different-but-overlapping markers (ITGB3, ITGAM, ITGAL, ITGB2); the overlap is exactly ONE marker (ITGAM). A full 4-marker ranking check against Guo's AUCs is therefore not possible with real data in hand — reported honestly, not faked.

Two checks instead: **Check 1** (primary): does the model reorder the four real-positivity markers relative to naive linear positivity, anywhere across the project's standing plausibility grid? Real finding: yes, 9/104 grid points (8.7%). **Check 2**: a real methodological pitfall caught while first running the code — applying the real ITGAM+/ITGAL+/ITGB2+ triple-positive cluster's 13.47%→1.13% healthy:CRC decline ratio (11.9x) to ITGAM's own individual positivity (14.62%) produces 174.2%, a mathematically impossible result, because a ratio computed on a rare intersection population does not validly transfer to a single marker's own marginal probability. The first draft silently clipped this to 99.9999%; this was caught and rejected before shipping. The final version reports the impossible arithmetic directly as the finding. 13 new tests, including a locked-in regression test for the impossibility finding itself.

Test suite after this session: 149 (136 + 13 new).

## Second session: full 20-category re-investigation

The user pushed back correctly that the third session (described above, though chronologically it happened first — session numbering here follows the CURRENT.md history) had only chased 2-3 leads instead of the diagnosis document's full 20-category appendix. A full re-investigation covered all 20 categories using parallel deep-research passes. Headline results:

- **The project's own existing source papers contained substantial unused real data.** Guo et al. (Aim 7/8's CRC source) reports real individual-marker AUCs (already listed above); Spitzberg et al. (Aim 6's MASEV source) reports real capture-efficiency numbers (already listed above). This also corrected the "no CRC dataset has a healthy-vs-cancer grouping" overclaim (same fix as above).
- SPR/QCM superselectivity lead (Dubacheva, Curk & Richter, *Acc. Chem. Res.* 2023) confirmed as the strongest K_INTRA-adjacent lead, reinforced by a second independent cross-domain literature (CAR-T antigen-density dose-response).
- A new category — capture-recapture ecology/statistics (Otis et al. 1978) — supplies an independent, formally-established statistical confirmation of the project's core reasoning (mean/fold-change is a biased predictor under heterogeneous capture probability), from a field with zero connection to biology.
- A real PDAC paper directly using PANC-1 (Aim 6's own cell line) was found with genuine per-marker ranking differences; a brand-new (Jan 2026) CRC head-to-head particle-ranking paper (supermeres vs. EVs) was found.
- Several categories were honestly confirmed dead ends after direct checking: clinical trial registries, FDA submissions, NCI/population biobanks, antibody vendor Kd datasheets.

## Fourth session: full re-derivation, SPR/QCM closure, and genuinely open-ended research

**Trigger:** a direct, important user correction — the third session (per the numbering used in the file at the time) had treated the diagnosis documents' "immediate next steps" list as the complete scope, rather than independently re-deriving a complete plan from their full content. The user: "I feel like you overfocused on a few parts... what I actually wanted you to do was fully read both files and thoroughly go through everything outlined in both of those files."

**Part 1 — complete re-derivation.** Every item in both diagnosis documents was researched, not just previously-acted-on ones:
- The SPR/QCM Langmuir paper (DOI 10.1021/acs.langmuir.5b01386, the most-hyped lead in both documents) was checked directly via PubMed's copyright tool: **confirmed to have NO open-access full text anywhere** — a real, structural access barrier, not a search failure. Two related papers (DOI 10.1021/acsami.8b18843; DOI 10.1021/acs.biomac.2c01338) checked the same way: also paywalled. **Restructuring Option B (calibrate K_INTRA via this lead) is closed for this specific reason — not because the idea was wrong.**
- The PDAC integrin paper (PMC11189911) was fully read: real Western blot bulk/relative-expression data, citable as qualitative support, not usable as quantitative model input.
- The CRC "supermeres" paper (PMC12902025, Jan 2026) was fully read: real head-to-head AUC comparisons, different particle type/method than this project models, added as a citation.
- The Shi et al. 2024 CRC systematic review (48 studies, 7,727 participants) confirmed via PubMed: EV protein panels reach AUC 0.92 vs. 0.87-0.90 for individual/RNA markers — independent confirmation that marker combination beats single markers.
- All five restructuring options evaluated on their merits: A (reframe around one headline claim) — implemented. B (SPR/QCM calibration) — closed, confirmed inaccessible. C (prediction registry) — implemented (`docs/Prediction_Registry.md`). D (active marker-ranking search) — substantially addressed via Shi et al./supermeres. E (prospective public prediction) — investigated, found not currently viable for a specific reason (see Part 2).
- The Ma/Wu five-cancer dataset was evaluated as an author-contact path; a ready-to-send (NOT sent) email draft was written (`docs/Ma_Wu_Author_Contact_Email_DRAFT.md`), flagging that the recipient address was not confirmed.

**Part 2 — Option E's verdict.** A viable prospective prediction needs (1) a free, legitimate public time-stamping mechanism — confirmed to exist (OSF Registries), and (2) a specific, not-yet-released dataset that will land within the project's timeline. Checked directly: the two most promising papers found this session were both already fully published — no future release left to predict against. Prediction 4 in the registry is written as a genuinely open, registered claim ready to check IF a matching dataset appears before the fair, rather than forcing a fake "prospective" claim against already-published data.

**Part 3 — genuinely open-ended research beyond both documents:**
- Hamzeh et al., Langmuir (Dec 2025), DOI 10.1021/acs.langmuir.5c03244 — real Monte Carlo phase diagrams showing superselective onset depends on nanoparticle mechanics; an EV's lipid bilayer resembles the "deformable" category that reaches superselectivity most easily. Citation-strength, not a calibration source.
- Xie et al., PNAS (Feb 2025), DOI 10.1073/pnas.2417159122 — describes a different selectivity mechanism (receptor clustering) but its own math defends this project's simplifying assumption (independent, non-clustering EV receptors) as physically correct at EV scale, not just convenient.
- A precision correction: `selectivity_alpha()`/`alpha_exact()` ARE, by exact definition, the field's own canonical formal superselectivity metric (Martinez-Veracoechea & Frenkel, PNAS 2011), not a loose analogy — stated precisely in `binding_model.py`.
- A CTC-chip receptor-density saturation finding: >90% capture yield across a ~250-fold EpCAM density range — a boundary-condition citation showing where density-sensitivity saturates (the opposite regime from EVs).
- Two brand-new (2026) EV multi-marker cancer-detection papers (Lobb et al., Cell Reports Medicine; Tamai et al./EV-Finder) both explicitly naming single-EV resolution as a gap their own bulk methods can't address.

**Part 4 — implementation.** All findings added to `binding_model.py` and `aim9_marker_ranking_order_check.py` as citations. `docs/Prediction_Registry.md` built new: four falsifiable predictions, two already checked using data already in hand.

**Part 5 — the two prediction checks, both honest, neither a clean win.** Prediction 3 (does AUC correlate non-monotonically with positivity?): checked using ITGAM/ITGB3's real overlapping data — **ITGAM has both the higher AUC and higher positivity, monotonic, the OPPOSITE of the prediction.** Reported as a real negative result. Prediction 1 (do low-copy markers show larger reordering effects?): the four real CRC markers' reconstructed mean receptor counts clustered too tightly (3.18-4.26) to test the contrast at all — an honest "this data can't test this," not a pass or fail.

Test suite unchanged at 149 (this session was citations/documentation, not new code).

## First research session (August 29, 2026): Aim 8 raw-data upgrade, the ISEF-competitiveness question, OMIX010139

This is the session that first asked "is this project actually good enough for ISEF" and found real, published summary-statistic-only validation (Aim 6/7) was the ceiling. Found SVAtlas's raw-data download page (`www.svatlas.org/download.php` — note the `www.`, the bare domain has a broken TLS cert), retrieved two real raw per-vesicle files via a real browser tab's native download (the sandboxed shell has no route to `raw.svatlas.org` directly), and built `src/aim8_raw_vesicle_correlation_check.py`:

- **Part A (MASEV tetraspanins):** real measured pairwise correlations (PANC-1 CD9-CD63=+0.066, CD9-CD81=+0.083, CD63-CD81=+0.153, all p<0.001) from Figure 6A's raw data — a DIFFERENT dataset from Figure 4's percentages Aim 6 already validated against (caught and corrected as a real provenance mistake mid-session, verified against the paper's own text). Lower than Aim 6's committed swept value of 0.2, but positive and significant.
- **Part B (CRC integrins):** real measured triple-positive fraction (ITGAM+/ITGAL+/ITGB2+) of 0.5% (n=120/24,000) vs. Aim 7's published cluster11 fraction of 4.58% — a real, reported gap, expected since a simple threshold isn't the same computation as the paper's clustering algorithm. Also found tumor location (colon vs. rectal) significantly affects all four markers' positivity (p<0.001 each) — a genuinely new result no summary table reports. This file has NO healthy controls — confirmed directly, do not describe as healthy-vs-cancer.

13 new tests. Test suite: 136.

A self-caught scoping mistake (a "search other databases" instruction had quietly narrowed to just EGFR/CA19-9 with no justification, caught when the user asked directly) led to a genuinely broader re-search, surfacing NanoFlow Repository (real, NCI-affiliated, but its one disease-framed dataset has no confirmed matching publication — a dead end) and OMIX010139 (real, confirmed: Wu, Yang, Cai et al., *NPJ Precision Oncology* 10:65, 2026, 85 HGSOC vs. 95 healthy per-vesicle ovarian cancer data — the first real dataset anywhere in this project with a confirmed healthy-vs-cancer per-vesicle grouping).

**OMIX010139 was pursued, then deliberately deprioritized.** Registration reached, but the data-access-request form required a "Lab/Personal Homepage" and "Representative Article Url" — fields a high schooler doesn't have. A fresh re-search specifically confirmed no better/easier alternative exists (two other real single-EV cancer-vs-healthy datasets found, both similarly or more restrictively gated). The user made a reasoned, informed decision to drop it — not a search failure, not unfinished business. Do not resume without a genuinely new reason (e.g., a mentor/institution willing to sponsor account access).

A full, standalone `docs/Novelty_Competitiveness_Diagnosis_and_Restructuring_Options.md` was then produced, finding (among other things) that fresh 2024-2026 ISEF winner research shows purely computational projects DO win First Award every year checked — the real competitiveness gap is not "no wet lab data" but "no genuinely new, previously-unstated claim." This finding drove the third/fourth sessions' Aim 9 work described above.
