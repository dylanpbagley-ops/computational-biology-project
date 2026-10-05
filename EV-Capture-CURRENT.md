# EV Capture / CAR-T Pivot Project — CURRENT STATE (living document)

**This is the file to attach to a new Cowork session.** It is kept small on purpose (see `PROCESS.md`). It contains only what is still actionable. Session-by-session narrative lives in archive files, pointed to below.

**➡️ READ `WHERE-EVERYTHING-IS.md` FIRST.** It is a one-page map of every folder, dataset, script and document in this project, with what each one is for. It exists so a new session never has to hunt.

**Consolidation note (September 8, 2026):** this file had reached 88 KB against `PROCESS.md`'s 30 KB trigger — ten "what changed this session" sections had accumulated back-to-back, exactly the drift PROCESS.md warns about. All fourteen dated narrative sections were moved to `ev_capture_model/docs/archive_sessions_sept3-8_2026.md` and every still-true fact folded into the sections below. Nothing was deleted.

**Size flag (September 10, 2026):** this file is back up to ~36 KB, past the 30 KB trigger again — mainly because this session's own "what changed" section is long. Not yet consolidated (the September 10 content is still current, not settled history to fold away) — due for the same treatment as the September 8 pass once this session's facts have made their way fully into "Current state" and "Known risks" and the narrative itself is no longer needed inline.

---

## Who I am and how to work with me

High school student, working solo, ISEF-level science fair project in **computational biology**. Regional fairs are typically February–March. Roughly six months out from early September 2026.

- **Plain language, minimal jargon.** Define terms the first time. Unfamiliar vocabulary hinders my understanding, doesn't impress me.
- **Explain fully.** Concise is good, but not at the cost of leaving something out.
- **I'm a beginner coder.** Python basics only. Walk me through both the physics and the debugging.
- **Stress-test before executing.** Evaluate from multiple angles and flag problems before implementing.
- **Be honest about weaknesses.** "This could fail because X" now beats finding out in January. Don't reassure me reflexively.
- **I catch errors and expect them acknowledged.** Acknowledge corrections directly, don't gloss over them.
- **When a decision genuinely needs more explanation before I can answer accurately, say so and re-ask rather than guessing.**
- **When a well-scoped search or attempt comes up empty or gives a surprising/negative result, report the honest state plainly and ask how to proceed**, rather than declaring defeat unilaterally, quietly working around it, or grinding indefinitely.
- **After editing a function's return value/keys/signature, run every affected script's own standalone entry point end-to-end, not just its unit tests.** A passing test suite is necessary but not sufficient.
- **Verify claim-by-claim against real code/data before trusting ANY prior summary** — a past session's own handoff, an external tool's critique, or your own code's first-pass output.
- **If you catch yourself narrowing a broad search to "the one thing that seemed most exciting," stop and ask whether that narrowing is actually justified.** This has happened for real multiple times in this project's history (see the archive files) — don't wait to be asked twice to broaden back out.
- **Don't treat any single document (including this one) as a fence around what's worth doing.** Read source material in full, not just its "recommendations" section, before deciding scope is exhausted.
- **Never use the built-in Cowork browser — it causes the Claude desktop app to randomly close.** Use the Claude in Chrome MCP tools instead for all browser work.
- **When a large file needs to move from a browser context onto the actual device filesystem, and in-session tricks (blob-download, base64 round-trip) stall or fail, ask Dylan to download it himself and drop it in the outputs folder — that has now reliably worked when nothing else did (see GSE193531, Sept 9).**

---

## The project in one sentence — the headline has CHANGED, read this

**Lead with the theorem, not the pipeline.** As of September 8, 2026 the project's strongest claim is a proved impossibility result, not a data analysis:

> **There is a hard ceiling on what any antigen-targeted cancer therapy can achieve, set entirely by the antigen's distribution across normal cell types — and the ceiling predicts which real drugs poisoned patients.**

Three checkable claims, in the order they should be presented:

**CLAIM 1 — the theorem (results §8, §17, §29).** For any kill rule that rises with antigen density, the achievable (tumour killed, harm) points lie in the convex hull of the ROC of the tumour distribution against a **weighted mixture** of normal cell types, where each type's weight is how much the body needs it. Verified: **0 violations in 642, 666, 1,953, 5,112 and 16,400 separate checks.** §29 proves the project's empirical rule (confinement × expendability) is the **leading term** of this theorem, so the two are one object rather than a theorem plus a heuristic.

**CLAIM 2 — the ceiling predicts real clinical toxicity (results §23, §26, §27).** Confinement across 154 normal cell types separates 19 tolerated from 10 on-target-toxic targets at **0.868, p = 0.00031**; the two-factor rule reaches **0.905, p = 0.00004**; and it holds with **every blood target removed** (solid tumours only: 0.792, exact p = 0.023). It replicates on an independent protein reference (antibody staining, 0.785, p = 0.010) and in an independent solid-tumour cohort from another lab (0.750, p = 0.046).

**CLAIM 3 — the counterintuitive one, replicated (results §22, §27).** How well a target separates tumour from normal does **not** predict whether its drug was tolerated (0.528, p = 0.44). How confined it is across normal tissue does (0.750, p = 0.046). Same targets, same labels, same statistic. This is the claim judges will disbelieve; it survived transplant to a solid tumour, an independent lab and a different technology.

**The EV foundation (Aims 1–9) is now the origin story, not the headline.** Claim: affinity capture is nonlinear in per-particle receptor count, so selectivity is governed by distribution overlap rather than fold-change; the CA19-9 / EGFR ranking flip on PANC-1 versus normal pancreatic cells is where the whole project started, and it is pancreatic, which is where §24 returns.

**Presentation risk, stated once because it dominates every ISEF outcome:** projects that open with "I analysed CAR-T data" land at Fourth Award — that is literally where 2026's ReinforCell and the PDAC multi-omics project landed. Open with the impossibility result.

---

## What changed, September 10, 2026 (OSF task walked through; panel grown to 44)

Definitive results: `ev_capture_model/docs/Aim13_Results_Two_Comparator_Ceiling.md`, now
**sections 1–35** (new "Sixth addendum").

- **Panel grown 42 → 44 targets.** Re-checked the ~30 remaining unresolved ADC trial-registry
  genes (De Schepper et al. 2025) for a verdict that has since settled. Found two: **ERBB3/HER3**
  (patritumab deruxtecan — FDA BLA voluntarily withdrawn 2026 after HERTHENA-Lung02 missed its
  overall-survival endpoint; a clean efficacy failure by the sponsors' own account) and **AXL**
  (enapotamab vedotin — Genmab discontinued 2020 for insufficient proof-of-concept efficacy, no
  toxicity concern cited). Both added as `efficacy_fail`, each verified against primary/near-primary
  sources. New composition: 19 safe, 11 toxic, **12** efficacy_fail, 2 excluded.
- **A wrinkle on the existing SLC34A2 row, investigated and documented, not acted on.** A related
  ADC on the same target (upifitamab rilsodotin) showed a bleeding-event signal, but it was never
  mechanistically confirmed as on-target, and the sponsor's own stated failure reason was efficacy.
  Label unchanged; the wrinkle is now recorded in the evidence field and in `Label_Review_Sheet.md`
  so a reviewer can look at it independently.
- **An unplanned side effect, found and reported honestly rather than glossed over: the growth
  moved the section 23.4 specificity control from borderline non-significant (0.700, p = 0.043,
  19 vs 10) to a small but real separation (0.746, p = 0.011, 19 vs 12).** Investigated in depth
  before deciding anything: the skew predates this session's two additions, has a coherent
  biological explanation (a bimodal split within efficacy_fail — a few immune-lineage high scorers
  like FCER2/CD23 and HLA-DRA versus ten broadly-expressed solid-antigen low scorers), and does
  **not** move any of sections 23.1–23.3 or 23.5–23.7 (bit-identical to their pre-edit values).
  Section 23.4 was rewritten to disclose this as a limitation of one control rather than presented
  as a clean pass — see results-doc §23.4/§34 for the full honest reading. `test_wholebody_
  restriction.py` and `verify_documented_numbers.py` both updated to guard the new numbers going
  forward. Full pytest suite for this slice of the pipeline (19 tests, 3 skipped) and a targeted
  section-23 number-by-number check against the results document both pass; the full repo-wide
  195-test suite and 22-check verifier were not re-run end-to-end this session (see Known Risks).
- **`Label_Review_Sheet.md` updated to match** (rows 38b/38c added, counts corrected, SLC34A2
  wrinkle noted) — rows 38b/38c are new since the last outside-reviewer pass and haven't been
  checked yet.
- **The OSF task was walked through this session** — two distinct sub-parts: (a) publicly
  timestamping the locked-predictions file, (b) rewriting `OSF_Preregistration_DRAFT.md` in
  Dylan's own words for formal submission at osf.io/registries. Both require Dylan's own OSF
  account and manual steps on osf.io, so completion is on him, not tracked in this repo.

Two long-open items from the September 8 assessment (§8 "next steps" #5 and #9) are now closed,
as far as public data allows either one.

- **§11a — a third real calibrated agent on the ROC.** Bluhm et al. 2018 (BCMA CAR-T,
  QuantiBRITE-calibrated, n*≈100) applied to the myeloma CITE-seq data via
  `bcma_operating_point.py`. **Honest result: the same absolute-calibration gap §11 found for
  CD19, replicated independently in a second disease and technology** — itself the reportable
  finding. Caveats and the superseded B7H3 candidate: see Known Risks and "What NOT to revisit"
  below.
- **§33 — an independent healthy bone-marrow reference for myeloma.** Closes the §16 asymmetry
  (CD19/DLBCL had one, myeloma didn't). Boiarsky et al. 2022 (GEO GSE193531, `extract_boiarsky_
  markers.py` + `boiarsky_healthy_reference.py`) implemented: BCMA detected in **82.87%** of
  9,329 genuinely healthy (NBM) cells, independent confirmation this project's own cohort isn't a
  cohort-specific artifact. Da Vià et al. 2025 is the scientifically stronger candidate but
  confirmed EGA-controlled-access, not pursued. Scope caveats: see Known Risks below.
- **Test suite still 195 passing; verifier still 22/22** (§11a/§33 are narrative findings, not
  wired into the checker, same as every prior non-tabular finding).
- **Outputs-folder cleanup, same session** (root clutter moved to proper homes; a `_to_delete/`
  folder needs Dylan to empty it himself): full list in `WHERE-EVERYTHING-IS.md`.

---

## September 6–8, 2026 (three sessions) — folded in, archived

The theorem (§29), the whole-body result (§23), the pancreatic thread (§24), the shuffled-label
null (§25), the protein cross-check (§26/§32), the solid-tumour replication (§27/§28), the
genome-wide discovery screen (§30), and the locked prospective predictions all landed in these
sessions. Their still-true facts are already folded into "The project in one sentence" (Claims
1–3), "Current state," and the results-doc section table in `WHERE-EVERYTHING-IS.md` — nothing
further to repeat here. **Full blow-by-blow, and all seven of that period's recorded self-caught
errors (a substring bug, an over-tight claim withdrawn, a wrong "third bug" diagnosis, a hull bug,
a NaN bug, a scoring inversion, and a wrong control that reversed §26's conclusion once
corrected):** `ev_capture_model/docs/archive_sessions_sept3-8_2026.md`.

---

## Current state: done vs. not done

### Done
- **Aims 1–9 (EV foundation)** — built, tested, settled. CA19-9 flip is the primary headline; superselectivity reachable across ~59.6% of the plausibility grid; identifiability proven.
- **Aims 10–12 (CAR-T pivot)** — built and run end to end. Aim 12 on the complete 41-patient DLBCL CyTOF cohort (6,865,361 cells) is the strongest of the original results.
- **Aims 13–15 (the ceiling work, §13–§33)** — the theorem, the whole-body extension, the pancreatic thread, two protein references, the solid-tumour replication, the discovery screen, the locked predictions, a third calibrated ROC agent (§11a), an independent myeloma healthy-marrow reference (§33).
- **Three real antigen calibrations** — CD19, CLDN18.2 (r = 0.964, RMSE 0.237), mesothelin.
- **All 12 write-up figures**, plus **fig9** (whole-body) and **fig10** (pancreatic thread), as 300 DPI PNG and editable-text SVG in `Poster, Paper, and Write Up/Figures/`.
- **Rigor infrastructure:** frozen label file, label review sheet, OSF pre-registration draft with outcomes recorded, shuffled-label nulls, a documented-numbers verifier, 195 tests.
- **Item 5 (a better-calibrated ROC agent) — closed.** Bluhm et al. BCMA added as a genuine
  attempt (§11a); result is an honest negative (the same calibration gap, replicated), not a
  plotted point — and that is the correct, honestly-reported outcome, not an incomplete one.
- **Item 9 (myeloma healthy-marrow asymmetry) — closed as far as public data allows.** Boiarsky
  et al. implemented (§33); Da Vià et al. verified as scientifically stronger but confirmed
  blocked by EGA controlled access, filed and not pursued further.
- **Panel growth toward 74 (September 10) — pushed from 42 to 44, and that is very likely close
  to the practical ceiling of this approach.** Systematically checked the ~30 remaining unresolved
  ADC registry genes; found exactly 2 with a settled clinical verdict (ERBB3, AXL). The rest are
  still in active trials with no verdict to label — see `Label_Review_Sheet.md` and the results
  doc's "Panel growth" notes. Getting to 74 is not realistically achievable by literature review
  alone; most of the gap is trials that haven't concluded yet, not judgements this project failed
  to make.

### Not done
- **The written ISEF paper and poster. This is the single most important remaining action.** ISEF bars generative AI from writing submitted abstract / research plan / poster text — Dylan writes the prose; Claude does structure, figures, feedback and drafts he rewrites.
- **One outside signature on `ev_capture_model/docs/Label_Review_Sheet.md`** — 44 rows now, each with a citation, written for a non-computational reviewer. Cheapest remaining credibility gain; rows 38b/38c (ERBB3, AXL) are new since the last check and specifically need a look.
- **Submit the OSF pre-registration** — walked through this session; needs Dylan's own account and final submission.
- **Re-run the full 195-test suite and the 22-check verifier end-to-end** after this session's panel-growth edits — only the affected slice (target_labels.py, wholebody_v2.py, the section-23 tests) was re-run and confirmed passing this session; the other ~175 tests covering earlier aims weren't touched by this edit but haven't been re-confirmed since.
- **ISEF forms, Adult Sponsor, AI-use log** — not assembled. **The AI-use disclosure log
  (`AI_Use_Disclosure_Log.md`) still stops at Aug 26-27 and does not cover the massive Sept 6-9
  work (sections 13-33) — a real, open, ISEF-required gap; needs Dylan's own review of the
  framing, not something Claude should write unilaterally.**
- **Patient-level antigen-density + response data** — exists but behind "available on reasonable request". Recommended AGAINST as a priority; send one email as a lottery ticket, build nothing on it.
- **`K_INTRA`, `z`, `k_geom` (EV side) remain uncalibrated placeholders**, permanently, by user decision.
- **CyTOF signal is confirmed NOT absolutely molecule-calibrated** — closed, not open.
- **The `ChatGPT_Transfer_Package.zip` + folder (~137 MB) in the outputs root** — flagged safe to
  delete since August, still there, still Dylan's call.

---

## Immediate next steps, in order

**0. Don't treat this list, or any document, as a fence around what's worth doing.** This project has repeatedly found real, missed items by reading source material in full rather than a prior list's recommendations.

1. **Write the ISEF paper.** Order: theorem (§29) → the prediction it implies (§23) → the evidence it holds (§26–28) → the locked forward predictions → the discovery screen (§30) → limitations (§23.4/§34, §26, §31, §32, and the recorded errors). See "Paper collaboration" note below for how Claude can help within the AI-authorship rule.
2. **Then the poster** — a distillation of the paper, not a separate artefact.
3. **Get the label review sheet signed by one outside person** — 44 rows now; rows 38b/38c (ERBB3, AXL) are the newest and least-checked.
4. **Finish submitting the pre-registration** at osf.io/registries, in Dylan's own words (walked through Sept 10; submission itself is Dylan's step).
5. ~~Grow the labelled panel toward 74 targets~~ — **done as far as literature review can take it** (42→44, Sept 10); the remaining gap to 74 is mostly unconcluded trials, not missed judgements. Revisit only if/when more of those trials report out.
6. Assemble ISEF forms / Adult Sponsor / AI-use log — **including catching the disclosure log up
   through the Sept 6-10 work**, which needs Dylan's own review of how to characterize it.
7. *Optional, low cost:* one email to the authors holding patient-level antigen-density cohorts.
8. *Housekeeping, whenever convenient:* empty the `_to_delete/` folder in the outputs root, and
   decide on the `ChatGPT_Transfer_Package` (zip + folder).

**ISEF odds, calibrated (see `Status_ISEF_Ratings_and_Expansion_v2.md` in the Claude project):** qualify from Arizona **85**, any CBIO Grand Award **72**, second-or-better **42**, First in CBIO **22**; top awards — Horvitz 8, Leighton & Berger 5, Coleman 5, Scripps 4, Young Scientist 3, Yancopoulos 2, Moore 2, Barrett 1. Anchored to AzSEF 2025 (786 students, 558 projects, 20 advance) and CBIO 2026 (2 First / 5 Second / 8 Third / 8 Fourth from ~85 projects).

Timeline fits — roughly five months to a February–March regional fair.

---

## Known risks

- **The EGFR healthy/PDAC-plasma flip does not survive strict pixel-exact-only filtering — retired.** CA19-9 is the primary EV-side headline case.
- **The Gillespie cooperative model does not reproduce k=3's selectivity optimum with genuine separation, for a precisely-identified physical reason** — not evidence the original paper is wrong.
- **`K_INTRA`, `z`, and `k_geom` (EV side) remain permanently uncalibrated** — impact is narrow (point predictions unavailable; range predictions, everything this project reports, unaffected).
- **A passing test suite does not guarantee a changed function's existing callers still work** — always run the actual affected script's `__main__` entry point after any change.
- **Superselectivity reachability, ranking-flip robustness, the identifiability proof, the stem-cell crossover's robustness are all SETTLED — do not re-litigate.**
- **The AND-gate finding (EV side) is built on a reconstructed, not real joint, EGFR/CA19-9 distribution.** The CAR-T AND-gate extension has the same limitation.
- **Aim 6's tetraspanin validation covers TETRASPANINS, not the project's original tumor-selectivity markers.** Never blur into "the tumor-selectivity claims are validated."
- **Aim 9's Check 1 finding (8.7% reordering) is a sensitivity/reordering finding, NOT a validation of Guo et al.'s AUC ranking** — only ITGAM has both a real AUC and real positivity data.
- **Aim 11/12's reordering findings are sensitivity/reordering findings, NOT validated clinical outcome predictions.**
- **Aim 12's CyTOF signal is not confirmed molecule-calibrated to the same scale as Aim 10's dose-response data** — read as a shape-sensitivity test, not an absolutely-calibrated clinical prediction.
- **The two prediction-registry checks already run (Predictions 1 and 3) both produced honest, non-clean-win results** — report them that way, both are "inconclusive with current data," not "false."
- **This project's own "we searched exhaustively" claims may only ever be exhaustive within one narrow category of source at a time** — always check a new search category before assuming a topic is fully closed.
- **A real, informed decision to stop pursuing a lead (OMIX010139) is a legitimate good outcome, not a failure.**
- **When an assumption produces an impossible or nonsensical result, report the impossibility directly as the finding** (Aim 9's Check 2 precedent) — don't silently clip or hide it.
- **The CLDN18.2 and mesothelin calibrations use figure-read/partial-data estimates, not full published data tables** — clearly labeled as such in the code and archive; disclose the same way in the write-up.
- **The CPTAC, BCMA/GSE234261, and pancreatic-search datasets are all real but underpowered or non-significant findings** — report as "directionally consistent, not proof of," never as confirmation.
- **The Boiarsky et al. dataset (§33) is CD138-sorted/plasma-cell-enriched, not whole bone marrow** — it can independently confirm BCMA's expression pattern within the plasma-cell compartment, but it CANNOT run a restriction-against-other-lineages test (too few other-lineage cells present). Don't describe §33 as "the myeloma equivalent of §16" without that caveat.
- **Da Vià et al. 2025 is scientifically the strongest myeloma healthy-reference candidate found (54,419 healthy plasma cells) but is EGA controlled-access** — not a quick download, needs a formal data-access-committee request if ever pursued. Not pursued this session; don't re-attempt without a genuinely new reason access has gotten easier.
- **NCBI/GEO hosts are completely unreachable from the device shell** (`curl -v` confirms a hard proxy-allowlist block on both the ftp and www hosts) — any future GEO download needs the Claude Browser pane, and if the file is large (tens of MB), asking Dylan to download it himself and drop it in the outputs folder is the most reliable path, not a last resort.
- **Section 23.4 (the specificity control) is now a disclosed limitation, not a clean pass** — as of the Sept 10 panel growth it separates safe from efficacy_fail at 0.746, p = 0.011, not the near-chance result originally predicted. This does NOT threaten the safe-vs-toxic headline (23.1–23.3, 23.5–23.7 are numerically unchanged) but it does mean the specificity claim needs to be stated honestly in the paper, not glossed over. Read §23.4/§34 before writing that part.
- **`device_bash`'s connected-folder mount was broken for an entire session (Sept 10)** — errored with a Plan9-drive-share failure, confirmed persistent by retrying twice. If this recurs, fall back to `device_list_dir`/`device_stage_files`/`device_commit_files` and do Python analysis work in the cloud workspace by reconstructing the needed subset of the repo tree there, rather than assuming the device shell is available.

---

## New risks added September 8, 2026

- **The measure is NOT technology-independent.** RNA-derived and antibody-derived confinement agree at only rho = 0.167 even at matched cell-type resolution (§32). Write "RNA-derived confinement" throughout the paper. Do NOT claim confinement is a property of the protein for which RNA is a proxy.
- **Small n is the project's real ceiling.** 19 vs 10, 8 vs 9, 9 vs 8 across the three outcome tests. The autopsy result (p = 0.046) is one relabelled target from crossing 0.05. Growing the labelled panel toward the registry's 74 targets is the highest value-per-hour work left.
- **The label file is the load-bearing soft spot.** Everything downstream rests on 40 judgements made by one person. `Label_Review_Sheet.md` exists to be signed; until it is, this is the first thing a hostile reviewer should attack.
- **The bulk-tumour contamination trap (§30.1) generalises.** Any future analysis that uses bulk tumour expression will re-derive blood-lineage genes as "cancer-elevated" because tumours contain blood and immune infiltrate. Use cancer cell lines, or a deconvolution, whenever tumour expression is needed.

---

## What NOT to revisit

Settled. Reopening these wastes time. Full reasoning for each, where not already summarized above, is in the archive files.

- ❌ Calibrating K_INTRA, z, or k_geom against any cross-marker ratio, or the calibration search itself in any form — permanently closed (EV side).
- ❌ Citing the EGFR healthy/PDAC-plasma flip as a second robust demonstration.
- ❌ Re-deriving `alpha_exact()` or the identifiability arguments from scratch.
- ❌ Modeling receptor counts as Poisson or Gaussian for real-data work — negative binomial is correct.
- ❌ Retuning the Gillespie model's `rho_n` until k=3 happens to win.
- ❌ Citing the Aim 3b combination-gain number as a single point value.
- ❌ Using GPC1 as a ranked target, or any "density category" numbers.
- ❌ The Mizenko Zenodo dataset (96.5GB of raw images) — checked, not usable.
- ❌ Re-testing the stem-cell crossover or superselectivity-reachability against yet another distribution family absent a genuinely new principled reason.
- ❌ Re-litigating the pivot decision itself.
- ❌ Treating the tetraspanin validation as validating the tumor-selectivity claims.
- ❌ Re-searching SVAtlas/PubMed/bioRxiv for a THIRD summary-statistic joint-distribution dataset beyond MASEV/CRC.
- ❌ Assuming `svatlas.org` (bare domain) works — broken TLS cert; use `www.svatlas.org`.
- ❌ Re-surveying Vesiclepedia, ExoCarta, EVpedia, MIFlowCyt-EV, EV-TRACK, or TcEVdb as raw single-vesicle data sources — all confirmed unusable, each for a specific, different, confirmed reason.
- ❌ Re-investigating the NanoFlow Repository's endometriosis dataset — confirmed a dead end.
- ❌ Resuming the OMIX010139 pursuit without a genuinely new reason to believe access has become easier.
- ❌ Re-attempting to read the SPR/QCM Langmuir paper (DOI 10.1021/acs.langmuir.5b01386) or its two related paywalled papers — confirmed no accessible full text exists anywhere legitimate.
- ❌ Assuming "no wet lab data" is this project's primary competitiveness gap — the real gap was "no genuinely new claim" (substantially addressed by Aims 9/11/12) and "no outcome-linked validation."
- ❌ Re-attempting a 4-marker ranking-correlation check against Guo et al.'s full AUC ranking without first finding genuinely new EPS8/CD9 positivity data.
- ❌ Transferring the ITGAM+/ITGAL+/ITGB2+ cluster's healthy:CRC decline ratio onto any single member marker's own individual positivity fraction — confirmed mathematically invalid.
- ❌ Re-implementing a log-log-slope "selectivity index" as if new — it already exists (`selectivity_alpha()`, `alpha_exact()`).
- ❌ Treating Predictions 1 and 3 in the registry as settled/concluded either way.
- ❌ Pursuing Option E (a public prospective prediction) against either the eLife pancreatic-proteomics paper or the CRC supermeres paper specifically — both already fully published.
- ❌ **Assuming the Mendeley 41-patient CyTOF cohort and the 2021 ASH abstract's 54-patient CyTOF cohort are the same patients just because both are Weng-lab CyTOF DLBCL work.** Checked directly, found to be an open, unconfirmed mismatch.
- ❌ **Describing the ASH abstract's CD57+/ABC outcome classifier as applicable to the Aim 12 CD19-reordering finding.** It's a T-cell-phenotype axis, not a tumor CD19-density axis.
- ❌ Re-checking the Jackson lab PDAC IMC dataset (bioRxiv 2025.08.15.670565) for MSLN, MUC1, or CEACAM5 — all three directly confirmed absent from the paper's full text.
- ❌ Re-pursuing the Kamal et al. 2025 myeloma IMC dataset (BPDC000148) or the Penn huCART-meso trial's scRNA-seq data — both ruled out for specific, confirmed reasons (wrong modality/wrong cell population respectively). See pancreatic-search archive.
- ❌ Re-attempting FlowRepository FR-FCM-Z5YT (Tislevoll AML CyTOF) downloads without checking if the site's infrastructure outage has resolved first.
- ❌ Citing an Aim 11 "17.8% broader synthetic sweep" number — confirmed unreproducible from any current code in a prior session, and removed from Figure 6's caption this session. Only the verified 5.6% Aim 11 number should be cited going forward.
- ❌ **Re-pursuing the B7H3 CAR-T functional-avidity paper (Nature Communications 2025) as a calibrated ROC agent** — confirmed too approximate (binder-dependent threshold, no single universal molecules/cell number), superseded by Bluhm et al. (§11a), which cleared the bar. Don't re-litigate; Bluhm was the better-matched find.
- ❌ **Re-pursuing Zhang et al. 2024 (GSE245108) as the myeloma healthy-marrow reference** — confirmed too thin (135 plasma cells of 72,179 total cells; the atlas targets CD34+ progenitors, plasma cells are incidental passengers). Boiarsky (§33) is the implemented answer; Zhang remains on file only as "has real protein/ADT data, just too few cells."
- ❌ **Re-pursuing Duan et al. 2023 (GSE230705) as the myeloma healthy-marrow reference** — real, large (17,347 plasma cells), but the protein/molecule-density layer (BCMA/CD138/CD38 flow data) is not in the downloadable GEO deposit, only in the paper's own summary figures. Superseded by Boiarsky for the same reason Zhang was: Boiarsky is the one candidate that's both large AND freely downloadable, even though it's RNA-only like Duan.

---

## Key literature (table)

| System | Data quality | Notes |
|---|---|---|
| PANC-1 vs. HPDEC (pancreatic) | Best matched pair — digitized, range-swept | CA19-9 is the primary EV-side headline finding. |
| Stem cells (AD/UC/BM) | Good — already range-swept | CD81/CD63/CD73 crossover, 4-seed-verified. |
| MASEV tetraspanins (Figure 4, summary) | Real per-vesicle joint SUMMARY data — validated (Aim 6) | TVD 0.011/0.019. |
| MASEV tetraspanins (Figure 6A, raw) | Real RAW per-vesicle data — directly measured (Aim 8A) | Correlation ≈0.08-0.15. |
| CRC plasma exosomes (Table S7, summary) | Real per-cohort cluster SUMMARY data — replication (Aim 7) | cluster11 replicates as CRC-enriched (2.08x). |
| CRC plasma exosomes (raw file) | Real RAW per-vesicle data — directly measured (Aim 8B) | Triple-positive fraction 0.5% (n=120/24,000). |
| CRC individual-marker AUCs (Guo et al.) | Real healthy-vs-CRC AUCs, 4 markers — ground truth for Aim 9 | ITGAM=0.9158 > ITGB3=0.8681 > EPS8=0.8498 > CD9=0.8242. |
| LBCL/myeloma CD19 (Aim 11) | Real published summary statistics (median/IQR) | Pre-treatment vs. relapse LBCL cohort; myeloma cross-disease comparison. |
| **Weng lab DLBCL CyTOF, 41 patients (Aim 12)** | **Real single-cell data, complete cohort, 6.87M real cells** | The project's strongest result: 59.3% meaningful reordering. No confirmed outcome linkage. |
| Ovarian cancer PBA, OMIX010139 | Real per-vesicle data — Controlled-access, DEPRIORITIZED | Settled — see "What NOT to revisit." |
| CLDN18.2 (Barrett et al. 2024 / AZD6422) | Real quantitative-FC dose-response, 8 cell lines (3 pancreatic) | Spearman r=0.964, k_intra fit RMSE 0.24. Real saturation/non-monotonic finding in the 3 pancreatic lines. |
| Mesothelin (Chu et al. 2024) | Real quantitative-FC titration, 2 clean anchor points | k_intra=1.31e-05. Same assay family as CLDN18.2. |
| CPTAC PDAC (Cao et al. 2021) | Real per-patient bulk protein + survival, n=135 | Directionally consistent (higher MSLN → worse outcome), not significant. Bulk data — can't feed reordering analysis. |
| BCMA CAR-T myeloma (GEO GSE234261) | Real per-cell CITE-seq + CR/non-CR outcome, n=10 patients | Directionally consistent, not significant at n=5 vs 5. Single-point calibration anchor found (222 molecules/cell, bb2121). |
| **BCMA CAR-T activation threshold (Bluhm et al. 2018)** | **Real QuantiBRITE-calibrated absolute-molecule threshold** | n*≈100 molecules/cell. Applied to this project's myeloma CITE-seq data in §11a; found the same calibration gap as CD19 — an honest negative, not a plotted ROC point. |
| **Boiarsky et al. 2022 healthy bone marrow (GEO GSE193531)** | **Real scRNA-seq, 9 genuinely healthy (NBM) donors, CD138-sorted** | 29,387 total cells, 9,329 from healthy donors. RNA-only (no protein layer). Used in §33: BCMA detected in 82.87% of healthy NBM cells, independent confirmation of this project's own myeloma cohort's BCMA behavior. |
| **Da Vià et al. 2025 healthy/malignant plasma cells (EGA EGAC50000000232)** | **Real scRNA-seq, 21 donors, 54,419 healthy plasma cells — scientifically the strongest myeloma healthy-reference candidate found** | **EGA controlled-access — confirmed, not a quick download.** Filed, not implemented. Revisit only if a data-access-committee request is ever pursued. |

Anchor paper (EV core model): **Mizenko et al. 2021**, *J Nanobiotechnology* 19:250, DOI 10.1186/s12951-021-00987-1.

**MASEV source:** Spitzberg, Ferguson, Yang, Peterson, Carlson & Weissleder, *Nat Commun* 14:1239 (2023), DOI 10.1038/s41467-023-36932-z, PMID 36870999.

**SVAtlas** (led to both MASEV and CRC raw data): *Nucleic Acids Research* (2025), DOI 10.1093/nar/gkaf1189, `https://www.svatlas.org` (`www.` required).

**CRC source:** Guo et al., PMID 37040507, DOI 10.34133/research.0041, journal *Research*.

**CAR-T calibration source (CD19):** Agilent/Antibody Analytics IndEx-2 application note, June 2024 (real titrated CD19 dose-response, Quantibrite PE beads).

**CAR-T calibration source (CLDN18.2):** Barrett et al. 2024, *Clinical Cancer Research* 30(23):5413-5429, DOI 10.1158/1078-0432.CCR-24-1853, PMC11609629.

**CAR-T calibration source (mesothelin):** Chu et al. 2024, *Cytotherapy* 26(4):325-333, DOI 10.1016/j.jcyt.2024.01.011.

**CAR-T calibration source (BCMA):** Bluhm et al. 2018, *Molecular Therapy*, DOI 10.1016/j.ymthe.2018.06.012, PMID 30078440.

**Aim 12 source:** Weng lab (BC Cancer Agency / UBC), "Mass Cytometry of Diffuse Large B-Cell Lymphoma (CyTOF cohort)," Mendeley Data, DOI 10.17632/rstwmcw94x.1 (Oct 2025, CC BY 4.0).

**Likely companion abstract (outcome data, unresolved cohort match):** Nissen, Wang, Sarkozy, Jiang, et al. (Weng, senior author), "Immune Profiling of Diagnostic DLBCL Biopsies Dramatically Improves upon Cell-of-Origin Risk Stratification," ASH Annual Meeting 2021 (Blood, conference abstract; 54-patient CyTOF/discovery cohort + 129-patient flow cytometry validation cohort).

**CPTAC PDAC source:** Cao et al., *Cell* 2021; cBioPortal study ID `paad_cptac_2021`.

**BCMA CAR-T (GSE234261) source:** Rade, Grieb, Weiss et al. 2024, *Nature Cancer*, DOI 10.1038/s43018-024-00763-8, PMID 38641734; GEO accession GSE234261.

**Boiarsky et al. source:** Boiarsky et al. 2022, *Nature Communications*, DOI 10.1038/s41467-022-33944-z, PMID 36396631; GEO accession GSE193531.

**Da Vià et al. source:** Da Vià et al. 2025, *Blood*, DOI 10.1182/blood.2024025643, PMID 40009503; EGA accession EGAC50000000232.

Full literature table with all citations, DOIs, and methodology quotes: `EV-Capture-Project-Handoff-v15.md` (oldest full archive). Full narrative history: the archive files listed in the file index below.

---

## Prior sessions' work — pointers only

- **September 9, 2026 (item 5 and item 9 closed):** `ev_capture_model/docs/archive_session_sept9_item5_item9.md`, and `ev_capture_model/docs/Work_Completed_Sept9.md` (Rounds 2–4).
- **September 3–8, 2026 (the ceiling work, ten sessions' worth):** `ev_capture_model/docs/archive_sessions_sept3-8_2026.md`. Contains every "what changed" narrative, all seven recorded errors with their diagnoses, and the reasoning behind each correction.
- **The CAR-T pivot and Aim 12:** `ev_capture_model/docs/ev_capture_model_narrative_archive_pivot_and_aim12.md`.
- **The pancreatic-cancer antigen search:** `ev_capture_model/docs/ev_capture_model_narrative_archive_pancreatic_search.md`.
- **The write-up figure sessions:** `ev_capture_model/docs/ev_capture_model_narrative_archive_writeup_figures_sessions.md`.
- **The novelty diagnosis and its follow-through:** `ev_capture_model/docs/ev_capture_model_narrative_archive_diagnosis_followthrough.md`.
- **Handoff versions v3–v15:** `archive/`. Historical record; never read by default.
