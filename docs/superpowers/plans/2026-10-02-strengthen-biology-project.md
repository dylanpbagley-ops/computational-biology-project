# Existing Biology Project: Strengthening and Extension Implementation Plan

> **For agentic workers:** Use `superpowers:executing-plans` if available to execute this plan task by task. Use the existing model/session, without spawning subagents unless Dylan separately requests them. Steps use checkboxes for tracking. Research gates override any generic instruction to implement every possible feature.

**Goal:** Preserve Dylan's existing computational-biology research, establish which findings are defensible, and identify and execute a feasible, substantive extension with a clear contribution beyond the repaired original result.

**Architecture:** Two coordinated tracks: evidence repair and evaluation of the existing work; bounded feasibility and novelty tests for substantive extensions. Run both tracks during the initial work period. Select one primary extension only after inspecting actual data and closest prior work. Do not turn repeated auditing into the project's endpoint.

**Tech stack:** Windows/PowerShell, existing Python project, TSV/JSON evidence records, Markdown advisory reports, primary scientific publications. Existing scientific Python includes NumPy/SciPy; use standard libraries for evidence tables where possible. No lab, new paid service, or external contact is assumed.

**Spec:** This document incorporates Dylan's approved scope and the latest conversation decisions. Current factual project state comes from `C:/Users/dylan.bagley/Desktop/outputs/PROJECT_STATUS_Oct2.md` and the versioned strengthening reports listed below. This is an execution handoff, not student submission prose, a preregistration, or an assurance of competition/publication outcomes.

## 1. Read this first: purpose and authorization

Dylan has invested substantial time in this project. He wants to remain in biology, preferably medically relevant, and aims for ISEF qualification/awards and publication. He authorized correcting weaknesses, expanding the project and investigating deeper questions. He does not want an unrelated project or a nonbiology pivot. Preserve the work that remains valid.

The previous assistant initially emphasized repairs and validation too strongly. Dylan correctly noted that those actions might improve trustworthiness without producing enough scientific originality. This plan explicitly corrects that omission. The success condition is not merely a cleaner original AUC or another limitations report. It is a defensible contribution supported by an appropriate evaluation, while candidly acknowledging that its competitive strength cannot be guaranteed.

School laboratory access is limited and equipment unspecified. ASU access is only a possibility, not secured. Do not make either a prerequisite. Do not email reviewers, contact labs, submit documents, register or publish anything without explicit authorization. Preparing analysis and reviewer materials is authorized. No reviewer has yet been contacted or returned judgments.

When this plan is provided with an instruction to execute, proceed within this scope without asking Dylan to approve every routine step. Ask early for the actual fair/manuscript deadline and practical weekly availability if still unknown; continue independent work while waiting. Do not invent these constraints. A request for this plan alone does not authorize executing it in the plan-writing turn.

## 2. Global constraints

- Working root, abbreviated **ROOT** below: `C:/Users/dylan.bagley/Desktop/outputs`. All relative file paths in this plan resolve under ROOT. Do not confuse this project with another repository.
- Keep the original label files, predictions, numerical outputs and historical reports intact. Write amended evidence, analysis and decisions into a new directory.
- Existing clinical score: larger values predict historically favorable outcomes. It is not a calibrated probability or a mechanistic law.
- Do not choose labels, source exclusions, sample-size cutoffs, uncertainty ranges or the primary endpoint to improve AUC. Record changes and their evidence-based rationale before computing a new endpoint's performance.
- Unknown mechanism is not absence of target-mediated injury. Approval, manageable effects and reversible depletion are not proof of no significant harm.
- Read historical instructions for settled research boundaries, but current user authorization and current corrected scientific scope take precedence over older contradictory statements.
- Stay with public clinical evidence, existing noninfectious biological datasets and abstract/statistical analyses. No need to introduce pathogen research or propose biological engineering. Do not attempt to evade platform safeguards.
- Dylan writes submitted abstract, paper, poster and formal research-plan prose. Provide analyses, figures, outlines, critique and factual assistance records. Do not imply that AI work was independently done by Dylan.
- Primary literature is required for scientific facts. Check current applicable official rules only when addressing actual submission or disclosure requirements. Do not rely on a generic memory of rules or make them the focus of this research phase.
- Never call assistant review independent, blinded adjudication. An outside person's score-withheld review improves independence of labeling; it does not make previously examined data an untouched test set.
- Do not make award probabilities, guaranteed publication claims, or patient-level treatment recommendations.

## 3. Orientation and authoritative files

Read in this order, selectively rather than loading the entire project:

| File under ROOT | Why it matters |
|---|---|
| `PROJECT_STATUS_Oct2.md` | Current entry point; supersedes older headline interpretations |
| `strengthening_Oct2/Strengthening_Report.md` | Repairs, fixed-rule diagnostics, remaining issues and expansion inventory |
| `strengthening_Oct2/negative_evidence_pass5/README.md` | Latest favorable-case review and remaining gaps |
| `strengthening_Oct2/negative_evidence_pass5/endpoint_contract.md` | Working observational endpoint rules, not adjudicated outcomes |
| `strengthening_Oct2/negative_evidence_pass5/favorable_case_triage.tsv` | All 19 historical favorable cases and exact unresolved requirements |
| `strengthening_Oct2/negative_evidence_pass5/source_notes.tsv` | Six freshly checked source records, with access/scope limits |
| `strengthening_Oct2/uncertainty_pass4/README.md` | Latest numerical sensitivity, explicitly conditional |
| `strengthening_Oct2/weakness_ledger.csv` | 23 issues; some statuses predate passes 2–5 and need cross-checking |
| `strengthening_Oct2/Reviewer_Instructions.md` | How to prepare genuinely score-withheld initial review |
| `strengthening_Oct2/mathematical_scope.md` | Correct mathematical claims and assumptions |
| `EV-Capture-CURRENT.md`, sections “Known risks” and “What NOT to revisit” | Settled searches/calibration limitations; its old headline is superseded |
| `audit_followup_Oct2/Review_of_Audit_Critique_Oct2.md` | Reconciliation of the original assessment and another AI's critique |

For disputed assessment details use `independent_audit_Oct1/Independent_Assessment_Oct1.md` and `Audit_Of_The_Audit_Oct2.md`. Do not redo their entire debate.

Historical `00_START_HERE_ChatGPT.md` and `00_MANIFEST.md` are mentioned in memory. If present, inspect their orientation guidance; they were not located by the recent targeted file search. Do not block on them or invent their contents. `HANDOFF_START_HERE_Oct1.md` explicitly points to the newer current entry point.

Relevant source/analysis files:

| File or directory | Use and restrictions |
|---|---|
| `ev_capture_model/src/aim13_scoping/target_labels.py` | Operational 44 tuples in literal `TARGETS`; historical outcome rationale, not corrected truth |
| `ev_capture_model/src/aim13_scoping/target_labels_frozen.tsv` | Older 40-row artifact; do not synchronize it by overwriting |
| `strengthening_Oct2/operational_labels_44_Oct2.tsv` | Dated reconciliation of operational rows |
| `ev_capture_model/src/aim13_scoping/wholebody_v2.json` | Scores: `reference_A_cell_type[gene]['top_share']` |
| `ev_capture_model/src/aim13_scoping/flag_cleanup.json` | Flags: `targets[gene]['M2_flag']` |
| `ev_capture_model/src/aim13_scoping/prospective_predictions.json` | Locked historical predictions; preserve exact rule/version |
| `strengthening_Oct2/agent_evidence_index.tsv` | 166 source rows; initial place to screen extensions |
| `strengthening_Oct2/evidence_pool_with_prediction_overlap.tsv` | Target-level overlap inventory; not eligibility determination |
| `strengthening_Oct2/clinical_review_pass2/` | Coding guide, initial extraction and sources |
| `strengthening_Oct2/clinical_review_pass3/` | Original-panel source coverage, extraction and remaining work |
| `strengthening_Oct2/negative_evidence_pass5/review_form_score_withheld.tsv` | Newest 30-target blank review template; expand into event/cohort rows |
| `strengthening_Oct2/diagnostics.py` and `diagnostic_results.json` | Baseline and modality diagnostics on historical labels |
| `strengthening_Oct2/frontiers.py`, `test_frontiers.py` | Reusable finite-distribution mathematics; inspect assumptions before reuse |
| `strengthening_Oct2/uncertainty_pass4/uncertainty.py`, `test_uncertainty.py` | AUC, exact arbitrary-label-flip bounds and zero-event illustration |
| `ev_capture_model/docs/Aim13_Results_Two_Comparator_Ceiling.md` | Detailed historical results with corrections; not all old claims remain valid |
| `ev_capture_model/docs/ISEF_Paper_Skeleton.md` | Current advisory outline; update only after substantive decisions |

## 4. Verified starting point and claims that must not reappear

The current clinical result is a retrospective association in a selected target panel. The empirical rule is `top_share + M2_flag`. The operational table contains 30 outcome targets: 19 favorable and 11 unfavorable, plus 12 efficacy failures and 2 excluded rows. Its older 40-row TSV omitted AXL, ERBB3, F3 and GPNMB. The corrected export preserves that history.

On original labels, rule AUC is 0.913875598086; specificity alone is 0.880382775120; blood flag alone is 0.710526315789. Higher score means favorable. AUC is pairwise ranking, not percent clinical accuracy. At fixed 154 cell types, tau is `(154 - 1/top_share)/153`, so tau and top_share rank identically. “A new specificity index” is not a defensible novelty claim. Never substitute a different historical lineage flag into locked predictions.

Pass 4 excluded FAP and hypothetically made GPRC5D and NECTIN4 unfavorable: AUC 0.901960784314. Across 1,024 assignments of ten disputed labels, the range is 0.516666666667–0.924242424242. The post-hoc 256-assignment subset holding CD47 and NECTIN4 unfavorable gives 0.642857142857–0.924242424242. These are mathematical scenario ranges, not confidence intervals. The other labels were held fixed; pass 5 shows uncertainty there too. Do not call either range comprehensive clinical uncertainty.

Important factual/interpretive corrections already established:

- FAP's cited injury evidence is in mice, not a human outcome.
- The Parkhurst CEACAM5 report is TCR-T, not CAR-T; a separate later CAR-T report exists and must not be conflated with it.
- GPRC5D's talquetamab report includes a grade 3 rash DLT; its historical blanket no-DLT rationale is wrong. Event-specific mechanism is not thereby proven.
- NECTIN4 skin injury and CD47 anemia illustrate why severe on-target injury and protocol DLT are different endpoints.
- FOLR1's cited corneal on-target explanation is unsupported by the inspected mechanistic evidence; that does not establish global safety.
- CD33's final ALFA-0701 report includes thrombocytopenia-related treatment cessation. Recovery alone cannot serve as a favorable endpoint.
- CD19 depletion with replacement care is not absence of normal-cell effects. A separate BCMA responder cohort cannot silently replace the historical ide-cel/cilta-cel evidence.
- MSLN supplies a context-specific explicit nonobservation; FOLH1 supplies a narrow zero for severe dry mouth. Neither automatically establishes gene-level absence of all injury.
- The stated CD45 exclusion reason is contradicted by human clinical research. Review its eligibility; do not simply add it as favorable.

The autopsy comparison reuses clinical labels and has uncertain paired superiority; it is not independent clinical replication. The small prior external-target evaluation and locked-prediction evaluations were not clean successes. Use their original records before quoting numbers. A newly found paper about an already known outcome does not become forward validation.

Mathematical scope: with a common response and specified weights, a weighted normal-reference mixture identity is valid. Monotone responses and arbitrary responses have different feasible frontiers. AUC does not determine the full frontier. Separate tissue optima need not be jointly attainable. This geometry does not derive the empirical score and is not a clinically calibrated safety ceiling. Classical ROC/Neyman–Pearson/convex-optimization principles are not new merely because applied here.

Do not repeat closed calibration searches, the retired EGFR healthy/PDAC-plasma flip, or treating reconstructed joint distributions as measured joint data. Check settled boundaries before proposing any data acquisition for Track B.

## 5. Execution directory, records and interfaces

Create `ROOT/strengthening_execution_Oct2/` for new execution work. Do not create it merely while reading this plan. If it already exists, inspect its progress and resume rather than overwriting completed records.

Required structure:

```text
strengthening_execution_Oct2/
  STATUS.md                         # compact current state, next actions, blockers
  decision_log.tsv                  # time, decision, evidence, alternatives, consequences
  input_hashes.json                 # inputs before any mutation
  assistance_log.md                 # factual assistance and student verification status
  evidence/
    sources.tsv                    # publication/access/prior-exposure inventory
    studies.tsv                    # agent, cohort and study identity/context
    observations.tsv               # event and explicit-nonobservation records
    screening.tsv                  # every screened candidate and eligibility reason
    overlaps.tsv                   # shared patients/reports; uncertain overlap explicit
    endpoint_v03.md                # versioned proposal, superseding only working rules
    reviewer_packet/               # scores/history/assistant judgments withheld
    amended_analysis_labels.tsv    # create only if label eligibility gate passes
  feasibility/
    candidate_cards.md
    prior_art.tsv
    data_audit.tsv
    pilot_selection.md
    path_decision.md
    candidate_a/
    candidate_b/
    fallback/
  analysis/
    analysis_protocol.md
    frozen_specification.sha256
    model_inputs.tsv
    evaluation.py
    tests/
    results/
  delivery/
    claim_evidence_matrix.tsv
    research_contribution.md
    figures/
    student_explanation_checks.md
    reproduction.md
```

Create files when they acquire content, not a large empty scaffolding. Tables should use UTF-8, stable IDs, explicit missing values and quoted TSV serialization where needed. Do not hand-split fields on tabs in code if they can contain quoted text.

Minimum evidence schema:

- `sources.tsv`: `source_id, primary_url, DOI_or_PMID, publication_date, access_date, access_level, inspected_locations, prior_project_exposure, search_id, notes`.
- `studies.tsv`: `study_id, cohort_id, source_id, gene, antigen_isoform, agent_id, agent_name, modality, disease_context, route, combination_context, treated_n, cohort_definition, observation_window, overlap_group`.
- `observations.tsv`: `observation_id, study_id, cohort_id, source_id, event_name, normal_tissue, event_count, denominator, denominator_unit, grade, injury_cessation, death, protocol_DLT, attribution_code, attribution_location, depletion, consequences, supportive_care, recovery, ascertainment, exposure_activity, missingness_censoring, observation_state, adequacy_limits, extractor, review_status`.
- `screening.tsv`: `candidate_id, source_id, study_id, eligibility, reason, prior_use_class, decision_date`.
- `overlaps.tsv`: `record_a, record_b, relationship, evidence, overlap_group, uncertainty`.

Use `NR` for not reported, `not_extracted` when the source has not been checked, `NA` only when not applicable. None means zero. Counts need denominator units: patients, events, treated patients, responders or evaluable subgroup. Prefer one event/grade/window per row; split mixed event-family summaries before rates are analyzed.

Allowed attribution codes follow pass 2: `direct_support`, `author_probable`, `author_possible`, `alternative_attribution`, `unresolved`, `no_event_observed`. Allowed observational states follow pass 5: `event_supported`, `explicit_nonobservation`, `indeterminate`, `ineligible`. A provisional assistant decision must remain identified as such. A study can have several different states for different events.

If an amended binary analysis becomes defensible, `amended_analysis_labels.tsv` must include `analysis_unit_id, gene, cohort_ids, endpoint_version, outcome, eligibility_reason, adjudication_source, review_date, uncertainty_status`. Use explicit `favorable`/`unfavorable` names, not unexplained 0/1. Do not create universal gene truth by pooling contexts automatically.

## 6. Task 0 — Establish a safe, short startup

**Inputs:** Orientation files and existing artifacts. **Outputs:** `STATUS.md`, `input_hashes.json`, initial decision log.

- [ ] Read the orientation files and applicable local `AGENTS.md` if present. Exclude `independent_audit_Oct1/packages/` from broad searches; prior scans there produced access-denied noise.
- [ ] Confirm the root, available runtime and whether version control exists. Do not assume a Git repository or create a worktree/commit workflow for a nonrepository.
- [ ] Hash existing source labels, frozen predictions, score inputs and pass 2–5 manifests. Record current versions rather than rewriting old manifests.
- [ ] Reconcile the latest evidence status with the 23-item ledger. Keep a new concise ledger rather than “closing” every open issue because a review form now exists.
- [ ] Ask for the deadline/weekly effort if unknown. State that the initial feasibility work can proceed without an answer.
- [ ] Record the critical distinction: repair is needed; a substantive extension is not yet selected or demonstrated feasible.

**Acceptance:** A resuming worker can identify the exact next task and protected inputs without rereading the entire project. Do not rerun the entire historical test suite just for orientation.

## 7. Task 1 — Define the next clinical evidence pass once

**Inputs:** Pass 2 coding guide, pass 5 contract/triage and original-panel coverage. **Outputs:** `endpoint_v03.md`, screening frame and initial populated evidence tables.

- [ ] Carry forward the common question: clinically significant human normal-tissue injury attributed to intended-target recognition in a specified context. Keep grade >=3, cessation and death separately inspectable even if a composite is later used.
- [ ] Explicitly distinguish the broader drug-toxicity outcome from the target-recognition outcome. Do not relax attribution simply to increase sample size.
- [ ] Specify what an explicit nonobservation actually covers. Zero severe dry mouth cannot become zero all-organ injury. No DLT cannot become no severe injury.
- [ ] Separate intended depletion from resulting impairment, serious consequences and required care. Apply the same severity standards to blood and nonblood tissues.
- [ ] Define required adequacy fields without inventing a numerical cutoff: endpoint-specific denominator, monitoring, exposure, observation window, censoring/missingness and supportive care. Record where outside clinical judgment is needed.
- [ ] Use the entire original 30-target panel for evidence repair and review the 14 other operational rows for eligibility/exclusion logic. Do not treat efficacy failures as automatically negative clinical examples.
- [ ] Copy previously verified observations with their actual provenance. Mark inherited abstract-level evidence as inherited; do not relabel it a fresh full-text review.

**Acceptance:** The same written rules apply to every lineage and agent, while unresolved clinical judgments remain visible. Protocol version is dated and hashed before scoring any amended labels. This is post-discovery specification, not retroactive preregistration.

## 8. Task 2 — Complete bounded, high-value source extraction

**Inputs:** Existing source pool, triage and remaining-work files. **Outputs:** event-level records, access/search log, corrected reviewer packet.

- [ ] Work first on CD19/MS4A1/CD22/CD79B/TNFRSF17/CD38/CD33/CD52 depletion and consequences; then antibody-drug-conjugate attribution (TNFRSF8, TACSTD2, MET, NCAM1, F3, FOLR1); then remaining original uncertainties. SLAMF7, DLL3, GPC3, CLDN18, GPRC5D and the unfavorable cases still require consistent treatment.
- [ ] Retrieve full primary safety sections/supplements where abstract evidence is insufficient. Record source locations, not only links. Use publisher, PMC, PubMed and legitimate author repositories as available.
- [ ] For each unresolved fact, attempt at most three materially different retrieval routes in the initial pass. If access still fails, log precisely what is missing and continue. Do not repeatedly search the same abstract or fill gaps from a review article as if verified primary evidence.
- [ ] Split mixed-agent historical rows: CD19 tisagenlecleucel versus axi-cel; BCMA ide-cel versus cilta-cel; MSLN regional CAR-T versus other agents; EGFR/EPCAM mixed settings. Record antigen variants such as CD44v6 and CLDN18.2.
- [ ] Keep initial/final/pooled reports from the same cohort linked. ELIANA/ENSIGN pooled follow-up and ALFA-0701 updates do not add independent original cohorts. Do not transfer final-report counts into initial-report denominators.
- [ ] Complete the first bounded pass after every priority item has either resolved evidence or a concrete gap/access record. Do not wait for all clinical uncertainty to disappear before running extension feasibility.
- [ ] Prepare a reviewer packet containing primary sources, neutral row identities and blank judgment fields. Keep scores, historical labels and assistant interpretation separate. One suitable outside review is useful; two are preferable if practical. Neither has been secured.

**Acceptance:** All 30 originals have an explicit evidence state and unresolved requirement, all 14 other operational rows have a scope decision or recorded uncertainty, and priority event rows are traceable. This is extraction completeness relative to this pass, not an exhaustive literature claim or proof of safety.

**Human dependency:** If review is unavailable, continue other tasks. Do not fabricate review, stop all work, or repeatedly ask for permission to prepare materials. Ask Dylan about reviewer availability only when a concrete packet and description of the expertise needed are ready; do not send it yourself.

## 9. Task 3 — Start extension feasibility before repair becomes endless

**Timing:** Begin after Task 1, alongside Task 2. By the second substantive execution session, produce candidate cards and an actual data/novelty feasibility update. A “session” here is a progress checkpoint, not a promised number of minutes.

**Outputs:** `candidate_cards.md`, `prior_art.tsv`, `data_audit.tsv`, `pilot_selection.md`.

- [ ] Assess candidates A and B below; do not assume A is adopted because the prior assistant found it compelling.
- [ ] For each write: biological question, proposed claim, retained prior assets, closest prior work, accessible data, independent units, minimal demonstration, evaluation design, main failure mode and stopping rule.
- [ ] Screen existing local data and sources before launching large downloads. Inspect a real table or file, not just a dataset description.
- [ ] Do one focused primary-literature search per candidate, then one gap-filling follow-up if a concrete uncertainty remains. Begin with roughly 5–8 closest primary works total, not an open-ended field survey. Expand only when a particular paper introduces a material unresolved overlap.
- [ ] Search both supportive and disconfirming terms. Record exact query, date, source and the feature it overlaps. For A, start with normal-tissue expression/off-tumor safety, target-specificity methods and cross-agent clinical safety; for B, response-distribution bounds, partial identification, uncertainty in single-cell measurements and decision robustness.
- [ ] For A check the closest already named methods: SCAN-ACT, BayesTS, Thomas et al. and Madan et al.; identify exact publications from existing critique/reference records before describing them. For B check the classical ROC/likelihood-ratio/convex-mixture antecedents as well as relevant biological applications.
- [ ] Fill a comparison table: prior paper's exact claim, same/different data, same/different question, proposed additional contribution, evidence needed. “Nobody combined these exact words” is not novelty evidence.

**Acceptance:** Each candidate has a specific potential contribution and identified accessible evidence, or a concrete no-go reason. A literature discussion without a data check is not completion.

## 10. Candidate A — When is a gene-level expression score insufficient, and does context add useful information?

This is a candidate within the current project, not an established novel result. The fact that different drugs behave differently is already unsurprising. Merely displaying discordant examples or giving them a new dataset name is insufficient.

### A1. Data feasibility and fair sampling

- [ ] Use `agent_evidence_index.tsv` to group by normalized gene, agent, modality, trial and cohort. Resolve duplicates before counting. Retain dual-target/isoform ambiguities as separate scope problems.
- [ ] Create a fixed pilot sampling rule before measuring candidate performance: inspect all sources for the first up to eight multi-agent genes in alphabetical order with at least two distinct study/agent candidates; add up to four single-agent controls from distinct modality groups using a documented deterministic order. This is a workload sample, not a statistically adequate sample size or confirmation set.
- [ ] Verify that apparently different agents truly differ and that reports do not reuse patients. Record all eligible and ineligible observations, including non-discordant cases.
- [ ] Check comparability of endpoint, tissue, grade, attribution, observation opportunity, supportive care, population and context. Do not count unequal observation as biological discordance.
- [ ] Produce a table of usable comparisons, unresolved comparisons, independent targets and actual outcomes. No binary performance metric if the negative category is not defensible.

### A2. What would be a substantive contribution?

There are two possible branches; select one only after inspecting the pilot:

1. **Context-aware empirical test:** Compare the fixed gene-only rule with a small, interpretable context model or stratified analysis using pre-outcome variables available across the dataset. Candidate variables include modality, route, combination setting and antigen variant. Choose variables for biological relevance and availability, not correlation with these outcomes. Do not use toxicity management initiated after injury, cessation, approval or authors' final safety language as predictors.
2. **Quantified applicability limits:** Establish reproducible conditions under which gene-only ranking fails to distinguish outcomes, with an external comparison and uncertainty assessment. A few cherry-picked discordant studies are a case series, not a general applicability rule. Within-gene comparisons test limitations of a constant gene score; they cannot validate ranking between genes by themselves.

No automatic complex machine-learning model. A small dataset with sparse modalities may justify only stratified descriptive evidence. If fitting is feasible, define target-grouped separation before any tuning. Evaluate on targets excluded from fitting; same-target new-agent evaluation answers a different question and must be reported separately. All observations for a target stay together for unseen-target generalization. Clinical expert review is still needed for endpoint validity.

### A3. Minimal pilot deliverables and gate

- [ ] Deliver normalized studies/observations, an overlap graph or table, comparability matrix, denominator accounting, and one interpretable diagnostic figure.
- [ ] Assess whether outcome diversity and independent target count support the intended inference. Use simulation/precision planning over plausible effect sizes if necessary; never declare “enough” based solely on a convenient row-count threshold. A pilot's favorable observed effect is not an unbiased power assumption.
- [ ] Compare the proposed contribution directly with closest work. State what a positive and a negative result would each teach.
- [ ] **Proceed** only if outcome definitions are supportable, the needed context fields exist, repeated targets are handled correctly, the evaluation is feasible, and a contribution beyond “context matters” can be specified.
- [ ] **Revise/narrow** if only a descriptive case series is supportable; identify whether it can be a useful supporting component, not pretend it is automatically the stronger core.
- [ ] **Stop A as the primary extension** if comparisons remain mostly unadjudicable, useful predictors are missing/leaky, or the claimed insight substantially reproduces existing work without stronger evidence.

## 11. Candidate B — What can existing biological measurements actually establish about selectivity?

This candidate uses the project's distributional/mathematical foundation. It must answer a real biological measurement or inference question. Another abstract ROC demonstration, an uncalibrated clinical harm curve, or a repeat of settled reference sensitivity is not sufficient.

### B1. Identify one grounded application

- [ ] Inventory local distribution-level data and earlier results through `WHERE-EVERYTHING-IS.md`, current limitations and the relevant source scripts. Candidate sources may include existing EV per-particle observations or single-cell expression/protein measurements; none is automatically adequate.
- [ ] Specify the measured quantity, biological groups, independent subjects/samples, units, normalization, assay detection limits and the particular comparison/question. Inspect raw or processed rows and metadata before choosing it.
- [ ] Do not mix RNA with molecule-calibrated protein, treat cells as independent donors, use reconstructed marker joint distributions as measured ones, or restart permanently closed absolute-calibration searches.
- [ ] Choose one dominant, evidenced uncertainty: missing normal compartments, finite donor sampling, measurement censoring, or plausible reference partition changes. Several arbitrary perturbations together do not constitute a biologically identified uncertainty model.

### B2. Define an inference that could add knowledge

Example form, subject to the data: “Under measurement uncertainty supported by this dataset, is conclusion X stable, contradicted, or not identifiable?” X must be a named biological selectivity/ranking conclusion or measurement-sufficiency question, not a generic statement that heterogeneity matters.

- [ ] Define the uncertainty set using measurement information or externally justified ranges; separate empirical support from assumed stress tests. Do not pick ranges to produce a striking reversal.
- [ ] Write the estimand before coding. A range across several sampled scenarios is not a proven bound over a continuous uncertainty set.
- [ ] If using `frontiers.upper_at(p,q,budget,kind)`, retain the monotone/arbitrary distinction and common-response assumptions. Normalized distributions do not establish clinical harm weights.
- [ ] Distinguish the worst/best oracle ceiling over uncertain distributions from a single response that works under every distribution. The order of optimization matters; do not swap min and max or label an oracle result robustly attainable.
- [ ] Compare against an appropriate simple baseline such as mean-expression comparison or the existing point-estimate frontier, on the same measured scale and biological units.
- [ ] Demonstrate why the uncertainty analysis changes an interpretation or identifies a concrete missing measurement. If it only redraws the existing result, it is not yet a substantive extension.

### B3. Minimal pilot and gate

- [ ] Produce an assumption table, one real-data example, an uncertainty visualization, and a comparison with the point-estimate approach.
- [ ] Check finite examples analytically or by exhaustive enumeration; for continuous claimed bounds use a valid optimization/proof rather than a grid dressed as a theorem. Confirm attainable examples where sharpness is claimed.
- [ ] Use donor/sample-level resampling where appropriate. If only one biological sample exists, label the example descriptive and seek another existing independent sample/reference before claiming generalization.
- [ ] **Proceed** only if the biological question is real, measurements support the chosen uncertainty model, existing classical results do not already exhaust the proposed contribution, and a convincing empirical demonstration is feasible.
- [ ] **Stop B as the primary extension** if it requires inaccessible calibration, invents normal-tissue harm weights, recreates already completed sensitivity work, or produces only a mathematical restatement without a new biological finding.

## 12. Task 4 — Make a bounded path decision

**Outputs:** `feasibility/path_decision.md` with evidence links and a recommended primary direction.

Evaluate A and B against the following hard questions; do not hide a failure in a weighted numerical “novelty score”:

1. What is the exact scientific claim beyond the corrected original project?
2. What closest primary work already does part or all of it?
3. What data have actually been opened and verified as usable?
4. What are the independent units, confounders and unresolved dependencies?
5. What observation would disconfirm the proposed claim?
6. What fair baseline and evaluation could distinguish real added value?
7. What is the smallest credible completed deliverable within Dylan's actual deadline?
8. What prior work is retained, and what new work is required?

- [ ] Select the strongest supported direction, not the one with the most exciting preliminary p-value. Record provisional versus resolved clinical dependencies.
- [ ] If A and B both fail, screen at most two additional adjacent questions within the existing assets. One allowable fallback is a reusable, uncertainty-aware clinical evidence benchmark plus a comparative result; a catalog alone is not presumed sufficient. Another may come from a specific unresolved biological inference exposed by the pilots. Apply the same data and novelty gates.
- [ ] Do not automatically accept a fallback because it preserves sunk effort. Do not pivot outside biology or promise a new project without discussing that change with Dylan.
- [ ] After this bounded comparison, report one of: a feasible primary extension; a conditional extension with a named external dependency; or no demonstrated strong extension from the tested paths. Explain why, including negative findings.

**Stopping rule:** Do not continue generating new audits or candidate lists without a path decision. One initial primary-literature/data pass and one targeted gap-filling pass per candidate is the default budget; extend only for a concrete solvable gap, record that justification, and keep Dylan informed. These are workload boundaries, not shortcuts around scientific adequacy.

## 13. Task 5 — Run the repaired original evaluation when its gate passes

This can support either extension. Do not wait to select a candidate merely because independent clinical review is pending, but do not call provisional labels validated.

- [ ] Reproduce the historical AUC from untouched inputs using a direct pairwise calculation with half credit for ties. Confirm score direction and the original 19/11 count.
- [ ] Complete source/endpoint review and explicit negative eligibility before freezing `amended_analysis_labels.tsv`. If this is not possible, retain observational states and report only clearly identified conditional analyses.
- [ ] Predefine one primary endpoint/comparison, inclusion set, uncertain-case handling, direction, resampling unit and missing-data policy in `analysis_protocol.md`; hash it before execution. Disclose historical exposure and any later changes.
- [ ] Evaluate the unchanged rule, specificity alone, tau, and blood flag alone on identical eligible records. Treat tau equivalence as known, not an additional independent success.
- [ ] Assess whether a relevant existing method can be compared fairly from primary documentation and available inputs. If its data or task differs, document the mismatch; do not invent a poor proxy and claim superiority.
- [ ] For repeated studies, resample by target and retain within-target dependence. Report target count, cohort count, class balance and eligible comparisons. Single-class samples have undefined AUC, not 0.5; track and explain invalid bootstrap draws.
- [ ] Report paired differences and their uncertainty, not just separate p-values. Describe target/modality strata without pretending tiny strata are a controlled causal analysis.
- [ ] Preserve uncertain cases in declared sensitivity analyses. Avoid deletion/relabeling until the result looks strong. Label exclusion-based subsets as selected subsets.
- [ ] Audit locked predictions using exact original versions, verifiable lock/publication dates and previously known outcomes. If timestamp or independence cannot be verified, say so.

**Acceptance:** A reproducible, honestly scoped assessment of what survives. No requirement that AUC improve. A negative corrected finding can shape the extension but is not automatically a major contribution by itself.

## 14. Task 6 — Execute the selected extension with a locked evaluation

**Inputs:** Passing candidate gate and path decision. **Outputs:** complete analysis, figure set, reproducible results and contribution assessment.

- [ ] Convert the selected pilot into a branch-specific implementation specification with exact datasets, endpoint/estimand, baselines, uncertainty assumptions, split strategy and acceptance tests. Do not code both candidates into full projects.
- [ ] Before new result inspection, freeze the full candidate inclusion rule and primary analysis. Record all later departures with reasons. Previously examined data remain development/exploratory data regardless of this timestamp.
- [ ] For a predictive A branch, avoid leakage through shared targets, repeated reports and post-outcome context variables. Keep fitting/tuning entirely inside training groups; if data are too small, use descriptive evaluation and state its ceiling.
- [ ] For a B branch, implement only the evidenced uncertainty mechanism. Show analytic/exhaustive benchmark examples and a biological application using actual measured quantities. Never equate a finite scenario envelope with an unproved universal guarantee.
- [ ] Perform one independent-dataset or appropriate held-out-unit check if feasible. Distinguish new biological samples, new reference measurements, new targets and genuine future clinical outcomes; they validate different things.
- [ ] Use precision/sensitivity analysis to determine whether the conclusion is supported. A non-significant difference is not proof of equivalence; a small p-value does not resolve selection or novelty.
- [ ] Preserve negative and inconclusive results. If the full analysis fails the pilot's premise, reopen the documented decision once with the new evidence, rather than silently redefining success.

**Completion:** The extension either supports a specific new claim, demonstrates a well-supported and relevant limitation, or fails. Deliver the actual outcome, not a promise that more computation will rescue it.

## 15. Verification and reproducibility requirements

Available scientific runtime, verified present when this plan was written:

```powershell
& 'C:/Users/dylan.bagley/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' -c 'import numpy, scipy; print(numpy.__version__, scipy.__version__)'
```

Standard-library tests for the existing uncertainty implementation:

```powershell
python -m unittest discover -s strengthening_Oct2/uncertainty_pass4 -p test_uncertainty.py
```

Existing standalone verifiers write generated outputs and manifests. Do not run them casually during read-only checks; inspect their behavior, use an isolated copy where necessary, and do not overwrite preserved results just to reproduce them. Historical recorded results were 10 new initial tests, 16,400 theorem checks and 221 documented-number checks; pass 4 recorded five tests; pass 5 checked 19 historical and 39 prior-manifest files. These are prior results, not checks you have just run.

For new code, require these meaningful checks where applicable:

| Failure mode | Test/expected behavior |
|---|---|
| Score reversed or ties mishandled | Scores `[1,2]`, labels `[favorable,unfavorable]` give AUC 0; tied scores give 0.5; scores `[2,1]` give 1 |
| Invalid clinical metric | Empty/single-class inputs return an explicit unavailable result or clear error; never invent AUC |
| Rates with wrong denominators | Reject negative counts and patient counts exceeding patient denominator; distinguish event counts which may exceed patients |
| Missingness mistaken for safety | `NR` and missing mechanism cannot become zero/no-event/favorable |
| Overlap/leakage | Duplicate publication rows do not increase independent units; train and test target sets have empty intersection for unseen-target evaluation |
| Review contamination | No scores, historical labels or assistant judgments in reviewer initial-judgment fields |
| Mathematical class confusion | Reproduce known monotone/arbitrary and separate/joint counterexamples; zero normal mass and ties handled as specified |
| False rigorous bound | For small finite uncertainty problems compare claimed extrema to exhaustive enumeration; disclose continuous untested regions |
| Input mutation | Hash all protected inputs before/after; fail on an unintended change |
| Manuscript drift | Each figure/table headline traces to a result file and exact cohort/endpoint version |

Implement tests with concrete independent expected values before relying on new numerical routines. Do not write tests merely mirroring internal implementation. Run every changed standalone analysis entry point end to end, not unit tests alone. Record commands, versions, input hashes, seed, output hashes and unresolved failures in `delivery/reproduction.md`. Do not run unrelated entire legacy suites repeatedly once affected checks pass.

## 16. Task 7 — Synthesize the scientific contribution and student deliverables

- [ ] Produce `claim_evidence_matrix.tsv`: claim, evidence, independent unit, assumptions, counterevidence, uncertainty, allowed wording and disallowed overstatement.
- [ ] Write an advisory `research_contribution.md` answering: what we knew before, what this project adds, the strongest result, the best competing explanation, and what remains unestablished. Compare explicitly with closest prior work.
- [ ] Prepare only figures that explain a result: evidence inclusion/uncertainty; baseline comparison; and the selected extension's primary biological finding. Keep diagnostic plots separate.
- [ ] Update the current outline around the surviving contribution. Retain the EV/distribution work and validated mathematics where they support the narrative; do not force every historical aim into the main claim.
- [ ] Prepare student explanation checks: AUC versus accuracy; cohort versus gene independence; attribution versus severity; supportive care versus no injury; theorem assumptions; the extension's actual novelty and limitations.
- [ ] Update factual assistance records. Dylan verifies sources and writes submitted prose himself. When preparing actual submission materials, check official applicable rules and involve the Adult Sponsor; do not fabricate compliance or student ownership.
- [ ] Give a candid final assessment: materially stronger, modestly stronger, or not demonstrated stronger, with evidence and remaining dependencies. Competition and publication results are not predictable guarantees.

## 17. Progress cadence, handoff and definition of done

After each substantive session, update `STATUS.md` with completed deliverables, concrete new findings, failed leads, decisions, next two actions and blockers. Update ROOT's `PROJECT_STATUS_Oct2.md` with a short pointer once a meaningful new stage is completed; retain historical reports. Avoid another large recursive handoff document.

Use this progress report structure:

1. What scientific question was resolved or narrowed?
2. What evidence/code/output supports that statement?
3. What remains unknown, and does it prevent the next action?
4. Does this strengthen validity, novelty, evaluation, or presentation?
5. What is the next concrete action or explicit stopping decision?

Do not claim project completion merely because the plan's boxes are checked. The desired end state is:

- Existing findings are scoped correctly and reproducible.
- Outcome uncertainty is explicit; clinical review is either obtained or clearly remains a limitation.
- A bounded search has selected or rejected substantive extensions based on data and prior art.
- One selected extension has been executed and evaluated, if feasible.
- The resulting contribution is materially distinguishable from just cleaning up the original association, or the failure to achieve that is stated plainly.
- Dylan receives usable figures, evidence, reproducibility instructions and an outline he can understand and write from.

If no feasible extension emerges, the correct outcome is a clear decision report and a preserved, honestly strengthened existing project—not endless additional audits, automatic reassurance, or an unapproved complete pivot.

## 18. Suggested first execution turn

1. Complete Task 0 and the working endpoint specification in Task 1.
2. Create the A/B candidate cards and inspect the evidence index plus one real distribution-level data source. Do not wait for all of Task 2.
3. Record the pilot sampling rule and closest-prior-art questions before calculating candidate metrics.
4. Complete one bounded source-extraction batch that serves both the original repair and candidate A feasibility.
5. Report concrete initial feasibility evidence and next actions. Do not end with only “I can do this” or another generic plan.

This plan intentionally does not supply a predetermined scientific answer, invented dataset adequacy or a guaranteed novel method. Those are the empirical decisions the execution must establish.
