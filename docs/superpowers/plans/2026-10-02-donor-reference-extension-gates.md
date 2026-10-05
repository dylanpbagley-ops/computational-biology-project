# Donor Reference Extension and Contribution Gates Implementation Plan

> **For agentic workers:** Use `superpowers:executing-plans` task by task. Dylan explicitly requests one `gpt-6-luna` executor, Max reasoning, with no inherited conversation. The parent drafts this plan and performs final review; the executor performs all research, implementation and corrections. Do not spawn additional agents.

**Goal:** Test and execute a medically relevant substantive extension beyond the original retrospective association, or establish candidly why the bounded tested directions cannot deliver it.

**Architecture:** Resolve specific biological-measurement alternatives and independent-data feasibility before promoting the existing donor-threshold result to the primary contribution. Continue a small, relevant clinical source-extraction batch alongside these tests. Use a written contribution gate, select at most one primary extension, then freeze and execute its evaluation; preserve failed leads and original artifacts.

**Tech stack:** Windows PowerShell, bundled Python 3.12/NumPy, standard-library unittest and TSV/JSON, existing threshold routines, primary publications/public repositories; existing local plotting dependencies if needed.

**Spec:** `docs/superpowers/plans/2026-10-02-strengthen-biology-project.md`, plus the user's subsequent authorization to investigate metadata, independent comparable data and novelty, reopen the path decision if needed, and delegate execution to Luna. This is an advisory execution plan, not student submission prose.

## Global constraints

- ROOT is `C:/Users/dylan.bagley/Desktop/outputs`; all relative paths below resolve there.
- No external messages, registration, publication, account access, paid services, lab access or biological engineering. Read public primary sources and public data only.
- Preserve original labels, predictions, results, historical reports and existing frozen specifications. Never create amended clinical labels without the original eligibility/review gate passing. No reviewer is identified.
- Previously examined data remain development/exploratory data. Parent review of agent work is technical review, not independent clinical adjudication or new biological validation.
- School fair is late January; state fair is mid-March; weekly availability is unspecified. No timetable depending on fixed student hours or guaranteed external access.
- Dylan writes all submitted paper, abstract, poster and research-plan prose. Supply evidence, analysis, figures, outline and explanation checks only.
- No new agent beyond the one Luna executor. The later request authorizes this single-agent delegation despite earlier no-agent instructions.
- New work lives in `strengthening_execution_Oct2/next_stage/`. Existing execution results are read-only. Only append compact status/pointer/assistance records after backing them up and retaining their history.
- No Git repository was found previously; verify cheaply, do not invent commits/worktrees or introduce repository setup.
- Do not repeat permanently closed calibration searches, controlled-access Da Via retrieval, retired EGFR flip, reconstructed joint distributions as measured data, or general literature hunts.
- RNA thresholds are measurement thresholds, not clinical harm/efficacy, protein density or killing thresholds. Donors, not cells/folds, are biological independent units. Zero UMI is not proved biological absence.

## Review focus

1. Missing metadata or confounded preparation cannot become an adjusted biological finding: retain missing/aliased covariates and refuse unsupported causal attribution.
2. External donor IDs, cohorts and source accession overlap must be checked: a republication/reanalysis of GSE193531 is not replication.
3. Normalization and threshold transport need valid units: do not silently transfer a numeric raw-count threshold between assays or normalize away genuine biology and call it corrected truth.
4. Protocol locking does not erase prior exposure; evaluation donors and resampling units must be explicit, and any fitting/transforms exclude evaluation subjects.
5. A runnable analysis or attractive figure is not a contribution gate pass: state the exact prior-art gap, disconfirmation criterion, useful result and remaining limitations.

## Read-before-change references and clean starting facts

Read this plan, the original plan's global constraints and Tasks 2–7/decision gates, `PROJECT_STATUS_Oct2.md`, `strengthening_execution_Oct2/STATUS.md`, `feasibility/path_decision.md`, `delivery/research_contribution.md`, `delivery/reproduction.md`, `analysis/analysis_protocol.md`, `feasibility/prior_art.tsv`, `evidence/source_gaps.tsv`, and EV-Capture-CURRENT's Known risks/What NOT to revisit. Check applicable AGENTS files; no broad scan through dependency/package folders.

Inspect reusable code: `analysis/run_b_pilot.py`, `analysis/thresholds.py`, `analysis/run_extension.py`, `analysis/tests/`, `evidence/build_batch.py`, `startup.py` (all under strengthening_execution_Oct2). Do not run legacy generators that overwrite preserved outputs.

Starting facts, to verify against files:
- Original 30-target outcome panel: 19 historical favorable/11 unfavorable; fixed top_share+M2_flag AUC 0.913875598086. These labels are unadjudicated, and higher score is favorable.
- GSE193531 / Boiarsky 2022 DOI 10.1038/s41467-022-33944-z: local marker counts and exact cell metadata alignment; 9 healthy NBM donors/9,329 cells, 8 MM individuals with neoplastic cells. CD138-selected plasma cells, not whole marrow; three genes TNFRSF17, GPRC5D, SLAMF7. Source expression-derived identity labels.
- Local data: `ev_capture_data/healthy_bone_marrow_boiarsky/extracted_markers.json` and the same directory's GSE193531 cell-level metadata file. Find exact filename through the existing loader, not assumptions.
- Existing primary illustration at 10% normal-cell fraction: pooled BCMA threshold reaches 24.6% in one healthy donor and 24.9% pooled MM coverage; observed-donor envelope reduces pooled MM coverage to 9.0%. Withheld donors sometimes violate envelope constraints. Raw/normalized scale changes results.
- Source labels, library depth, immunoglobulin composition, preparation/batch, age/sex may explain differences; this stage must test what is actually observable. Missing information stays missing.
- Local myeloma ADT inspected previously is an assay/population/gating mismatch, not comparable independent replication. Da Via is controlled access; do not revisit absent concrete new access evidence.
- Prior work: Boiarsky heterogeneity/same data; Smith 2019 and Verkleij 2021 normal/malignant antigen expression; classical ROC and distribution geometry. Generic heterogeneity, normal-expression selection and a new ROC theorem are not novelty claims.
- Clinical fresh prior batch: CD19 pooled ELIANA/ENSIGN, CD33 final ALFA-0701 and CAIX cohorts; unresolved attribution/negative adequacy. No reviewer contacted/identified.

## Output map

Create only populated files under `strengthening_execution_Oct2/next_stage/`:
- `STATUS.md`, `decision_log.tsv`, `assistance_log.md`, `input_hashes.json`, `output_hashes.json`, `before/`: compact state and recoverable preservation.
- `evidence/metadata_sources.tsv`, `sample_metadata.tsv`, `metadata_limits.md`: sample ID join and inspected provenance.
- `evidence/external_candidates.tsv`, `retrieval_log.tsv`, `prior_art_delta.tsv`, `clinical_sources.tsv`, `clinical_observations.tsv`, `clinical_overlaps.tsv`, `clinical_gap_delta.tsv`: actual searches, access, data and event evidence.
- `data/`: new publicly retrieved inputs only, with URL/access date/checksum/size/source overlap recorded; never overwrite source data.
- `analysis/diagnostic_protocol.md`, `diagnostic_protocol.sha256`, `metadata_diagnostics.py`, `tests/test_metadata_diagnostics.py`, `results/`: bounded metadata/processing tests.
- `feasibility/contribution_gate.md`, `fallback_cards.md` if needed, `path_decision_v2.md`: specific pass/conditional/fail decisions with inspected evidence.
- `analysis/selected_protocol.md`, `.sha256`, `run_selected.py`, `tests/test_selected.py`: only if one selected branch passes its execution gate.
- `delivery/claim_evidence_matrix.tsv`, `research_contribution.md`, `advisory_outline.md`, `student_checks.md`, `reproduction.md`, `figures/`: usable final evidence and technical handoff.

## Task 1 — Safe startup and bounded diagnostic specification

- [ ] Read the references and applicable instructions. Verify root/runtime and existing loader interfaces. Do not reread all historical assessments.
- [ ] Verify every hash in the existing `strengthening_execution_Oct2/input_hashes.json` without modifying it. Snapshot its protected inputs plus current execution source/results/protocol/report files and root/status records into the new hash manifest; exclude dependencies/caches and any previously created next_stage directory. Back up status/pointer records before later authorized updates.
- [ ] Record authorization, deadlines, no reviewer, and the distinction between a completed narrow supporting result and an unestablished substantive extension in new STATUS.
- [ ] Before new diagnostic results, write/hash diagnostic_protocol.md. Keep the same three markers, all eligible donors, 10% primary and 5/20/50% secondary budgets; no result-driven donor exclusion or threshold retuning. Primary diagnostic compares donor-tail variation on original scale against available depth/preparation/composition information. Diagnostic findings are exploratory; do not declare a corrected biological effect.
- [ ] Define descriptive cell and donor diagnostics: per-donor cell count, marker detection fraction, n_counts/n_genes/frac_mito distribution, marker raw/normalized tails, associations between donor-level depth and marker summaries. If composition or preparation is unavailable, say so. With 17 individuals do not fit a many-covariate causal model. NBM/MM are separate strata; do not let disease group confounding manufacture a depth effect.
- [ ] Define a depth sensitivity only when source units support it. If target UMI count c and total captured UMI N are verified integer counts with 0<=c<=N, use hypergeometric thinning to a common captured depth m (m<=N), with NBM/MM selection rules fixed by depth before marker outcomes. Otherwise omit thinning and document why. Do not treat observed gene counts or filtered depth as a validated full library. No RNA molecule calibration or universal dropout correction.

Acceptance: inputs preserved, protocol versioned, exact diagnostic assumptions/selection/seed recorded before numerical execution. Seed 20261003 for new stochastic diagnostics, independent from historical seed. No new clinical metrics.

## Task 2 — Retrieve source metadata and run the specified alternatives test

- [ ] Inspect Boiarsky paper Methods, supplements, GEO sample records and author repository `https://github.com/getzlab/Boiarsky-etal-2022`. One initial pass plus one targeted gap-fill; at most three materially distinct retrieval routes for a named blocked fact. Log attempted routes and exact inspected locations.
- [ ] Populate sample_metadata.tsv with original sample_ID, accession if established, disease state, prep/fresh-frozen, batch/library, age/sex and composition availability. Every nonmissing value needs primary provenance; NR is unknown, not a common batch. Check one-to-one versus one-to-many relationships and incompatible metadata before joining.
- [ ] Write new diagnostic code importing read-only existing threshold routines if appropriate. Use csv readers for quoted TSV, exact sample-ID joins, explicit missing values, finite numeric validation, source units and donors separated.
- [ ] Before relying on it, write meaningful unittest cases: duplicate/conflicting metadata rejects merge; absent covariate stays unknown; cell IDs reorder only by explicit keyed join; unknown sample IDs fail; patient grouping does not split a donor; thinning count never exceeds m and known cases (c=0, c=N, m=N) are exact; invalid c/N rejected; no thinning if units not verified. If thinning omitted, test its refusal condition instead. Do not mirror implementation with tautological assertions.
- [ ] Run unit tests then new standalone diagnostics end to end into new results. Generate a small depth/donor figure if it resolves a question. Do not overwrite original B results or protocol.
- [ ] Write metadata_limits.md separating supported diagnostic associations, unidentifiable technical/biological effects, and any premise weakened by the test. No association between depth and counts proves causality. Retain all negative/inconclusive results.

Acceptance: actual retrieval + reproducible measured diagnostics, not just another proposed analysis. Report whether technical explanations can be tested or separated; lack of metadata is a gate limitation.

## Task 3 — Test independent comparable data and exact prior-art gap

Run alongside the limited clinical task below in the same single executor, alternating tasks as needed; no extra agents.

- [ ] Question: does healthy-reference choice materially alter antigen comparisons, and does a donor-aware threshold policy improve transport to independent comparable healthy/malignant plasma-cell measurements? This is a candidate question, not a novelty assertion.
- [ ] Search public primary study/data records for up to six distinct candidate datasets, ordered by healthy and malignant plasma-cell donor access/comparability. One initial search pass and one gap-fill for the best concrete lead. A study paper alone is not usable data. Stop on controlled access or unavailable critical healthy reference; preserve candidate-specific reasons.
- [ ] For each candidate record accession, authors/cohort provenance, assay, normal comparator, malignant label method, donor IDs/n counts if known, marker availability, prep/normalization, actual files inspected, data size/access route, GSE193531 overlap and eligibility. Unknown overlap blocks independent-replication claims. Normal marrow without enough identified plasma cells is not automatically comparable.
- [ ] Inspect metadata and a real expression/marker slice for eligible candidates. Prefer processed marker tables or small matrices. Do not download whole raw sequencing projects; total new downloads at most 500 MB, each at most 250 MB. If a viable dataset exceeds this or is tool-blocked, report the precise file/dependency rather than invent data or spend indefinitely on retrieval. Network read failures may be attempted through approved available tools; no account creation or contact.
- [ ] Read closest primary studies for cross-donor/cross-cohort reference robustness, transportability, and the selected antigen comparison. Extend prior_art.tsv with exact methods/results inspected, shared versus proposed work, and the specific unanswered comparison. At most eight additional primary papers across this stage; metadata/source Methods access is not a separate general prior-art hunt. A failed search does not prove nobody did it.
- [ ] Write an explicit external-data evaluation feasibility verdict. Separate independent samples, shared-source reanalysis, independent assay of a different question and clinical validation. Different scales may allow a within-study standardized policy comparison, not transport of the exact RNA numeric threshold. Require a prespecified valid cross-study estimand before any final evaluation.

Acceptance: named actual file-level evidence or named missing file; no abstract-only promises. A plausible dataset remains conditional until inspected.

## Task 4 — Clinical repair batch with fixed scope

- [ ] Use existing endpoint_v03 and source gaps. Prioritize TNFRSF17/ide-cel source PMID 33626253 and GPRC5D/talquetamab PMID 36507686; third source only if it resolves a specific material ambiguity in these two. This connects to the antigen extension and keeps audit work bounded.
- [ ] Inspect primary Methods/safety tables and exact denominator definitions for depletion, hypogammaglobulinemia, infection, skin/oral/nail events, severe grades, cessation/DLT, follow-up and supportive care. Split agents, cohorts, event families/windows. Do not replace historical cohort with a convenient responder subgroup or later independent agent.
- [ ] Populate new clinical_sources/observations/overlaps using original plan's schemas and attribution/state vocabulary. Record exact source location and access scope. NR/not_extracted/NA never become zero. Distinguish intended-target attribution from severity, payload/combination effects and author speculation.
- [ ] Record overlaps with existing evidence and prior project exposure. Assistant extraction is provisional, never independent review. Produce a concise clinical_gap_delta; do not amend labels, fit context models or compute new clinical AUC.
- [ ] Check denominator/count units, invalid patient count>denominator, source joins, blank initial review fields if supplying addendum material. Retain no-reviewer dependency, but do not make it a blocker to computational work.

Acceptance: actual limited primary-source gap repair or exact access failures, no general clinical sweep.

## Task 5 — Contribution gate and one bounded path decision

- [ ] For B state eight original gate answers: exact new claim, closest work, opened usable data, independent units/confounders, disconfirmation, fair baseline/evaluation, credible January deliverable, retained assets/new work.
- [ ] Pass only if there is a specific empirical contribution beyond known heterogeneity/ROC, a supported estimand, enough matched metadata or a comparable independent evaluation to address dominant ambiguity, and a feasible evaluation with honest precision. Do not pass solely because computation runs or a preliminary difference is large. Biological attribution may remain limited even when independent measurement transport is informative; state exactly which claim passes.
- [ ] Conditional if a particular feasible analysis requires a named public file or unresolved metadata; fail as primary if only supporting sensitivity remains or prior work exhausts the proposed finding. Preserve current narrow B result regardless.
- [ ] If B cannot supply a substantive primary extension, reopen once and screen at most TWO adjacent fallbacks, without fully implementing both: (F1) assay detection/depth-limited identifiability of one named antigen selectivity conclusion, using verified counting units/actual references and a comparison that prior work does not already answer; (F2) uncertainty-aware clinical evidence benchmark with a real comparative applicability result, using the existing evidence pool, supportable objective observations and target/cohort independence. A catalog alone, generic dropout discussion, or unadjudicated toxicity predictor does not pass. Independent review may leave F2 conditional.
- [ ] For each fallback use one small actual-data pilot and closest primary comparison within the remaining source budget (or record a specific solvable gap justifying one extension). No third fallback, unrelated pivot, or award promises.
- [ ] Select exactly one primary branch, conditional named dependency, or no demonstrated strong extension. Record path_decision_v2 before full implementation. Never silently redefine substantial improvement as a new limitation report.

Acceptance: decision is falsifiable and evidence-linked; distinguishes validity gains, useful supporting extension and demonstrated substantive contribution. If no branch passes, complete the evidence/figures/handoff and report failure candidly; do not invent a selected full analysis.

## Task 6 — Freeze and execute the selected extension if its gate passes

- [ ] Write selected_protocol.md with exact input hashes/sample IDs, primary estimand, baseline, uncertainty assumptions, normalization, evaluation role/exposure, direction, primary 10% and secondary budgets where applicable, missingness/exclusion rules, donor-level resampling, seed and disconfirmation/precision criterion. Hash before final results; no claim that developmental analysis becomes pristine.
- [ ] For the B transport branch compare pooled, equal-donor and observed-donor-envelope policies, a simple mean-expression contrast, and measured coverage/normal fraction violations on identical individuals. Document whether external evaluation uses locked transferable thresholds or a distinct within-study policy comparison. Source-specific fitting stays inside development; no outcome-informed cohort/gene choice.
- [ ] For F1 identify a specific measurement-supported uncertainty set and whether exact finite bounds or simulation sensitivity is warranted. No grid as continuous proof, arbitrary clinical weights or invented true expression. For F2 use target-grouped comparisons and objective/context-specific observations; endpoint attribution validity remains separate, and no amended clinical classification without gate.
- [ ] Write independent toy expected-value tests for boundaries/ties/all-none response, missing data, donor weighting, leakage/overlap and reproducibility. Reuse established numerical code only after reviewing its interfaces. Include a failure test demonstrating an evaluation donor cannot affect fitted threshold/transform.
- [ ] Run tests and affected new standalone entry points to completion. Record donor counts, cohort counts, uncertainty/precision and failed constraints. Never bootstrap cells or dependent folds as new patients. Large uncertainty may make result inconclusive; non-significance is not equivalence.
- [ ] Produce result tables and a small figure set that answer the selected scientific question; each headline traces to a result and protocol. Preserve failures and protocol departures in decision_log.

Acceptance: complete evaluated branch with its actual supported/negative/inconclusive claim. Parent later verifies numerical/source traces; execution does not stop at writing a future work list.

## Task 7 — Deliver, preserve and hand back for parent review

- [ ] Write claim_evidence_matrix with prior knowledge/new result, independent units, assumptions, counterevidence, allowed wording and prohibited overstatement. research_contribution answers what the project adds, strongest result, best competing explanation, remaining uncertainty and whether substantial improvement is actually demonstrated.
- [ ] Supply an advisory outline and student explanation checks, not submission prose. If the contribution gate fails, say so in title/first paragraph and provide the exact dependency/stop decision, not reassurance.
- [ ] Write reproduction with commands, versions, hashes, source access, random seeds and every unresolved failure. Serialize TSV correctly and JSON without nonstandard NaN/Infinity (all-none threshold represented explicitly).
- [ ] Verify original protected hashes and baseline-stage manifest remain unchanged. Backed-up root/status records are allowed append-only pointer/history changes. Compute new output manifest excluding dependency/cache and itself.
- [ ] Append compact old STATUS and ROOT PROJECT_STATUS with the new checkpoint; maintain existing content/history. Update assistance log honestly distinguishing Luna implementation, parent technical review pending and student source verification pending.
- [ ] Return to parent: completed tasks, passing/conditional/failed gate, exact strongest supported finding, precise remaining dependencies, files to review, commands to rerun and preservation count. Do not claim final parent verification yet.

## Verification commands and executor interface

Use bundled runtime (confirm first):
```powershell
& 'C:/Users/dylan.bagley/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' -m unittest discover -s strengthening_execution_Oct2/next_stage/analysis/tests -v
& 'C:/Users/dylan.bagley/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' strengthening_execution_Oct2/next_stage/analysis/metadata_diagnostics.py
```
Run run_selected.py only when implemented following a passing gate; document its exact command. Original 14 tests may be run read-only once as preservation checks; do not run historical output-generating programs. Existing plotting dependencies may need approved sandbox access; do not reinstall everything or hide access failures.

Reusable read-only interfaces to inspect before using:
- `run_b_pilot.load(scale='normalized')` returns donor groups containing marker arrays (confirm signature/units in source).
- `thresholds.select_threshold(tumor_arrays, normal_arrays, budget, policy)` supports pooled/equal_donor/donor_envelope and common empirical thresholds.
- `thresholds.fraction_above(values, threshold)` and startup TSV/hash helpers; inspect actual signatures and import paths.
No changes to these legacy interfaces are authorized by this plan; add new code in next_stage instead.

## Parent final review (parent executes after Luna returns)

- Read contribution gate, decision, source/metadata tables, actual implemented code and changed reports; check original-user constraints and all scope departures.
- Independently recompute the main headline from source/result rows; check denominator, donor independence, held-out exclusion, normalization, access/overlap and exact supporting primary source. Check figures visually where generated.
- Rerun new meaningful tests and every affected new standalone entry point, inspect locked protocol/source hashes, verify original protected inputs and allowed status changes.
- If there are defects, send concrete corrections to the same Luna agent; parent does not implement. Review/rerun affected checks after corrections.
- Only then report completed work, negative findings, substantive-strength verdict and exact dependencies. Technical review cannot manufacture clinical independence or scientific novelty.

## Plan self-review

This plan covers the four authorized next steps, clinical repair alongside feasibility, original two-track gates, bounded searches/two fallbacks, artifact preservation, independence distinctions and student authorship. It deliberately makes selected implementation conditional on actual contribution/data gates. Missing access or novelty is an outcome, not a placeholder. No repeated permission gate is added; user already instructed execution with one isolated Luna at Max. Parent writes no execution code or scientific outputs before handoff.

## Targeted access exception — issued during execution, before external evaluation

Luna identified GSE305967 as a public lead with a 423.8 MB RAW.tar listing and healthy/MM subjects. The 250 MB per-file limit above was a planning workload default, not a user restriction or a scientific gate. For this single lead, allow one targeted access pass and a per-file exception up to the existing 500 MB total budget if the archive contains usable expression tables rather than raw sequencing reads. Prefer individual GSM supplement tables or author-supplied processed data first. Record actual file contents, count units, normal-cell identification, patient/cohort overlap, transfer estimand and access failures; do not treat the archive listing as inspected data. The exact donor-policy question must be compared to prior work, not rejected merely because generic target discovery already exists. This exception does not authorize unbounded downloads/searches, contacts, controlled access or evaluation before protocol locking. Parent changes only this plan; Luna performs retrieval and implementation.
