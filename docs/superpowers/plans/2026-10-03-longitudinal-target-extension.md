# Longitudinal Target Complementarity Implementation Plan

> **For agentic workers:** One gpt-6-luna executor at Max reasoning, with an isolated context fork. Execute the tasks and their decision gates in order. Do not spawn further agents. The parent performs the final review and verification; all analysis, code, and report implementation belongs to Luna.

**Goal:** Test whether a feasible, original biological extension can establish when expression-based target complementarity persists or fails across verified serial myeloma samples, and execute it if its gates pass.

**Architecture:** Resolve the research gap and paired-data feasibility before substantial computation. Run one small, locked pilot, decide, then execute only a passing extension. Preserve the completed retrospective project and negative extension results; do not equate a completed feasibility assessment with achieving the scientific goal.

**Tech Stack:** Windows/PowerShell; bundled Python 3.12, NumPy and openpyxl; standard-library TSV/JSON; existing plotting dependencies only. Scientific packages may require approved read access. No installation is assumed.

**Spec:** Dylan's October 3 request authorizes the preceding five-step proposal: research-gap check, paired-data feasibility, minimal pilot, hard decision, and conditional full execution. Broader disease/measurement changes must be discussed before a full alternative project is built.

## Global constraints and starting facts

- ROOT is `C:/Users/dylan.bagley/Desktop/outputs`; all relative paths below resolve there.
- Create all execution artifacts in `longitudinal_extension_Oct3/` (STAGE). Do not alter any previous analyses, manifests, labels, predictions, clinical observations, plots, protocols, or status files. This stage's STATUS is the resume entry point. No Git/worktree or commit is assumed.
- No external messages, registration, paid access, publication, credentials, or laboratory work. Public primary-source browsing and bounded public processed-data retrieval are authorized. Do not bypass access controls.
- Dylan authors submitted paper, abstract, poster, and research-plan prose. Provide evidence, figures, outlines, explanation checks, and factual assistance logs only.
- School fair is late January; state fair mid-March; weekly hours are unspecified. Target one working week for the first gate, not an instruction to sleep or schedule an automation.
- The original clinical finding remains a selected-panel retrospective association, unadjudicated. No clinical reviewer exists; independent attribution is not solved by assistant review. Do not repeat clinical extraction or change classifications.
- The completed RNA threshold extension failed its strong-contribution gate. Stricter constraints mechanically changed detection; GPRC5D was invariant. Depth and clinical fallback screens also failed as primary extensions. Do not recycle these as new discoveries.
- GSE305971 workbook: 416,184 annotations, 311 Hospital.Code values, 366 code/biopsy keys, 45 codes with repeated biopsy sequence values. This does NOT verify patient identity, chronology, treatment exposure, or useful pairs.
- Previous first-biopsy expression parse: 456 files, 62,811 exact joins, 2,697 unmatched IDs; W/WMC aliases unresolved. All expression already inspected is development data, including new endpoints computed from it.
- Study/participant/biopsy/cell are different units. Repeated patients do not create independent subjects. RNA zeros are not proven antigen absence; joint RNA detection is not protein coexpression or treatment response.

## Review focus

1. False pairing: require documented participant identity, specimen identity and temporal order; biopsy sequence alone is insufficient for treatment timing.
2. False novelty: published pair selection, antigen escape and dual targeting already exist; name an unanswered contrast, not a generic claim that heterogeneity matters.
3. Leakage: select samples from metadata, disclose all prior expression exposure, and keep patient groups together for any validation split.
4. False joint measurement: only markers measured in the same cells may define joint detection. Reject duplicate cell keys, wrong batch joins, fractional molecule counts and unverified aliases.
5. False clinical meaning: distinguish measured transcript detection, inferred protein availability, observed response, and causal treatment resistance.

## Task 0 — Preserve and orient

- [ ] Read these current references selectively: `strengthening_execution_Oct2/next_stage/feasibility/contribution_gate.md`, `path_decision_v2.md`, `fallback_cards.md`, `prior_art_delta.tsv`, `external_data_assessment.md`; `strengthening_execution_Oct2/next_stage/delivery/reproduction.md`; `strengthening_execution_Oct2/next_stage/analysis/results/external_annotation_identity_audit.json` and `external_annotation_by_hospital_biopsy.tsv`.
- [ ] Read the known risks and what-not-to-revisit sections of `EV-Capture-CURRENT.md`; inspect applicable local AGENTS files without scanning dependency trees. User scope takes precedence over older restrictions on extension work.
- [ ] Create `STAGE/STATUS.md`, `logs/decision_log.tsv`, `logs/assistance_log.tsv`, `preservation/input_hashes.json`, and a preservation verifier. Record hashes of the original 63 protected inputs, previous 65-file snapshot, previous 155-output inventory, manifests themselves, and directly used prior files. Preserve the prior inventory as a read-only source; do not regenerate it.
- [ ] Verify the existing output inventory before relying on it. All previous files remain unchanged throughout this stage. Use a durable ledger instead of commits in this nonrepository.

## Task 1 — Establish the precise gap

Create `feasibility/question.md`, `prior_art.tsv`, `access_log.tsv`, and `gap_gate.md`.

- [ ] Check these primary starting sources and their actual accessible methods/results: Zada et al. DOI `10.1038/s41588-026-02725-5` (GSE305971); antigen escape DOI `10.1038/s41591-023-02491-5` (GSE226336); combinatorial single-cell targeting DOI `10.1038/s41587-023-01686-y`; myeloma multiomic target discovery, Cancer Research 83:1214 (2023); `https://pmc.ncbi.nlm.nih.gov/articles/PMC13004696/` (2026 GPRC5D escape). Add at most three directly relevant primary papers if needed, one targeted gap-filling pass. Do not claim full-method access from an abstract.
- [ ] For each source record URL/DOI, accessed sections, exact overlap, cohort, assay, temporal design, published comparator, candidate remaining gap, and access limitations. Author code/source tables are preferred to proxies; no exact-method superiority claim without faithful replication.
- [ ] Candidate question: do target pairs selected from initial samples retain measured complementarity in later samples, including named malignant states, and does a prespecified selection method anticipate losses better than an appropriate simple baseline? Target family initially fixed to TNFRSF17, GPRC5D, SLAMF7; no unrestricted search for striking pairs.
- [ ] Separate two possible contributions: an evidenced biological pattern in an underexamined serial setting; or a validated improvement in prediction/selection. An expected OR coverage increase or a stricter cutoff is neither.
- [ ] Gate G1: continue only for a precise potentially unanswered contrast that the available sources do not clearly already answer. If access prevents establishing it, label novelty unresolved and constrain execution to feasibility. A failure is recorded, not repaired by renaming the same calculation.

## Task 2 — Verify paired-data feasibility

Create `feasibility/datasets.tsv`, `sample_inventory.tsv`, `pair_inventory.tsv`, `data_gate.md`; `data/download_manifest.tsv` and metadata copies only as needed.

- [ ] First inspect the existing GSE305971 workbook `strengthening_execution_Oct2/next_stage/data/GSE305971_cells_annotation_for_geo_ver_20250306.xlsx`, its GEO sample records, source supplements and author metadata. Existing workbook hash is `2fc3f9e19ab82097f6a50689de7bf4083aeb81aba442f263ea81afe8b150b94d`.
- [ ] Inspect at most two additional directly relevant public series, starting with the expression subseries of GSE226336 and the 2026 escape study's deposited expression data. Genome/CNV measurements are not RNA or protein. Check individual processed-file availability before assuming that a large archive is needed.
- [ ] Sample schema: dataset, source_url/location, patient_key, identity_evidence, specimen_key, temporal_order, order_evidence, sampling_time, treatment_class, exposure_evidence, disease_state, tissue, assay, source_malignant_label, eligible_cell_count, file_accession, batch_key, prior_expression_exposure, usable_for_serial_question, unresolved_reason. Missing values are NR, not inferred dates or treatment.
- [ ] Pair schema: pair_id, dataset, patient_key, initial_specimen, later_specimen, identity_verified, order_verified, intervening_treatment_verified, assay_comparable, cell_label_supported, file_available, new_expression_exposure, gate_status, exclusion_reason. A repeated sample code is only a candidate pair.
- [ ] Use metadata alone for candidate inclusion. Choose first and last comparable eligible specimens per verified participant; preserve other timepoints for declared sensitivity. Do not select pairs by marker change.
- [ ] New downloads capped at 500 MiB total compressed processed data, no single file above 250 MiB unless individually justified within the cap. Reuse verified existing files read-only. No FASTQ processing or 33-GB archive. Missing access/metadata is a scientific dependency, not a reason to keep downloading indefinitely.
- [ ] Gate G2: at least three verified paired participants with accessible comparable measurements can support a minimal descriptive pilot, not a generalized or competitive-strength claim. Treatment-specific wording also requires treatment timing. Fewer pairs, no valid joint measurement, or no identity/order evidence closes this candidate. Full execution requires a precision/validation justification stronger than this workload screen.

## Task 3 — Lock and execute one minimal pilot, only if G1/G2 pass

Create `analysis/pilot_protocol.md` + SHA sidecar before new expression parsing; `analysis/load_paired_counts.py`, `run_pilot.py`, `tests/test_paired_analysis.py`; `analysis/results/` tables; one pilot figure and its value table.

- [ ] Specify the exact question, assays, sample keys, primary endpoint, markers/pairs, cell labels, normalization, baseline, missingness handling, depth sensitivity, analysis units, exposure history and disconfirmation criteria. Hash downloaded inputs and record lock timing. Never call already viewed data pristine validation.
- [ ] Endpoint starts as per-participant change in same-cell transcript complementarity: for ordered markers A/B, fraction of eligible cells detected for B among those not detected for A. Also report joint four-category counts, OR fraction and double-undetected fraction. A zero conditional denominator is unavailable, not zero. These are transcript measurements only.
- [ ] Primary transcript detection uses strictly positive observed integer molecules where the assay supports this. If the actual source scale cannot support it, record a pre-result amendment or stop. Never transfer an absolute cutoff across assays. Report raw joint counts and source normalization; analyze depths without equating downsampling to causal correction.
- [ ] Compare against initial-sample mean-expression selection and any faithful accessible published comparator. Define selection/ties before inspecting later samples. If selection methods cannot meaningfully differ with the chosen candidates, do not invent a poor baseline. No method superiority claim from descriptive profiles alone.
- [ ] Use all nonempty eligible specimens; flag small denominators and quantify one-cell resolution. No after-result cell-count exclusion. If state labels are unavailable, do not infer states using the endpoint marker genes; keep the pilot specimen-level and explicitly lower its scope.
- [ ] Tests use independent small examples: same marginals but different joint distributions; absent conditional denominator; ordered timepoints; patient train/test disjointness; reordered exact cell joins; duplicate/wrong-batch rejection; invalid count rejection; missing observations never becoming zeros; deterministic pair selection/tie handling. Example joint counts `(A+B+=2,A+B-=1,A-B+=3,A-B-=4)` imply OR=6/10 and B detection among A-undetected=3/7.
- [ ] Run tests, then the standalone pilot on all selected specimens. Keep participant-level values visible; no cell-as-person p-values. Depth/missingness explanations must be tested where identifiable, otherwise listed as alternatives. Compare the observation to the source paper's reported findings.

## Task 4 — Make the hard contribution decision

Create `feasibility/path_decision.md` addressing: exact incremental claim; closest prior work; opened data; independent units/confounders; disconfirmation; fair baseline; deadline-feasible deliverable; retained assets/new requirements.

- [ ] Proceed to full analysis only if G1/G2 pass, the pilot provides a testable result beyond trivial pair coverage/known antigen escape, and a credible biological interpretation plus evaluation is feasible. Require enough participants and measurement resolution for the named conclusion; justify from actual pilot dispersion/precision, not an arbitrary universal count.
- [ ] If the result is exposure-influenced, it can support a disclosed biological analysis but not untouched predictive validation. A new holdout must be genuinely unexamined for the relevant outcome.
- [ ] If the candidate fails, do not rerun the old threshold audit. Prepare at most two adjacent disease/measurement alternatives as concise gap/data cards, with inspected primary metadata and specific questions. Do not build a full alternative project without discussion. Return the failed gate and recommended scope change candidly.

## Task 5 — Conditional full execution and validation

- [ ] Only for a passing decision, create `analysis/full_protocol.md` + hash and a self-contained `run_full.py` reusing tested loaders. Lock participant/specimen inclusion, endpoint, baseline, uncertainty, splits, interpretation and stopping conditions before additional results.
- [ ] Include eligible paired participants, preserving within-patient dependence. For prediction, fit all selection/tuning on training patients and test on disjoint unexposed patients; otherwise describe paired longitudinal measurements. Do not label one patient's later specimen an independent patient replication.
- [ ] Run one appropriate validation check using held-out patients or verified independent cohort if feasible. If not feasible, explicitly limit claims and reassess whether the strong-extension gate survives. Cross-assay comparisons require compatible estimands, not shared numeric cutoffs.
- [ ] Report per-patient effect, coverage, change/rank stability, paired baseline differences, missingness and depth sensitivity. Any bootstrap resamples patients, preserves paired specimens and records undefined draws. Small-sample intervals are descriptive and do not establish causal therapy effects.
- [ ] If the full data contradict the pilot, reopen the decision once, preserving both results. No new success endpoint chosen from the failed run.

## Task 6 — Solidify and deliver

- [ ] Create `delivery/claim_evidence_matrix.tsv`, `research_contribution.md`, `advisory_outline.md`, `student_checks.md`, `reproduction.md`, `output_hashes.json`; numeric-axis SVG + 300-dpi PNG figures with source-value TSVs if analysis runs. No empty figure placeholders on a failed gate.
- [ ] Create `delivery/remaining_dependencies.md`: clinical reviewer remains unidentified; source-method access and identities relevant to this stage; student explanation checks; future adult/school review. The existing score-withheld clinical packet is retained, not re-extracted or sent.
- [ ] Update only STAGE STATUS/logs, including precise stop/proceed outcome, next two actions, assistance roles, command exits, seeds, hashes, deviations and failed leads. Keep previous stage and root status unchanged.
- [ ] Run meaningful new tests and every implemented standalone entry point; verify all prior protected hashes, new output inventory and strict finite JSON serialization. Preserve temporary failed/provisional scientific outputs with labels; do not recursively delete history.
- [ ] Return a compact executor report linking gate evidence, actual analyses, verification commands, preservation outcome, and any unmet goal. Parent independently verifies counts/numerical results, source scope, plots, hash preservation and standalone behavior. Return corrections to the same Luna agent; parent must not implement them.

## Execution and verification commands

Runtime: `C:/Users/dylan.bagley/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe`.

Run new tests with `-m unittest discover -s longitudinal_extension_Oct3/analysis/tests -v` when that directory exists. Run implemented inventory, pilot and full entry points by their exact saved paths and record them in reproduction.md. No test suite or full generator is fabricated when its branch did not pass. Verification tools must be stage-local and must not regenerate old outputs. Required invariant: every record in the previous output manifest and original protected-input manifest retains its original hash.

## Plan self-review

This plan covers all five authorized steps, persistent status, preservation, student ownership and final parent verification. Gates deliberately leave the scientific answer open. It does not assume that 45 repeated codes are valid pairs, that a new endpoint is novel, that RNA establishes target availability, that three pilot patients justify generalization, or that negative feasibility work satisfies the substantial-extension objective.
