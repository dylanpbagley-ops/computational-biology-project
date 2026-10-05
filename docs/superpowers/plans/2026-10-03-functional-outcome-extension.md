# Functional-outcome extension selection and pilot plan

> **For agentic workers:** One isolated gpt-6-luna executor at Max reasoning. Execute sequentially with superpowers:executing-plans if useful. Do not spawn further agents. The parent reviews sources, calculations, preservation and the recommendation after execution. All implementation and research execution belong to Luna.

**Goal:** Select, or decisively reject, a substantial medically relevant extension based on real experimental outcomes rather than another myeloma transcript comparison.

**Architecture:** Examine two distinct questions using existing source leads. Verify the closest published work and actual input/outcome matching first. Run at most one minimal pilot per eligible candidate; select one route only if a specific incremental contribution survives fair comparisons. Preserve the existing project as a completed component. Do not equate negative feasibility work with achieving the extension goal.

**Tech stack:** Windows/PowerShell, bundled Python 3.12, NumPy, standard library, existing openpyxl and optional project-local scientific dependencies. Reuse packages read-only; any essential installation is stage-local, official, pinned, hashed and without changing the bundled runtime.

**Spec:** Dylan's current request authorizes this plan, one isolated Luna/Max executor, the bounded investigation and minimal pilots, and parent final review. The preceding accepted strategy is a biological question with measured experimental outcomes, at most two distinct candidates, and the smallest decisive test before infrastructure. No full alternative project is authorized in this selection stage.

## Global constraints and files

- ROOT: `C:/Users/dylan.bagley/Desktop/outputs`. STAGE: `functional_extension_Oct3/`. All implementation, downloads, logs, outputs and status go inside STAGE. Only this plan is outside STAGE.
- Preserve original labels, scores, predictions, models, previous analyses, manifests and all old status files. No overwrite, relabeling, deleting history, restoring old snapshots, commits or worktree creation are needed.
- No external messages, publication, registration, credentials, paid/controlled access requests, lab activity or scheduled jobs. Public primary-source browsing and bounded processed-data retrieval are authorized. Do not bypass access controls.
- Dylan writes submitted paper, abstract, poster and research-plan prose. Deliver scientific evidence, analytic figures if warranted, a factual recommendation and student explanation points, not submission prose.
- School fair: late January; state fair: mid-March; available weekly hours unspecified. Judge a computational deliverable against those deadlines without assuming lab or reviewer access.
- No clinical reviewer is identified. Clinical attribution remains unresolved; do not spend this stage re-extracting clinical evidence.
- Preserve the previous longitudinal outcome: three verified venetoclax RNA pairs, 44,981 barcodes, strong-contribution gate failed. No further myeloma RNA, healthy-reference threshold, donor marker-panel or Cell Painting acquisition pilots in this pass.
- Prior direction exploration was inconclusive. Antibiotic and HIV files are leads, not selected projects. Old expression/data already inspected remain exposed development material; disclose exposure and never call it pristine validation.
- Generic heterogeneity, nonlinear averaging, pair coverage, better fit on training data, a stronger AUC or classical probability bounds alone are not sufficient contributions. A failed simple baseline is not evidence of superiority to published methods.
- Use `PYTHONDONTWRITEBYTECODE=1` or `sys.dont_write_bytecode=True` before old-module imports; put Matplotlib configuration/cache in STAGE. Old directories must receive no writes.

Keep the package small: `STATUS.md`, `selection.tsv`, `sources.tsv`, `data_manifest.json`, `protocol.md`/`protocol.lock.json` if any pilot passes, small loader/pilot and meaningful tests, `results/`, `decision.md`, `reproduction.md`, and `output_hashes.json`. Add only artifacts that affect a scientific decision or reproduction. Do not build a documentation generator, general framework or separate script for every prose change. Update STATUS immediately after each gate.

## Review focus

1. Real outcomes versus simulated/model-derived outputs: prove provenance and measurement meaning; a file named trial analysis is not automatically clinical data.
2. Prediction timing and circularity: predictor information must precede the outcome; do not predict an outcome with its own post-treatment measurements.
3. Independent units: identify culture/run, donor, virus lineage, isolate, plate and technical replicate; cells, concentrations and timepoints are not interchangeable independent experiments.
4. Increment over the closest source: check the source's own model, supplements and available code before calling an application new; access limitations mean unresolved novelty, not a pass.
5. Fair evaluation: shared information, training-only fitting, grouped splits, assay scale, censoring, detection limits and comparator fidelity must support the exact claim.

## Task 0 — Orient and preserve (small setup, no new audit campaign)

- [ ] Read `longitudinal_extension_Oct3/STATUS.md` and `delivery/path_decision.md`; `direction_landscape_Oct1/search_protocol.md`, `Stronger_Project_Options_Oct1.md`, and the relevant antibiotic/HIV inventories/provenance. Read only the known-risks and what-not-to-revisit sections of `EV-Capture-CURRENT.md`. Inspect applicable AGENTS files in root/STAGE ancestry.
- [ ] The old `00_START_HERE_ChatGPT.md` and `00_MANIFEST.md` were not found in the root or a filename-only search of available project folders. Do not invent them or block on them; current status and source files are the operative references. Record this orientation limitation once.
- [ ] Create STAGE and a compact STATUS with the two question slots, passed/failed gates, current blocker, actual runtime versions and next action. No old STATUS append.
- [ ] Preserve by recording current hashes of this plan, directly read old sources and the following inventories. Verify all entries in `strengthening_execution_Oct2/next_stage/delivery/output_hashes.json` (155 old outputs) and `longitudinal_extension_Oct3/delivery/output_hashes.json` (201 stage files, manifest SHA `0a5a069983a6f32d6f3ab09e50c94a20738f310ea3c6a2a5fbcabe4723a6dcb7`). Run the old longitudinal preservation verifier read-only to cover 63 original protected inputs, the old snapshot and reconciled pre-existing status suffixes. Do not restore two legitimate historical suffix updates. Save a compact before/after record within STAGE.
- [ ] Existing model core may be read from `ev_capture_model/src/binding_model.py` and `aim10_car_t_kill_probability_model.py`. Transfer the idea of distribution-dependent biological response, not uncalibrated physical parameters or the assertion that unrelated systems share one response equation. Do not reopen settled EV calibration/theorem work.

## Task 1 — Two specific candidate questions and source challenges

Investigate exactly these two substantive question slots. They are hypotheses to test, not claims of established feasibility or novelty. A source may answer a question already; record failure promptly. Do not replace failed endpoints silently or add a third candidate.

### A: Antibiotic survival from population history

Candidate question: **Can measured pre-exposure growth/lag heterogeneity, or a explicitly defined early killing observation, improve prediction of independently measured later survival/regrowth beyond population means and the closest published survival model?** Distinguish these two possible input types in metadata; choose only one precise input/outcome contrast before outcome calculations. Do not mix an early-curve forecasting task with pretreatment prediction. If the source already answers the chosen contrast, reject it rather than rebrand it.

Starting evidence, all local and read-only:
- `direction_landscape_Oct1/antibiotic_aging_article.txt` and its HTML/provenance: "Observation of universal ageing dynamics in antibiotic persistence", DOI `10.1038/s41586-021-04114-w`.
- `aging_source_fig1.xlsx` through `aging_source_fig4.xlsx`, their provenance and `numeric_inventory.json`. These are source workbooks, not guaranteed matched individual-history/outcome data.
- `aging_metadata.json`: author code release `NQBLab/RCCN` at `https://github.com/NQBLab/RCCN`; the listed archive is about 261 MB. Inspect a small exact code/data route rather than retrieving that archive.
- `bacterial_dryad.txt`, `bacterial_file_metadata.json`: "Observation of persister cell histories reveals diverse modes of survival in antibiotic persistence", data DOI `10.5061/dryad.s1rn8pkb1`; listed zip `DataForDryad_Uploaded.zip`, 1,098,934 bytes. Validate current primary methods, conditions, independent runs and actual released columns before using it.

Required distinctions: death versus non-growth versus censored persistence; culture/run versus cell; pre-exposure versus post-exposure history; actual CFU/viability/regrowth versus simulated trajectories; prior paper's ageing/hazard law versus a new predictive contrast. Relevant comparator must include the source's suitable conventional/age-dependent/biphasic model if the question overlaps it. Mean-only comparison is insufficient.

### B: HIV antibody susceptibility and finite sampling

Candidate question: **Does limited sampling of baseline viral susceptibility change a fixed antibody-treatment decision in a way that predicts withheld experimentally measured neutralization or clinical rebound better than a strong published susceptibility/combination model?** Choose one actual measured endpoint and one fixed decision problem from inspected data before analysis. Measured isolate-level neutralization can support an assay-level claim; it cannot validate patient rebound. If only simulation outputs or published-model predictions are available, fail the real-outcome gate for that endpoint.

Starting evidence, read-only:
- `direction_landscape_Oct1/hiv_combinations_article.txt` and HTML provenance, primary PLOS Pathogens DOI `10.1371/journal.ppat.1005520`; supplement DOCX and extracted text/table inventory.
- `hiv_prediction_repo.json`, `hiv_tree.json`, `hiv_escape_sample.csv` and provenance.
- `hiv_trial_combo.h5` and provenance: source `https://raw.githubusercontent.com/StatPhysBio/HIVTreatmentOptimization/main/H5Output/trialanalysiscombo.h5`, 267,024 bytes. The H5Output location makes provenance inspection essential: do not treat this as measured patient data without proof. Public author code: `https://github.com/StatPhysBio/HIVTreatmentOptimization`.

Required distinctions: measured neutralization versus predicted inhibition; donor versus isolate/lineage/pseudovirus versus assay replicate; baseline sample versus post-treatment rebound sequencing; susceptibility and concentration censoring; correlated resistance and published Bliss-Hill/additive models. Combination optimization and heterogeneity are already established. Require a specific unresolved finite-sampling decision/evaluation contrast, not another optimization exercise.

### Shared bounded search and eligibility gates

- [ ] Inspect old lead files first. Make `selection.tsv` with candidate_id, precise_question, biological_decision, measured_outcome, input_timing, independent_unit, direct_old_asset_transfer, closest_prior_art, exact_increment, dataset, access, status and reason.
- [ ] Primary-source browsing only for scientific support. For each candidate: source paper/data/code above plus at most three additional directly relevant primary papers/datasets; up to six focused search queries and one targeted gap-fill. At most one alternate legitimate access route per blocked source. Stop repeated CAPTCHA/network retry loops. Do not build a new broad landscape catalog.
- [ ] Put source title/DOI/URL, accessed sections, overlap, access limitation and exact data/code evidence in a single `sources.tsv`. Check current primary source and repository metadata, not inherited summaries alone. Distinguish source text, listed files, opened files and reproduced results.
- [ ] Inspect only the minimal processed files needed to prove schema, units, measurement scale and linkage. New data downloads capped at 100 MiB total and 50 MiB per file. No raw images/FASTQ, controlled datasets, large repository archives or training a large model. Small existing files may be reused read-only with exposure disclosure.
- [ ] Novelty gate: a specific potentially unanswered contrast survives the closest paper's methods/results, not just its abstract. If missing supplements/code prevent evaluating overlap, mark unresolved and do not recommend the candidate as established original.
- [ ] Data gate: actual accessible records support the exact input/outcome contrast, meaningful scale, provenance, temporal order, joins and independent-group evaluation. Quantify the eligible group count and missingness. No universal arbitrary sample-count rule: state what this count permits and whether grouping/censoring leave enough variation for the named claim. Fewer than three independent groups cannot earn a robust multi-group result. A tiny cohort can be a data inspection only.
- [ ] If either gate fails, record the reason and move to the other candidate. Stop A's data/model work when it fails; do not spend an hour reproducing a source that already answers it. If both fail, proceed directly to Task 3. No fabricated protocol/tests/figures for a pilot that did not run.
- [ ] Budget: aim to decide these two eligibility gates in about 45 minutes of active execution, with a hard checkpoint at 60 minutes. At that checkpoint unresolved access/matching stays unresolved and no larger retrieval begins. This is a research-work budget, not an instruction to wait/schedule; use actual timestamps and do not backdate. No hard calendar promise is made. Send the parent a brief checkpoint immediately when both gate decisions exist.

## Task 2 — Small decisive pilot, only for eligible candidate(s)

- [ ] For each passing candidate (maximum two), write a short protocol section and lock its SHA/timestamp before comparing target outcomes. State existing exposure, inclusion, group IDs, input timing, endpoint, biological question, competing predictions, fitting, fair baselines, exact metric, evaluation groups, practical effect criterion and failure rule. No invented untouched holdout. A published/reused dataset with known results may support disclosed retrospective evaluation, not prospective validation.
- [ ] Resolve the observed schema in the protocol: input_id, outcome_id, group_id, condition, time, measurement unit, censoring and normalization. Reject duplicate/wrong-condition joins and ambiguous groups; missing values stay missing. For an early-curve task explicitly bound information available at forecast time, and keep an entire experiment grouped. For HIV do not let related variants from the same person straddle a claimed person-level split.
- [ ] Implement one focused loader and one standalone pilot, ideally less than 300 lines combined before necessary tests. Avoid adding classes/frameworks/general pipelines. Standard Python/NumPy and existing packages suffice for small source tables; if an essential package is missing, use an official pinned stage-local wheel with no dependency changes, or report the runtime dependency. Missing optional software is not a biological infeasibility finding.
- [ ] Use meaningful tiny test fixtures for: same mean/different distribution; exact group split with no overlap; predictor-time restriction; duplicate/wrong-condition join rejection; missing/censored values not treated as zeros; and an independently computed metric on known values. Add system-specific checks for the actual source scale. Do not write tests that merely mirror code or create artificial populations as biological validation.
- [ ] Keep preprocessing/model fitting inside training groups. Compare population-mean response, a simple biologically appropriate baseline, and the closest feasible published model with comparable inputs/training data. Source calibration transferred from unrelated experiments must not be labeled independent outcome validation. If the strongest comparator cannot be faithfully used, restrict the claim and mark superiority unresolved.
- [ ] Run the smallest grouped evaluation that can distinguish the proposed increment. Report every held-out group, denominators, censored/undefined outcomes, calibration where pertinent, and the actual difference from each baseline. Use a participant/culture/experiment as the resampling unit if uncertainty is estimated. Do not pool cells/timepoints as independent replicates. Do not tune endpoints, treatment panels, splits or thresholds after seeing a failure.
- [ ] A new method need not be built if a well-defined biological model comparison tests the question. Reproducing an established result, a trivial nonlinear-versus-mean gap, a gain only in simulated data, or a gain over a deliberately weak baseline does not pass.
- [ ] Candidate pilot implementation budget: about 45 minutes active work per eligible candidate, one focused repair if a routine bug appears. If larger infrastructure or inaccessible source replication is necessary, record the dependency and stop rather than silently expanding this stage. Protocol lock is an integrity record, not preregistration.
- [ ] Save small numeric value tables and source hashes. One analytic figure only if it clarifies an actual result; numeric axes, explicit units, descriptive scope, SVG plus 300-dpi PNG with existing plotting dependencies. No poster/abstract prose or decorative placeholders.

## Task 3 — Choose a defensible path or recommend a broader pivot

- [ ] Create `decision.md`, at most about two pages. Compare A and B against biological importance, exact prior-art gap, actual data, independent evaluation, fair baseline performance, transfer from the original project, remaining dependencies and January feasibility. No arbitrary aggregate score that hides a fatal gate.
- [ ] Use precise outcome labels: `reject`, `unresolved dependency`, `pilot only`, or `recommend scoped extension`. State what observed result would disconfirm the recommendation. Distinguish a promising pilot from achieving the substantial extension.
- [ ] Recommend one candidate only if its exact increment, measured outcome and credible evaluation survive. Provide a concrete next implementation block for that single route: dataset, hypothesis, endpoint, baseline, held-out units, next validation and stop rule. Do not execute a full new project in this stage.
- [ ] If neither candidate earns recommendation, do not propose a third nearby variation or recycle the closed myeloma/CITE-seq/imaging pilots. Recommend an explicit broader pivot decision, naming the binding obstacle (question already answered, inaccessible matched outcomes, unsuitable independent units or no useful pilot increment). Preserve the original project as completed work. State honestly that no demonstrated extension is selected.
- [ ] Update STATUS with result, current files, failed leads, next two actions and assistance roles. Write short `reproduction.md` with exact commands, seeds, dependency/access needs and observed exits. Record original/memory-derived versus freshly verified facts. Student explanation checks can be a few bullets in decision.md.
- [ ] Verify preservation against captured old hashes, all original protected files and both old inventories. Validate finite JSON and TSV widths. Create `output_hashes.json` for scientifically relevant new files and downloaded data/dependency provenance; explicitly state any excluded volatile caches and exclude the manifest itself. Do not hash hundreds of irrelevant installed/cache files into a scientific deliverable.
- [ ] Return to parent with both gate decisions, actual pilot results if any, favored path or candid pivot verdict, exact reproduction commands, source links, preservation outcome and unachieved goal. No further agents. Parent will independently check sources, raw counts/units, calculations, standalone behavior, figures and preservation; all corrections return to this same Luna.

## Runtime and verification interface

Python: `C:/Users/dylan.bagley/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe`.

Useful existing h5py, if needed: `longitudinal_extension_Oct3/dependencies` (h5py 3.16.0, NumPy 2.3.5). Existing plotting/scientific packages: `strengthening_execution_Oct2/dependencies`; filesystem sandbox may require approved read access. Never interpret a partial module import caused by AccessDenied as a scientific gate failure. Do not modify old dependency trees.

For a implemented pilot, save runnable `functional_extension_Oct3/pilot.py` and tests in `functional_extension_Oct3/tests`. Required commands from ROOT:

```powershell
$env:PYTHONDONTWRITEBYTECODE = '1'
$taskPython = 'C:/Users/dylan.bagley/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'
& $taskPython -m unittest discover -s functional_extension_Oct3/tests -v
& $taskPython functional_extension_Oct3/pilot.py
& $taskPython longitudinal_extension_Oct3/preservation/verify_preservation.py
```

If a pilot branch does not pass its gates, do not create or claim a test suite/pilot command. Supply the actual saved verification commands instead. Each executed script's output scope must be only STAGE; preservation verifier reads old files without changing them.

## Plan self-review

The plan covers the accepted broader strategy and exact delegation method. It fixes two distinct functional questions while leaving assay/clinical endpoints conditional on actual inspected data. It has search/download/active-work caps, fair-comparator and independent-unit gates, small pilots before infrastructure, explicit rejection/dependency outcomes and a final path decision. It preserves previous projects, closes failed directions, avoids assumed novelty and simulated-outcome substitution, and leaves submitted prose to Dylan. No extra permission question is needed: the user explicitly authorized planning and this bounded delegation. No full alternative project is built before the path review.
