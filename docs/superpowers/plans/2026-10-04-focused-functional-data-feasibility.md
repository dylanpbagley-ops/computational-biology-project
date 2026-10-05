# Focused Functional Data Feasibility Implementation Plan

> **For agentic workers:** Execute sequentially using the existing isolated gpt-6-luna Max executor. No additional agents, external contact, or repeated user confirmation. This stage implements the parent's reviewed next steps, not a full new project. Research gates override implementation enthusiasm.

**Goal:** Determine whether either of two biologically meaningful questions has a defensible gap and accessible matched measured data sufficient for a small decisive pilot.

**Architecture:** First test early-response/later-regrowth data and prior work; then test directly measured antigen-presentation/recognition data and prior work. Read decisive primary sources deeply, inspect actual small processed records, and return a specific proceed/revise/stop decision. Parent reviews before writing a pilot plan.

**Tech stack:** Primary literature, official repository/public author files; Markdown/TSV/JSON; existing Python standard library and already installed scientific runtime if genuinely needed for schema inspection. No model training, full pipelines, raw datasets or new dependency installs.

**Spec:** Dylan approved iterative detailed plans, isolated Luna Max execution, and parent review between stages on October 4, 2026. Cycle1 mapping is in `opportunity_search_Oct4/cycle1`. Parent reviewed the candidate map, independently checked key source claims and reran its integrity and original preservation checks successfully. Medical biology is preferred; the proposed mechanism need not be evaluated with patient clinical outcomes if the claim is ex vivo.

## Paths and constraints

- ROOT: `C:/Users/dylan.bagley/Desktop/outputs`.
- PLAN: `docs/superpowers/plans/2026-10-04-focused-functional-data-feasibility.md`.
- STAGE: `opportunity_search_Oct4/cycle2`. Write all new scientific artifacts, scripts, downloads/cache here. Leave cycle1 and all previous files read-only.
- One additional writable handoff file: `opportunity_search_Oct4/STATUS.md`, the new compact current-search pointer. Create/update it with current cycle, plan, gate decisions and next step; do not rewrite old project status files.
- Python: `C:/Users/dylan.bagley/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe`; set `PYTHONDONTWRITEBYTECODE=1`. Existing old dependencies remain read-only; prefer stdlib zip/XML for small XLSX/DOCX inspection and CSV/JSON files. Matplotlib not needed.
- No additional agents, people contacted, publication, registration, lab work, paid services or controlled-access requests. Do not request credentials. A publicly accessible author deposit is allowed, asking the author for files is not.
- Dylan writes submitted prose. Evidence records and analysis/protocol outlines only. Preserve old labels/predictions and every previous stage manifest.
- School fair late January, state mid-March; year not confirmed. Weekly hours unknown; no independent reviewer identified. Do not assume an experiment or clinical collaboration will become available.
- No reopenings of old EV calibration/source searches, RNA cutoff/myeloma variants, three-participant venetoclax analysis, Umetani non-growing-fraction contrast, or HIV finite-panel pair-selection proposal.
- Do not force the old confinement/binding model onto these questions or tune an endpoint to obtain an effect. Biological value and novelty are separate from access and software correctness.

## Review focus and known findings

Cycle1 mapped 6 reviews, screened21 primary/preprint/resource papers and8 questions, and made10 resource cards. Only one primary had a full main-text read; most records are targeted excerpts or abstract-level. HOPE was the only actual processed-row preview, with no local downloads. This is a useful lead map, not an adequate deep review or a feasibility pass.

Parent verified these distinctions:
- TESLA primary DOI `10.1016/j.cell.2020.09.015`, PMC7652061, describes measured patient-matched T-cell recognition but its label “presented” is a model filter using affinity/abundance/stability. It is not direct mass-spectrometric observation of natural tumor presentation. Binding stability is predicted in this paper's methods. Do not interpret this as a measured two-stage failure dataset.
- PROSPERO, PMID37171617/DOI `10.1007/s00018-023-04772-1`, is2023, and already compares acute single-cell markers with later viability. That alone is neither washout-regrowth evidence nor an open novelty gap. Main datasets are on request in the checked availability statement. Small XLSX/CSV supplements and public author pipeline are possible legitimate next routes.
- Generic organoid AUC/PK combination scoring, generic single-cell resistance signatures and broad CD123-DBP prediction are already covered; do not use them as the extension.

Critical failure modes: early marker positivity is not necessarily irreversible death; late ATP is not reproductive survival; continued-drug exposure is not washout; mass-spectrometry nondetection is not proof of absent presentation; synthetic decoys are not measured biological negatives; exogenously loaded peptide recognition is not proof of natural tumor-cell recognition; many peptides/cells/wells do not create independent patients.

## Task1 — Freeze references and deepen decisive reading

**Create:** `STAGE/STATUS.md`, `STAGE/preservation_before.json`, `STAGE/sources.tsv`, `STAGE/search_log.tsv`.

- [ ] Read PLAN, cycle1/opportunity_map.md, candidates.tsv rows C02/C05, their cited paper/data rows and final cycle1/STATUS.md. Do not reread the full historical project.
- [ ] Hash PLAN, the cycle1 manifest and the named current cycle1 references. Start true UTC clock and create both new status files. Current Phoenix date remains October4; derive actual timestamps, not from old filenames.
- [ ] Read relevant main methods/results and supplements of approximately4–6 decisive primary sources across the two questions. Prefer full readable sources and clearly record what was actually read. Include at least one meaningful competing method per branch. If fewer are readable, say so; do not relabel indexed snippets as full reading.
- [ ] Official journal/institutional author manuscripts, repositories and legitimate open-access representations such as Europe PMC fullTextXML are acceptable reading routes. Use at most one alternate legitimate route per blocked source within this cycle; no CAPTCHA solving, retry loops, or access-control bypass.

Sources TSV columns: `source_id, branch, title, year_and_status, doi_or_accession, url, reading_level, sections_read, measured_variables, derived_variables, exact_gap_or_overlap, data_links, access_limit`.
Search TSV: `search_id, utc, branch, purpose, query_or_url, result_ids, next_action`.
Use exact tabs, no tab/newline inside a cell, `|` between IDs, IDs A01/B01 for sources, D01 resources, S01 searches; missing values explicitly `not_inspected`, `not_reported`, `none`. Verify title/year/DOI on primary pages.

## Task2 — BranchA: early drug response versus durable regrowth

**Initial question:** In cancer models exposed to a specified treatment, what explains discordance between acute cell-resolved response/death measurements and later ability to repopulate after drug removal?

This revises C02's data requirements openly: allow patient-derived non-organoid models or established cancer lines if that preserves meaningful biology. Do not require all ATP, GR, apoptosis and washout measures simultaneously merely to consider a useful question. Claims of superiority to ATP/GR need those comparators measured or legitimately computable; otherwise restrict the question/claim.

- [ ] First inspect the actual PROSPERO methods, data statement and smallest publicly exposed processed supplements, then author repository `https://github.com/augpath/PROSPERO`. Primary article `https://pmc.ncbi.nlm.nih.gov/articles/PMC11071868/`; possible public small supplements are an XLSX roughly14KB and CSV roughly4.6KB. Determine actual fields/units; do not download160MBPDF or24MB DOCX merely because listed. Repository code alone is not a measured dataset.
- [ ] Verify growth-rate inhibition/GR primary work (Hafner2016), fractional-killing/death metric work, and a relevant delayed-response/washout/rechallenge primary study. Search exact endpoint, not generic heterogeneity. Record what each already measures and the exact increment, if any.
- [ ] Identify up to3 distinct public processed-data resources for paired acute response and later reproductive fate/washout/regrowth. At least one may be an experimentally simpler cancer-line system rather than a clinical organoid cohort. Read experimental schedules before attempting files.
- [ ] Inspect the smallest actual rows for the strongest resource: model/donor ID, independent run, treatment, dose, exposure period, early-readout timing/definition, washout schedule, regrowth endpoint/time, controls and repeated wells. Record original keys and count independent models/biological replicates separately from rows/cells. If only continued-drug viability exists, label it and do not rename it reproductive regrowth.
- [ ] Decide whether the initial question is supported. If exact washout pairing is absent but actual data expose another meaningful response mechanism, propose at most ONE narrower revision derived from those actual measurements. Name the mechanism and precise contrast, check closest competing work, and explain how it differs from generic early-versus-late viability prediction already in PROSPERO. Do not invent fate or call technical repeats independent lines.
- [ ] Stop this branch's access work if all3resources lack requisite matched records. Record specific failed gate and whether a different data-backed revision survives. No lab assumptions, raw image/FCS processing, classifier or pilot yet.

Acceptable A gate: a clear biological contrast, practical measured paired inputs/outcomes, independent units and valid timing/control design, fair feasible comparison, plausible unaddressed increment supported by a close competing source. A dataset passing access alone is insufficient.

## Task3 — BranchB: natural presentation versus T-cell recognition

**Initial question:** Which experimentally observed stage limits recognition of naturally presented tumor peptides, and can presentation and downstream T-cell recognition be distinguished in matched data?

- [ ] Read TESLA primary methods/Figs3–4 and relevant supplemental columns from an accessible official source. Institutional author archive `https://authors.library.caltech.edu/records/y57ap-h8p92` is a legitimate supplement lead identified by parent. Synapse `syn21048999` is known, but do not spend time on login or raw WES/RNA. Verify whether processed peptide/patient assay tables are genuinely public. Record natural presentation, measured binding, predicted stability, expression, loaded-pMHC T-cell recognition and tumor-cell recognition as distinct fields.
- [ ] Read later primary comparator BigMHC: `https://www.nature.com/articles/s42256-023-00694-6` (2023). It trains on MS ligands and transfers to immunogenicity. Determine whether its cohorts are matched across stages or separate training sources; no claim of novelty from separating those stages generically.
- [ ] Check at most2 additional legitimate resources for actual matched naturally presented peptide/HLA measurements and downstream recognition, preferably processed immunopeptidomics plus patient-matched assay tables. One targeted later-work lead is epiVIP2026 preprint `https://www.biorxiv.org/content/10.64898/2026.05.05.722287v1.full`; its abstract distinguishes imputed abundance from measured biology. Verify its claims/status if used; do not assume its predictions provide the required ground truth.
- [ ] Inspect actual processed keys: tumor/patient/model ID, peptide sequence, HLA allele, mutation, source specimen/time, direct MS evidence/intensity/limit, T-cell source, assay procedure and measured endpoint. Count independent biological units; verify joins before proposing modeling.
- [ ] Distinguish identified ligands from synthetic random negatives and unobserved MS candidates. Nondetection alone cannot support the label “presentation failed.” A loaded peptide assay can test recognition of that ligand but does not show the patient's natural tumor presents it. Record selection/censoring and whether these limits prevent the initial causal-stage distinction.
- [ ] If the initial question fails, at most ONE narrower revision may be derived from an actually inspected dataset (for example a specific processing/recognition perturbation contrast). It must have a different biological increment than another peptide-ranking score and a closest-primary-work comparison. Do not create a new candidate from predictions alone.

Acceptable B gate: linked measured natural-presentation/recognition or a clearly justified revised perturbation contrast, valid biological units, selection/censoring handling, and a plausible specific biological increment. Model “presented” labels or peptide-random holdout do not pass.

## Task4 — Data cards, gate decision and minimal-pilot specification

**Create:** `STAGE/data_cards.tsv`, `STAGE/gates.tsv`, `STAGE/decision.md`, `STAGE/data_manifest.json`.

Data TSV columns: `data_id, branch, source_ids, official_url, accession_or_file, access_level, file_path, bytes, sha256, schema, matched_measurements, timing, unit_and_counts, grouping_and_join_keys, control_and_censoring, measured_vs_derived, limitations`.
Access_level: `identified`, `metadata_verified`, `file_accessed`, `rows_inspected`, `blocked`, `restricted`, `summary_only`. If browser preview only, record that explicitly without claiming a downloaded hash/full row count.

Gate TSV: `branch, original_question, revised_question_if_any, source_ids, data_ids, biological_value, novelty_gate, data_gate, evaluation_gate, verdict, precise_next_test, stop_reason`.
Gate states `pass`, `unresolved`, `fail`; verdict `pilot_candidate`, `revise_with_data`, `stop_current_question`. No forced winner.

- [ ] For each surviving question, specify the smallest next pilot in decision.md: one fixed contrast/endpoint, actual variables, inclusion/exclusion, independent unit/grouped splits, feasible baseline, falsifying result, expected uncertainty, and why the finding would matter biologically. Do not fit or simulate it yet. No numerical effect-size gate invented after seeing outcomes; justify any precision/sample threshold against its intended claim.
- [ ] Select at most one primary pilot candidate by scientific value and evidence, not merely file availability. A secondary can remain conditional. Parent decides after review and writes exact pilot plan.
- [ ] If both fail, distinguish original-question failures from all-field impossibility. Identify the most actionable actually inspected resource and ONE specific data-backed follow-on search/revision, with competing work and why it addresses this cycle's bottleneck. If no such lead exists, report that rather than manufacture a third nominal candidate. Do not default back to old clinical repair/logistics analysis.
- [ ] Report exact deep-reading and row-inspection counts and access routes; separate missing data from confirmed absence in an inspected file. Update compact statuses.

## Task5 — Preservation and validator

**Create:** `STAGE/validate_artifacts.py`, `STAGE/output_hashes.json`.

- [ ] Use the existing cycle1 validator pattern, simplifying to this stage's tables/IDs. Check exact TSV widths/unique IDs/foreign references; finite JSON; actual new files against complete manifest; downloaded file source/size/hash; preserved reference hashes. Include every scientifically relevant file; exclude only manifest itself and documented cache/bytecode directories. No new framework or test suite.
- [ ] Rerun cycle1 validator read-only (which checks its9files plus365prior-manifest files), and `longitudinal_extension_Oct3/preservation/verify_preservation.py`. No restoring historic status suffixes.
- [ ] Hash new outputs last, run final stage validator and return decision plus precise sources/data and next pilot/follow-on suggestion. Parent independently reviews/reruns. Root current-search STATUS should point to this finalized cycle and intended next plan; it is an authorized live handoff file outside cycle manifests, not a historical project output.

## Bounded execution and parent next step

- Sequential A then B; no extra agent. Target60–90minutes active; send meaningful source/data milestone at30minutes and freeze research by90minutes. Hard120minute checkpoint: finalize candid records and return. Measure real UTC/wall time; active time only if instrumented.
- Up to12focused search queries total (6perbranch), up to4targeted citation/access gap fills total. Source-link reading does not require extra general searches; log it. Maximum3resources A and3resources B including TESLA. No breadth padding.
- Download budget75MiB total,25MiB per file. Small processed CSV/TSV/JSON/XLSX/DOCX and readable papers only; no raw sequencing/MS/images/FCS or whole repositories/archives. Size-check first. Existing public CSV source is not a request to download all supporting material.
- If a measured endpoint is absent in an inspected file, conclude only about that file. If a source is blocked, record access unresolved, not data nonexistent. Do not spend the whole cycle retrying access.
- Parent will read final source/gate/data records, independently verify decisive source/provenance claims and rerun checks. Any errors go back to same Luna. Parent then writes and delegates a minimal pilot plan for an evidence-supported candidate, or a clearly different bounded data-led plan justified by this cycle. No positive result or substantial contribution is guaranteed; keep proceeding only when a concrete next test can change the scientific decision.
