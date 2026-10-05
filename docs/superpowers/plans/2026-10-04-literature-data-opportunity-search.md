# Literature and Data Opportunity Search Implementation Plan

> **For agentic workers:** Execute this research plan task by task. Use `superpowers:executing-plans` if useful, but do not spawn agents. Dylan explicitly approved one isolated gpt-6-luna executor at Max reasoning and parent review between stages; that authorization overrides generic confirmation or additional-agent workflows. Mark completed steps in the new status file, not by rewriting this plan.

**Goal:** Identify several important, potentially original computational-biology questions supported by accessible measured data, then give the parent enough evidence to select the strongest questions for a separate feasibility plan.

**Architecture:** Sequential review-led field mapping, primary-paper screening, concurrent data-access inspection during reading, and a ranked candidate synthesis. This first cycle produces an opportunity map, not a selected extension or a model. Parent verifies the map and writes the next narrowly scoped plan; the same Luna executor handles subsequent authorized execution.

**Tech stack:** Primary publications and official repositories; Markdown, TSV, JSON; Windows PowerShell and existing Python standard library for records and integrity. No new packages are required for this cycle.

**Spec:** Dylan approved the staged workflow on October 4, 2026: detailed plans, isolated Luna Max execution, parent scientific review between stages, and continuation while evidence supports useful further work. He seeks substantial biological improvement, preferably medical, rather than another cleanup pass. This plan incorporates the scope and stopping conditions explicitly.

## Global constraints and exact paths

- ROOT: `C:/Users/dylan.bagley/Desktop/outputs`. All paths below are relative to ROOT unless absolute. Current local date October 4, 2026, America/Phoenix. Record true UTC timestamps for execution logs; do not invent elapsed times.
- PLAN: `docs/superpowers/plans/2026-10-04-literature-data-opportunity-search.md`.
- STAGE: `opportunity_search_Oct4/cycle1`. Create this directory and write all execution outputs/downloads/cache there. Never edit previous project outputs, source data, label files, predictions, reports, status files, or old dependencies.
- Existing Python: `C:/Users/dylan.bagley/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe`. Set `PYTHONDONTWRITEBYTECODE=1` before imports from old trees. Standard-library scripts suffice; no installation, raw sequencing, training, or large archives.
- One executor only. Do not spawn agents or contact people. No emails, external messages, publication, registration, paid services, controlled-access requests, or lab work. Public-source reading and bounded downloads are authorized.
- Dylan writes submitted paper, abstract, poster and formal research-plan prose. Produce evidence, analysis plans, advisory outlines and feedback only. Do not ghostwrite submission text.
- School fair late January; state fair mid-March. Year not explicitly confirmed. No fixed weekly hours and no clinical reviewer. Do not assume either dependency will appear.
- Preserve existing mathematics and computational skills as possible assets; do not force the old confinement score or binding equations onto a new biological question. Disease, modality and method may vary within biology. Prefer a connection to existing work when scientifically natural; distinguish an extension from a broader pivot.
- Use reviews for orientation, primary papers for scientific claims, and official data repositories for access/schema. Do not declare a future-direction sentence proof of an open problem. Locate later work, including preprints with clear status, before declaring a gap.
- Do not equate assistant review with independent clinical adjudication or prospective validation. Keep measured responses, calculated responses, simulations, inferred states, fitted parameters, technical replicates and independent biological units separate.
- No test suite or document generator for this stage. A small integrity validator is appropriate; no pipeline/framework building before selection.

## Review focus

1. Future directions may already be solved: each shortlisted gap needs a later-work search and named closest competing study.
2. Public access statements may hide summary-only data: distinguish identifier/link found, metadata read, file accessed, and actual rows inspected.
3. Many cells or wells may represent few independent specimens: identify biological grouping and holdout limitations; never turn subsamples into patients.
4. Attractive method ideas may lack biological value: require a concrete decision, mechanism or biological question and explain the consequence of resolving it.
5. Novelty/data claims may be inherited or misattributed: verify titles, DOI/URLs, exact figures/table sections and which quantities were actually observed.

## Task 1 — Orientation, preservation, and boundaries

**Create:** `STAGE/STATUS.md`, `STAGE/preservation_before.json`.

- [ ] Read PLAN completely. Check applicable AGENTS.md in ROOT and ancestors, without traversing dependency directories. Follow direct user authorization if historical instructions conflict.
- [ ] Read `functional_extension_Oct3/STATUS.md` and `functional_extension_Oct3/decision.md`, `longitudinal_extension_Oct3/delivery/path_decision.md`, and `PROJECT_STATUS_Oct2.md`. For `EV-Capture-CURRENT.md`, read only `Known risks` and `What NOT to revisit`; older headline claims are superseded by the newer status.
- [ ] Check whether ROOT has `00_START_HERE_ChatGPT.md` and `00_MANIFEST.md`. These were not found previously; record absence if still absent and use the present references. Do not broad-scan dependencies or invent missing files.
- [ ] Hash the plan, the above reference files, and the three prior output manifests before research: `functional_extension_Oct3/output_hashes.json`, `longitudinal_extension_Oct3/delivery/output_hashes.json`, `strengthening_execution_Oct2/next_stage/delivery/output_hashes.json`. Record path, bytes and SHA256 in preservation_before.json. Capture real start UTC and current stage in STATUS.md.
- [ ] Summarize exclusions in STATUS.md: no renewed myeloma healthy-cutoff/coexpression/RNA longitudinal variants; no reanalysis of the three venetoclax participants; no Umetani pre-exposure non-growing-fraction versus survival question; no repeat finite-panel HIV pair-selection proposal; no closed EV source/calibration searches. These are specific closed proposals, not proof that entire fields are closed.

Existing result: the original clinical score is a selected-panel retrospective association with uncertain labels, not a general clinical predictor. RNA/protein mismatch, uncalibrated binding parameters, selected outcome panels and limited independent units constrain its extensions. Previous narrow screens did not find a substantial extension. Do not spend this cycle auditing that conclusion again.

## Task 2 — Review-led map across distinct biological problems

**Create:** `STAGE/search_log.tsv`, `STAGE/papers.tsv`.

- [ ] Find 4–6 substantial reviews, preferably published 2022–2026, spanning at least three distinct biological problem families. Initial search directions: functional cancer drug response/resistance; antigen recognition and therapeutic selectivity; dynamics/lineage-associated escape or response; signaling/apoptosis or another medically relevant noninfectious response process. These are search categories, not nominated projects. At least one family must extend beyond antigen-expression ranking.
- [ ] Prefer reviews whose limitations and references lead to measurable, unresolved questions and public processed data. Do not choose papers merely because they use heterogeneity/distributions/single-cell terminology. At most two reviews from any one narrow family.
- [ ] Read relevant substantive sections, not abstracts alone. Record actual reading level. Extract important unresolved claims, unexplained/contradictory results, measurement limitations, static-versus-dynamic mismatches and authors' future directions. Record why each may matter biologically and what evidence would resolve it.
- [ ] Follow useful references and citation chains into primary work. Search later publications explicitly for each emerging gap. Reviews orient the search but cannot establish novelty.

Search-log TSV columns (tab-separated, UTF-8):
`search_id, utc, purpose, query_or_url, results_used, next_action_or_stop_reason`.
Paper TSV columns:
`paper_id, type, title, year, doi_or_accession, url, problem_family, reading_level, sections_read, measured_inputs, measured_outputs, independent_units, key_limitation_or_gap, candidate_links, access_limit`.
Reading_level allowed: `metadata_only`, `abstract`, `targeted_sections`, `full_main_text`, `main_plus_relevant_supplement`. Use `not_reported` or `not_inspected` rather than guesses. Verify bibliographic identities against the primary page.
For all TSVs, use actual tabs as field separators, no tabs/newlines inside cells, and `|` between multiple IDs. Use IDs R01/P01 for review/primary papers, D01 for resources, C01 for candidates, S01 for searches. Reference checks split multi-ID cells on `|`; missing values use `not_reported`, `not_inspected`, or `none` as appropriate.

## Task 3 — Focused primary-paper screen and data inspection

**Create:** `STAGE/data_cards.tsv`, optional small source files under `STAGE/sources/` and processed samples under `STAGE/data/`.

- [ ] Screen approximately 15–25 distinct primary papers across the mapped problems, including useful older foundational work. A paper counts as screened only if its question, measured endpoint, and contribution were examined. Closely read approximately 8–12 of the most relevant, including competing work; quotas are ranges, not permission to pad or fabricate evidence. Document if access or relevance yields fewer.
- [ ] During reading, follow official data availability and repository records for at least 4–6 promising resources. For the strongest resources, inspect the smallest legitimate processed file or example rows that reveal the actual schema. Dataset-first opportunities may emerge here; formulate biological questions from what is measured, not solely the authors' suggestions.
- [ ] Record identifiers, file names/size, current access, modality, baseline versus post-treatment timing, outcome type, sample IDs, grouping, joins, replicates, censoring/detection limits, interventions and controls. Identify whether a separate study/cohort could evaluate the proposed finding. Distinguish independent-cohort, same-study grouped holdout, temporal holdout and exploratory reuse.
- [ ] Verify measured endpoint status directly. Drug sensitivity inferred from RNA is not measured drug response; predicted combination values are not combo assays; inferred trajectories are not observed lineages; technical wells are not independent specimens. Mechanistic/causal claims require intervention/design support; computational predictions may be useful without such claims.
- [ ] Stop a blocked source after one legitimate alternate route (official mirror, repository or author manuscript). No repeated download troubleshooting. Summaries can support a descriptive question but do not substitute for linked individual records if that is what the question requires.

Data-card TSV columns:
`data_id, paper_ids, accession, official_url, access_status, file_examined, bytes_if_known, biological_system, input_measurements, outcome_measurements, timing, independent_unit, independent_unit_count, grouping_keys, join_keys, censoring_or_limits, controls, evaluation_resource, measured_vs_derived, candidate_links, unresolved_access_or_schema`.
Access_status allowed: `identified`, `metadata_verified`, `file_accessed`, `rows_inspected`, `restricted`, `blocked`, `summary_only`. Give uncertainty separately; never upgrade access from a citation alone.

## Task 4 — Form 6–8 concrete candidate questions, then rank

**Create:** `STAGE/candidates.tsv`, `STAGE/opportunity_map.md`.

- [ ] Form 6–8 genuine candidate questions across at least three problem families from the inspected evidence. Do not produce eight minor variations of one question. If fewer survive, state why; retain rejected leads briefly to prevent repeat work.
- [ ] For each candidate, write one precise sentence describing inputs, outcome and biological question. Explain what knowledge or decision could change if the answer were positive OR negative. Specify connection to the old project as `natural_extension`, `related_pivot`, or `broader_biology_pivot`.
- [ ] Name the closest competing papers and the exact increment beyond them. Search for that increment being already addressed. A different dataset, higher AUC, bigger table, generic tail statistic, resampling exercise or previously known combination optimization is not sufficient by itself.
- [ ] Identify one primary endpoint, fair published/simple comparator when applicable, the independent unit, accessible development/evaluation data, and a minimal decisive experiment. Label novelty `overlap`, `unresolved`, `plausible_gap`; data `inadequate`, `partial`, `schema_verified`; biology `minor`, `potentially_meaningful`, `well_motivated`. A plausible gap is not validated novelty.
- [ ] Rank strongest 2–3 using written reasoning, not a weighted score that creates false precision. Hard gates: clear important question; plausible distinct contribution; practical measured data; valid group/timing design; a feasible experiment that could falsify the hoped-for benefit. Do not rank an inaccessible-data candidate first merely for excitement.
- [ ] Explain disconfirming evidence, risks and what each subsequent feasibility test would resolve. Do not let a missing preferred second cohort automatically imply all research is impossible: identify the defensible claim scope, and distinguish descriptive/mechanistic/predictive requirements. Conversely do not call internal holdout independent external validation.

Candidate TSV columns:
`candidate_id, family, question, biological_value, relationship_to_existing, paper_ids, closest_prior_art, exact_increment, later_work_search_ids, novelty_status, data_ids, data_status, independent_unit, primary_endpoint, comparator, evaluation_design, minimal_decisive_test, success_gate, revise_or_stop_gate, main_risks, rank_or_rejection_reason`.
Opportunity-map Markdown: brief overall recommendation; family map; 6–8 concise candidate cards; comparative shortlist; exact next feasibility steps; negative findings/access limitations. Aim <=2,500 words excluding TSVs. No submission prose or invented results.

## Task 5 — Integrity and parent handoff

**Create:** `STAGE/data_manifest.json`, `STAGE/output_hashes.json`, `STAGE/validate_artifacts.py`; update only `STAGE/STATUS.md`.

- [ ] Save a compact JSON manifest of actual downloaded files with source URL/accession, bytes, SHA256, content type and measured/derived scope. Empty downloads are legitimate; do not claim rows inspected if none were.
- [ ] A short standard-library validator must: parse finite JSON; check consistent TSV widths, unique IDs and reference links; verify downloaded and stage-file hashes; compare preservation_before.json; verify all files in the three previous output manifests using their correct base directories. Include all stage outputs except output_hashes.json itself and explicit cache exclusions. Validate counts against listed files and actual directory contents.
Prior manifest bases are respectively `functional_extension_Oct3`, `longitudinal_extension_Oct3`, and `strengthening_execution_Oct2/next_stage`; rows use `path`, `bytes`, `sha256`, and manifest `file_count`. New stage rows use the same keys and paths relative to STAGE. Hash validation must not write to any old directory.
- [ ] Run validator using the existing Python runtime with bytecode disabled. Run `longitudinal_extension_Oct3/preservation/verify_preservation.py` read-only to verify historical original protections and documented status suffixes. Do not restore or modify historical files to make them match an older snapshot.
- [ ] Check titles/DOIs again, particularly shortlist sources. Confirm main text versus supplement reading accurately recorded. Confirm no empirical findings/pilot results fabricated. Refresh output hashes only after edits complete.
- [ ] Final STATUS <=900 words: true start/end UTC; completed scope and counts; shortlist; what remains unresolved; preservation/validation; exact next actions. Return concise parent handoff naming strongest candidates, disconfirming facts and actual source-reading/access limits.

## Bounds and stage decision

- This first cycle is deliberately broader than the prior two-candidate screen. Use up to 24 focused search queries plus 6 targeted gap-fill/citation searches, and citation-link traversal as needed. Log searches; no claim of exhaustive field coverage.
- Target completion within roughly 90 minutes active work; at 60 minutes send parent progress with coverage and remaining key gaps. Hard checkpoint at 120 minutes: finish records with candid limits and return for review. Record real elapsed time and distinguish access waits; do not claim a clock was exceeded before it was.
- Downloads: <=100 MiB total, <=25 MiB per file, processed/tabular/PDF sources only. Check size first when available; unknown-size large archive is not authorized. No raw sequencing/images, giant matrices or model weights. If a promising processed file exceeds the cap, record exact URL/size and proceed with metadata, not silent rejection of the whole question.
- No minimal pilot in cycle1: parent must first review the opportunity map and write a separate plan. No new full project or infrastructure before that gate.
- If nothing earns a feasibility test, report it without manufacturing a winner. Identify whether the bottleneck is novelty, data, biological value or validation, and which genuinely different search would address it. Do not automatically repeat the same search. Parent decides a subsequent bounded plan within the already authorized workflow.

## Parent responsibilities after return

Parent reads the map, candidate/data/source tables and validates key primary-source claims for the top 2–3, including later work and exact schemas. Parent reruns integrity/preservation checks. Corrections are assigned to the same Luna, not implemented by parent. Once corrected, parent writes an exact feasibility plan specifying shortlisted questions, access tests, fixed contrasts, valid comparators/units, smallest decisive experiments and proceed/revise/stop gates. Same isolated Luna executes it. Parent then reviews experiments and writes either a development plan for an evidence-supported question or a clearly different search plan. Continue toward a meaningful extension without guaranteeing a positive result or accumulating analysis on unsupported candidates.
