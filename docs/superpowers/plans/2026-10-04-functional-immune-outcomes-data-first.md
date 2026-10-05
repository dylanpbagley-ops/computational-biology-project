# Functional immune-outcome data-first search implementation plan

> Executor: one fresh, isolated gpt-6-luna subagent, Max reasoning. No additional agents. Execute the steps in order; return to the parent for scientific review before any pilot fitting.

**Goal:** Find a substantive, testable immunotherapy extension using publicly inspectable functional perturbation outcomes, without recycling established gene lists or inventing a winning question.

**Architecture:** Inspect measured outcome tables first, then evaluate at most two precisely stated biological questions against closest prior work. Separate access, originality, and evaluability. Parent reviews the candidate and writes a separate locked minimal-pilot plan if a candidate passes.

**Tech stack:** Existing Windows PowerShell/Python, PubMed connector, official publishers/repositories, small processed XLS/XLSX/TSV/CSV tables; no installs or raw sequencing.

**Spec:** Dylan's current authorization to continue the search using isolated Luna Max; original constraints in `docs/superpowers/plans/2026-10-02-strengthen-biology-project.md`; current negative decisions in `opportunity_search_Oct4/STATUS.md` and cycle1–3 decisions. This is an advisory execution plan, not student submission prose.

## User-specific success criteria and time limit

Dylan reiterated: biological research, strong preference for medical relevance, mainly or entirely computational, and no more unfocused three-hour searches. A worthy extension is a clear and important biological question that the original retrospective association did not answer, with a materially distinct contribution, measured supporting evidence, fair comparisons and a finishable analysis for the school/state fairs. More audits, a gene list, another generic predictor, novelty by recent publication date, or a paper outline do not satisfy this. This stage stays in cancer immunotherapy; an outcome-backed efficacy/context extension is related to the existing targeting project, but is not validation of its original clinical safety score. No lab dependency, biological engineering or unrelated complete pivot is allowed. Competitive strength is assessed candidly, never promised.

Budget correction: executor must deliver by60 minutes wall from its captured start, reserving last10 minutes for closeout, so parent review fits within75 minutes TOTAL for this round. The20-minute actual-row checkpoint is mandatory; if no rows are available, return the specific access problem rather than expanding search. Parent may allow ONE explicit, bounded next step within the75-minute total only when it can change the decision. No silent extensions. Parent plans no further broad rounds during this task without a concrete measured-data lead.

## Context, scope and parent field assessment

ROOT is `C:/Users/dylan.bagley/Desktop/outputs`. New STAGE is `opportunity_search_Oct4/cycle4_functional`. Original normal-tissue score is a selected-panel retrospective RNA association, not a calibrated general safety predictor. Preserve the original labels, predictions, locked rules and every prior stage. No outside reviewer exists. School fair is late January and state fair mid-March; weekly time and deadline years were not confirmed. Student writes submitted paper/abstract/poster/research-plan prose.

Parent inspected current primary metadata/abstracts and targeted Results: generic immune-resistance lists, conserved pathways, IFN effects, and general context dependence are already established. Specific evidence:

- Lawson et al., Nature 2020, PMID32968282, PMC9014559, DOI10.1038/s41586-020-2746-2: six mouse tumor models, CTL vs untreated differential-fitness screens; 182 core genes, genetic co-similarity, fitness overlap and cytokine/in-vivo checks are already reported. GEO GSE149933 is an identified official processed-data lead; actual files have NOT been opened by the parent.
- Larson et al., Nature 2022, PMID35418687, DOI10.1038/s41586-022-04585-5: IFNγR loss makes solid tumors resist CAR-T killing; liquid-vs-solid contrast already published. This is prior art, not the proposed discovery.
- Leukemia-intrinsic determinants, Nature Communications 2023, PMID38052854, PMC10698189: IFNγR/JAK/STAT loss has opposite in-vitro/in-vivo effects; NK/Qa-1b mechanism already tested. Parent read first10,000/101,277 body characters, including introduction and initial Results; not whole-paper verification.
- CELLFIE, Nature 2025, PMID40993398, PMC12545207, DOI10.1038/s41586-025-09507-9: multiple effector-cell readouts, prioritized in-vivo screens, combinatorial pairs and RHOG/FAS findings already reported. Do not mix T-cell perturbations with tumor-cell perturbations or relabel proliferation as killing.
- Recent 2026 primary record PMID42826199, "Tumor cell-intrinsic stress states drive sensitivity to CAR T cell therapy in pancreatic cancer": KEAP1/SLC33A1/NRF2 conclusions already in abstract; only metadata/abstract seen, processed access unknown. New date does not establish novelty.
- Dufva et al., "Integrated drug profiling and CRISPR screening identify essential pathways for CAR T-cell cytotoxicity", PMC7098811: matched pharmacological/functional readouts are a possible alternative within the same immune-killing scope. Published drug mechanisms and general immunomodulation are not new. Resolve PMID/DOI from source, do not guess.

Parent inference: a question about the limits of transporting an immune-specific tumor sensitization effect across measured settings may be useful, if a distinct mechanism/contrast remains and real controls allow testing it. This is NOT an established gap. A second allowable question, only if actual drug tables support it, is a specific tumor-benefit versus effector-impairment conflict. Do not force either, invent a hit, or require a positive effect. Do not introduce pathogen, wet-lab, engineering, clinical treatment, or unrelated project work.

Closed prior paths: EV calibration, myeloma healthy cutoffs/RNA coexpression, venetoclax cell-count association, Umetani nongrowing fraction, generic HIV combination optimization, early viability versus regrowth already studied in CTC models, JCI/Science/TESLA presentation-immunogenicity proposals, and cycle3 inflammatory peptide safety. Reopen only with truly new evidence; these are not fallback tasks.

## Authorized changes and preservation

- Read original plan selectively plus `PROJECT_STATUS_Oct2.md`, current `opportunity_search_Oct4/STATUS.md`, and cycle3 corrected decision. No full history reconstruction.
- Write only new STAGE files and the compact live `opportunity_search_Oct4/STATUS.md`. Parent alone writes execution plans. Do not change root PROJECT_STATUS, original analysis files, earlier stage files, memory, or Git state.
- At start record hashes of this plan, cycle1/2/3 output manifests and prior live status; record live status as authorized mutable, not protected. Preserve prior stage manifests and their referenced artifacts. Use existing preservation validators for original outputs; no duplicate audit framework.
- Public sources only. No messages/contact, accounts, controlled access, publication/registration, paid service, installations or raw FASTQ/H5/scRNA/MS. Unsupported browser binary reader is not evidence of restricted access. Follow observed official public download links; at most one legitimate publisher mirror or repository alternate per resource. Never bypass CAPTCHA/authentication.

## Bounds and checkpoints

- Record actual UTC start/end and query timestamps; never reconstruct missing times. Target45 minutes, hard60-minute executor closeout; total round including parent review is capped at75 minutes as stated above. At20 minutes send parent actual resource/row progress; reserve last10 minutes for compact records and validation. If time threatens the cap, freeze what is inspected and return the exact next decisive step; no silent overrun or broad automatic continuation.
- At most8 new publication-discovery queries (including exact-title/DOI gap fills),40 unique title screens,8 targeted primary-paper Methods/Results reads,4 actual dataset resources. Every page count/cursor and exact screened scope must be recorded. No obligatory pile of reviews.
- At most25MiB per processed file and75MiB aggregate. Check official size headers if available; stream with a byte cap if unknown. Do not download an entire repository archive unnecessarily. Inspect actual schemas and rows from at least one resource early; file-not-found is not a scientific novelty failure.
- At most2 candidate cards and one data-backed revision of the strongest. No models, tuned rankings, performance sweeps, figures or manuscript prose during this feasibility stage. Tiny counts, coverage, fixed joins and assay-direction checks are allowed.

## Task1 — Find and open measured outcomes first

- [ ] Create STAGE and concise STATUS with start clock, next action and protected-input hashes.
- [ ] Start with Lawson2020 official supplementary differential-fitness tables (SupplementaryTable7 and untreated-fitness tables as needed). Inspect publisher and GSE149933 records; download smallest pertinent processed tables, not raw reads. Resolve exact URLs and actual bytes. Parent has not verified file availability.
- [ ] Inspect a genuinely complementary CAR-T tumor-side study next: Larson2022 or leukemia2023, prioritizing actual processed outcome/control rows and compatible effect contrasts. Recent2026 PDAC may substitute if publicly processed and relevant; do not exceed4 resources. CELLFIE is primarily prior art, not default tumor-side dataset.
- [ ] If first resources cannot support a distinct transport contrast, inspect Dufva's matched drug outcomes as the ONE within-scope alternative. This replaces a remaining resource slot; it does not authorize another broad search branch.
- [ ] For each opened file record source/publisher link, retrieval UTC, size/hash, sheet/table/column names, measured vs computed fields, outcome and direction, treated/control definitions, perturbation side (tumor or effector), species, cancer model, antigen/CAR/CTL identity, assay time, independent preparations/donors/models, and nested guides/wells/timepoints. Inspect exact rows and missingness; do not infer fields from paper descriptions.
- [ ] Treat differential enrichment as a derived statistic of measured survival/selection, not direct percent killing. Untreated fitness may already be incorporated: establish normalization before subtracting it again. A selected hit list without all-tested genes/control information limits evaluation and generalization.

## Task2 — Specify the scientific question from the files

- [ ] Write at most2 concrete cards: biological question, practical biological importance, exact measured contrast, compatible comparator, unit of inference, controls, and source of uncertainty. Distinguish within-model mechanism, cross-model reproducibility, cross-study transport and human/clinical inference.
- [ ] For transport, specify exactly what moves across what settings. Do not call thousands of genes/guides thousands of independent tumors. KO/CRISPRa/drug exposure and murine/human ortholog effects are not automatically interchangeable. Record ortholog mapping source and ambiguities if needed; avoid speculative expression-based response labels.
- [ ] Do not demand a separate clinical cohort for a carefully scoped experimental model claim. Conversely, one line with four cultures cannot support patient/population generalization. A single well-controlled multi-model dataset can support a biological contrast if the precise question remains unpublished; fair evaluation is claim-specific.
- [ ] Candidate must change biological interpretation or an experimentally measured choice. Generic atlas construction, enrichment, more genes, replication alone or a marginal better prediction score is insufficient. A sharply defined negative result may pass if evaluable, important and distinct; do not tune until positive.

## Task3 — Close the closest-prior-art gap

- [ ] For strongest card, read relevant Methods/Results of origin studies plus closest later primary reanalysis identified by small PubMed queries. Start with `(CRISPR AND ("CAR-T" OR "CAR T")) AND (resistance OR sensitivity)` and `(CRISPR AND "T cell") AND ("cancer intrinsic" OR "tumor intrinsic")`; adjust only to solve a named gap within query budget. Broad key concepts are better than overconjunctive strings. Query count0 does not prove absence.
- [ ] Record what authors already did, exact comparison left untested, why it matters, and competing explanations. In particular conserved sensitizers, untreated-fitness filtering, IFN sign changes, antigen-retention resistance, NK-induced resistance, gene modules, effector-side multiphenotype screens, and drug immunomodulation have known precedents above.
- [ ] Retrieve truncated full texts only through returned continuation offsets if needed. Mark abstract, targeted sections, full main text and inaccessible supplement distinctly. No novelty declaration from abstracts alone or a paper's future-directions statement.

## Task4 — Candidate gate and decisive pilot specification

- [ ] Use separate columns/statuses for actual data access, distinct contribution, biological importance, and evaluation feasibility. Mark unresolved rather than rejected where evidence is incomplete. Do not silently change criteria after seeing rows.
- [ ] If a card passes, return an exact minimum pilot specification: files and required columns, fixed eligibility, primary contrast/direction, fair baseline, grouping/split appropriate to biological claim, handling of untreated fitness/assay selection/missingness, uncertainty method, leakage checks, and proceed/revise/stop rules. Include which result would falsify premise or make contribution modest. Do NOT execute pilot until parent review and new plan.
- [ ] If no card passes, return the best evidence-backed near-miss with ONE exact unresolved test that could change decision, or state no such test remains in this scope. No manufactured winner or generic promise to search more.

## Task5 — Compact delivery and verification

Create only: `STATUS.md`, `search_log.tsv`, `sources.tsv`, `data_cards.tsv`, `candidate_gates.tsv`, `decision.md`, `data_manifest.json`, `preservation_before.json`, `output_hashes.json`, small read-only `validate_artifacts.py`, relevant files under `data/`, and a minimal inspection helper if needed. Avoid elaborate reporting or redundant validators.

TSV schemas (tab delimiter, UTF-8, one header, escape embedded tabs/newlines):
- search_log: id,utc,query_or_url,total_returned_and_screened,cursor_or_limit,purpose,source_ids
- sources: id,pmid_or_doi,url,title,reading_scope,measured_outcomes,already_reported,remaining_question,limitations
- data_cards: id,source_ids,url,path,bytes,sha256,schema,conditions,units,measured_vs_derived,access_and_limits
- candidate_gates: id,question,importance,exact_increment,closest_sources,data_ids,data_gate,novelty_gate,evaluation_gate,unit,comparator,next_decisive_test,verdict

Manifest: relative stage file paths, bytes,SHA256,official URL for downloads; hash every relevant output last except manifest itself. Validator rejects malformed TSV widths, invalid/nonfinite JSON, missing/unlisted downloads, wrong file sizes/hashes, unsafe paths, and protected-input changes. Do not hardcode zero downloads or use check success as evidence of scientific validity.

Python runtime: `C:/Users/dylan.bagley/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe`; set `PYTHONDONTWRITEBYTECODE=1`. Existing installed libraries only; use UTF-8 stdout or ASCII-escaped diagnostics where necessary. Run new stage validator and prior cycle3 validator once at closeout; parent independently reviews source/schema/gates and reruns affected checks. Earlier original preservation checks need repetition only if older-file changes are suspected; original changes are not authorized.

Update compact live `opportunity_search_Oct4/STATUS.md` with stage, actual finding, next two actions, gate decision and blocker. Return concise results, paths, actual wall time, checks run, and exactly what remains unresolved. Parent independently decides next step; executor may not declare the whole substantive-extension goal achieved.
