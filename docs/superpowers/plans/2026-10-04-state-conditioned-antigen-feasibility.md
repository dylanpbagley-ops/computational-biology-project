# State-conditioned antigen selectivity: bounded feasibility plan

Executor: the existing isolated gpt-6-luna Max subagent. No additional agents. Parent reviews cycle2 and this plan before delegation. This is a feasibility test, not an adopted extension or a claim of novelty.

## Purpose and decision

Test one different, medically relevant question: **does a defined inflammatory perturbation cause nonmalignant cells to present peptides classified as tumor-restricted using resting-tissue evidence, and can a specific treatment-conditioned targeting rule be tested against independently measured presentation?**

Generic IFN-gamma effects on HLA, healthy-tissue filtering, peptide ranking, and tumor-versus-normal expression are established concepts. They cannot be the contribution. The exact potential increment must concern measured change in target selectivity across cell states, with a falsifiable contrast that the closest primary work did not already test. No causal toxicity or clinical-safety claim follows from peptide detection alone.

If this precise question is already answered or cannot be evaluated with accessible processed measurements, stop it. Do not manufacture an extension or continue a broad paper collection. One revision grounded in actually inspected measurements is allowed; speculative revisions are not.

## Paths, constraints, and clean context

- ROOT: `C:/Users/dylan.bagley/Desktop/outputs`.
- PLAN: `docs/superpowers/plans/2026-10-04-state-conditioned-antigen-feasibility.md`.
- STAGE: `opportunity_search_Oct4/cycle3`. All new downloads, records, helpers and outputs go here. One additional authorized live file is `opportunity_search_Oct4/STATUS.md`.
- Read cycle2 `decision.md`, `data_cards.tsv`, `sources.tsv`, `output_hashes.json`; do not read full older histories. Hash these references and PLAN before work.
- Existing Python: `C:/Users/dylan.bagley/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe`. Set `PYTHONDONTWRITEBYTECODE=1`; no installs or changes to older dependencies.
- Preserve every old stage, label, prediction and workbook. Original clinical association remains retrospective; this proposed antigen-processing direction would be a new biological extension, not a relabeling of original clinical validation.
- No contact, lab work, registration, publication, credentials, paid services or controlled-access applications. Dylan writes submitted prose. Only evidence, feasibility records and analysis outlines here.
- School fair late January, state mid-March; weekly hours unknown, no reviewer identified. Do not rely on future access or collaboration.

Verified starting evidence:

1. JCI 2024, PMID37906288, DOI10.1172/JCI166740, https://pubmed.ncbi.nlm.nih.gov/37906288/: untreated breast tumors; Table5 has74 selected antigens, peptide sequences, tumor-sample field, HLA, category, mTEC-RNA mean and HLA Ligand Atlas flag. DataValues has tissue expression-frequency summaries. These RNA summaries are not measured healthy peptide abundance. TMT was not used quantitatively. Some sample-list cells are decimal-coerced; do not invent sample IDs. Eleven peptides have peptide-pulsed healthy-donor assays; two donors, no donor key in workbook; clone rows are nested. Generic safety filtering and CAF origin are already studied.
2. Science 2025, PMID40339010, DOI10.1126/science.adk3487, https://pubmed.ncbi.nlm.nih.gov/40339010/: PDAC organoid immunopeptidomics, healthy-tissue immunopeptidomic/Ribo-seq filtering, 56 selected cryptic-peptide priming assays and selected engineered-TCR endogenous PDO recognition. The paper acknowledges unexamined normal cell types/states and limited healthy Ribo-seq tissue coverage. Its generic healthy filter and high/low-avidity recognition contrast are already published. Its PDOs were IFN-gamma treated for immunopeptidomics. A limitation or future direction is a lead, not evidence of novelty or feasible data.
3. Cycle2 found no pilot candidate. Do not reopen acute-regrowth/p27/AKT, old clinical cutoff/myeloma/venetoclax, EV calibration, Umetani or HIV pairing work.

## Step1: closest prior work before a candidate

- Start true UTC clock; create compact STAGE/STATUS.md and `preservation_before.json`.
- Prefer enabled PubMed connector for study searches, metadata and rights-eligible full text. Discover through ALL_TOOLS; use returned pagination offsets. Retrieve single PMCIDs in10k-character chunks if larger calls fail. Read relevant Methods/Results, not just abstracts; record reading level honestly.
- Use at most8 focused literature queries, then at most2 targeted citation/access gap fills. Focus on processed immunopeptidomics of **nonmalignant** cells before/after IFN-gamma or a defined inflammatory stimulus, and state-dependent loss of tumor-antigen selectivity. Also examine the closest prior method that already incorporates inflammatory healthy-cell presentation into target selection. Record exact queries/timestamps.
- Closely inspect2–4 decisive primary sources. Record what already measures the question, what remains untested, and whether a proposed rule merely reproduces authors' conclusions. Broad IFN-gamma remodeling is not enough.
- Do not add sources just to reach a quota. Preprints must be labeled; predictions cannot supply experimental labels.

## Step2: two-resource data test, before analysis

- Select at most2 new public processed resources. Prefer one genuinely paired nonmalignant perturbation experiment and one independent challenge/replication resource if available. Baseline healthy Atlas already used for candidate filtering is not a new independent validation dataset.
- Verify nonmalignant status, stimulus/dose/duration, HLA, sample/culture/donor ID, paired controls, biological versus technical repeats, measured peptide sequence and quantitative field if any, normalization, identification limits and original sample keys from actual methods/files.
- Download only smallest relevant processed CSV/TSV/XLSX/JSON tables; size-check first. Limits25MiB/file and75MiB total. No raw mass spectra, sequencing, images, whole repositories or portal dumps. Follow legitimate public publisher/archive redirects, including public signed downloads; no login or access-control bypass. One legitimate alternate route per blocked file, then record unresolved access and stop retrying.
- Inspect actual headers and bounded rows immediately. Record file source/size/hash; count biological units separately from peptide observations. Scan needed sheets/records only, not all million-row tables.
- Use JCI Table5 at `opportunity_search_Oct4/cycle2/data/JCI166740.sdt5.xlsx` read-only for an exact peptide-sequence join if scientifically applicable. Do not repair sample cells or infer HLA presentation from binding prediction alone. A positive direct healthy MS hit can challenge restriction; a missing hit cannot prove safety. I/L ambiguity, HLA compatibility, culture conditions and MS depth must be recorded.
- A small descriptive overlap/count to establish whether evaluable records exist is authorized; no fitted model, significance fishing, broad gene ranking or large analysis. Do not silently change the endpoint if the overlap is zero. Another exact biological question can be proposed only once, from inspected paired measurements, and needs closest-prior-work comparison.

## Step3: contribution and evaluation gate

Return `decision.md` with one verdict: `pilot_candidate`, `revise_with_data`, or `stop_current_question`. Treat originality, access and evaluability separately.

A pilot candidate needs all of these:

1. A precisely stated untested biological or targeting decision, beyond generic IFN-gamma remodeling or another peptide score.
2. Actual processed experimental rows that support it. Clinical outcomes and autologous patient T cells are not required for an appropriately limited ex-vivo presentation claim.
3. A fixed contrast and defensible baseline, for example resting-only evidence versus a defined state-aware constraint evaluated on **held-out measured** normal-cell presentation. Training/selection and evaluation must not reuse the same condition or source as independent evidence.
4. Valid unit/grouping design. Multiple peptides and wells do not create donors. A single cell model may support a model-specific mechanism but cannot establish population safety. State available units, missing units and justified scope; do not invent an arbitrary universal minimum sample size.
5. An explicit falsifying result, selection/observability handling, realistic uncertainty and biological importance. Do not use MS nondetection as a true negative or label detection as proven injury.

If a candidate passes, specify the smallest next pilot: exact files/columns/joins, inclusion/exclusion, endpoint, independent groups, fixed comparator, expected failure modes, and what conclusion could follow. Parent must review and write a separate pilot execution plan before fitting/running it. If it fails, say whether it is prior-art overlap, absent measurements in inspected files, unavailable files, weak independent replication or lack of a meaningful contribution. Do not default to clinical repair or call a bigger table an extension.

## Step4: compact artifacts and preservation

Produce only what supports the decision:

- `STATUS.md`: timings, verdict, key evidence, exact next step; at most500 words.
- `sources.tsv`: `source_id,title,year_and_status,doi_or_accession,url,reading_level,sections_read,measured_variables,derived_variables,prior_art_overlap,exact_unanswered_contrast,data_links,limits`.
- `data_cards.tsv`: `data_id,source_ids,official_url,access_level,file_path,bytes,sha256,schema,stimulus_and_timing,units_and_counts,join_keys,paired_controls,measured_vs_derived,limits`.
- `search_log.tsv`: `search_id,utc,purpose,query_or_url,source_ids,next_action`.
- `decision.md`, `data_manifest.json`, `preservation_before.json`, `output_hashes.json` and a small `validate_artifacts.py` adapted from cycle2. Record finite JSON, exact TSV widths/IDs/references, all relevant new output hashes and source download hashes. Exclude manifest itself and explicitly documented bytecode/cache only; no zero-data hardcodes or claims derived only from hardcoded counts.
- Any necessary tiny schema/join helper and its output stay in STAGE and are inventoried. Do not build a research software framework.
- Update authorized root `opportunity_search_Oct4/STATUS.md` with latest cycle and honest verdict.
- Rerun cycle2 validator and original `longitudinal_extension_Oct3/preservation/verify_preservation.py` read-only. The two documented old status-suffix reconciliations are pre-existing; do not restore old snapshots.
- Freeze final outputs, hash last, run stage validator, return exact results and constraints. Parent independently reviews source claims, actual schemas/counts, hashes and preservation before selecting any next plan. Parent assigns corrections back to this same Luna.

## Bounds and stopping rule

Target35–50minutes wall; research freeze at50minutes, hard75-minute checkpoint to finalize candid evidence. Log actual UTC/wall; active time only if measured. Send one meaningful checkpoint by25minutes.

Stop resource hunting after2 resources or when closest work already answers the exact question. Do not spend the budget repeatedly retrying blocked sites. If no candidate survives and no different data-backed next test exists, explicitly recommend stopping this search direction. Iteration is justified only when a concrete test can change the decision; successful file checks and more reading alone are not a substantial extension.
