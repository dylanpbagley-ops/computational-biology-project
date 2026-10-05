# Targeted review of the five computed transport flags

This is a short, explicitly bounded post-pilot review, authorized by the ongoing request to assess actual candidates. It replaces no frozen pilot rule. Parent owns this plan and independent verification; Luna performs literature implementation. No further numerical exploration or fitting is authorized here.

## Context and objective

Workspace: `C:/Users/dylan.bagley/Desktop/outputs`.
Stage: `opportunity_search_Oct4/cycle4_functional`.
The frozen exploratory pilot found five exact-symbol genes depleted in at least three of six mouse CTL tumor contexts but enriched with reported FDR <=0.05 in both Ramos in-vitro CAR-T conditions: Med16, Nprl3, Cbfa2t3, Etv5, Irf1. Parent independently checked 720 raw-count numerical values and all 1,021 workbook flag decisions; results agree. These are leads, not a substantial extension or a causal modality comparison.

Primary files: `pilot/protocol.json`, `pilot/summary.json`, `pilot/transport_flags.tsv`, `pilot/gene_effects.tsv`; read only. Earlier/later sign-consistent sensitization: Med16 in EMT6/4T1/Renca-HA (three), Nprl3 and Etv5 in two, Cbfa2t3 in zero. Irf1 has known pathway roles and opposing effects within the CTL panel. Prioritize Med16; do not pick a favorable gene merely because it passes a flag.

Parent verified Lawson main-text Results explicitly calls Med16 an in-vivo EMT6 sensitizer: https://pmc.ncbi.nlm.nih.gov/articles/PMC9014559/ . Ramos main-text searches found no Med16/Nprl3/Cbfa2t3/Etv5 mentions: https://link.springer.com/article/10.1038/s41467-023-43790-2 . Main-text absence is not proof of novelty; the workbook already publishes the CAR-T effects. Cross-study comparisons confound background, antigen, timing, effector modality, protocol, and scoring. No patient or clinical-treatment claim.

## Exact execution steps

1. Capture actual start UTC. Read the five flag rows, corresponding early/late effects, and existing prior-art review. Preserve all four frozen `prior_art/` files, pilot protocol and data.
2. Conduct at most TWO new PubMed queries, each screening at most TEN returned titles. Query A: MED16 plus CRISPR, T-cell, CAR-T or immune-response terms. Query B: NPRL3/CBFA2T3/ETV5 plus CAR-T or cytotoxic-T-cell terms. Record exact strings, timestamps, result totals, offsets, and screened counts. This deliberately extends the completed eight-query discovery cap to at most ten for specific new numerical leads; it does not authorize a new broad search.
3. Read at most FOUR relevant primary abstracts and ONE additional relevant primary Methods/Results article or excerpt. Reuse existing source evidence wherever adequate. No reviews as sole support, no inaccessible full-text inference, no bypassing access challenges. If necessary full text is unavailable, record it as unresolved.
4. For each flag, classify separately: source-study effect already published; exact cross-study comparison directly tested or not found in this bounded search; biological interpretation known, plausible but unverified, or unsupported; robustness/timing limitations. A missing search hit never proves novelty. Do not exclude every possible IFN/autophagy-related extension merely because a generic pathway is known.
5. Recommend whether MED16 merits one later, prespecified validation stage, or whether these data do not justify further extension work. The recommendation must identify what evidence would discriminate an artifact/known restatement from a material biological contribution, and whether that evidence is computationally accessible. Do not claim accessibility without inspecting it; do not download anything in this task.
6. Save only `cycle4_functional/flag_review/search_log.tsv`, `sources.tsv`, `review.md`, `STATUS.md`. Provide a concise message with counts, strongest result and remaining uncertainty to parent. No changes to live status, hashes or old files; main executor will close its inventory separately, and parent will request one final refresh if needed.

## Bounds and stopping rules

Use the existing isolated gpt-6-luna Max reviewer. No new agents. Deadline: **2026-10-05 05:42:00 UTC**, and stop earlier if exact prior art resolves the question. The whole round still ends by **05:52:35 UTC**. Do not start a new dataset hunt, predictor, pathway-enrichment exercise, simulation, experiment, contact or paper prose. Preserve original clinical labels and predictions. Report negative or ambiguous results candidly. A useful next-stage candidate is not an established substantial extension.
