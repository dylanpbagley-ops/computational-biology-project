# Original-panel clinical evidence review: pass 3

October 2, 2026. Advisory analysis, not student submission text.

**The representative source pass now covers all 30 original outcome-panel targets.** This pass adds 25 target/source records to the five original targets covered in pass 2. Human reports were located for 29 targets; the cited FAP injury evidence is preclinical. This does not mean all agents, studies, events, or mechanisms have been exhaustively reviewed. Many records still rely on abstracts or selected primary-text sections.

The most useful new result is a demonstrated inconsistency in the historical outcome rationales. Both historically favored and unfavorable examples need the same definition of severity, dose limitation, target attribution, and adequate observation. Approval and descriptions such as manageable or tolerable cannot supply that definition.

## Concrete findings

| Target | Verified observation | Consequence for the research |
|---|---|---|
| GPRC5D | Talquetamab's phase I report includes a grade 3 rash classified as dose-limiting | The historical blanket no-DLT skin-effect rationale is incorrect; causal attribution still needs review |
| NECTIN4 | EV-101 reports severe skin events and describes skin toxicity as likely on-target; formal DLTs were different events | A severe-injury endpoint and a protocol-DLT endpoint give different questions, even for an approved treatment |
| CD47 | The trial reports on-target anemia, while its DLT definition expressly exempts grade 3 anemia | Severe on-target harm cannot be inferred from, or excluded by, DLT status alone |
| GPC3 | The clinical report includes one fatal CRS event | The old tolerability summary is incomplete; fatal CRS alone does not establish normal-tissue target recognition |
| MSLN | One regional-treatment study explicitly reports no observed on-target/off-tumor toxicity and no DLTs | This is a useful context-specific negative observation, subject to adequacy review |
| FOLR1 | Human corneal IHC lacked FRalpha; supporting animal findings were consistent with off-target corneal effects | The old corneal on-target rationale is unsupported; this does not make the gene universally safe |
| FAP | The cited cachexia/bone injury occurred in mice | It cannot be a human clinical outcome in the same panel |
| CEACAM5 | The colitis study used TCR-T, not CAR-T | Agent identity and modality comparisons require the correction |

Primary sources: [GPRC5D](https://pubmed.ncbi.nlm.nih.gov/36507686/), [NECTIN4](https://pmc.ncbi.nlm.nih.gov/articles/PMC7106979/), [CD47](https://pmc.ncbi.nlm.nih.gov/articles/PMC7186585/), [GPC3](https://pubmed.ncbi.nlm.nih.gov/32371538/), [MSLN](https://pmc.ncbi.nlm.nih.gov/articles/PMC8563385/), [FOLR1](https://pubmed.ncbi.nlm.nih.gov/30413525/), [FAP](https://pubmed.ncbi.nlm.nih.gov/23712432/), [CEACAM5](https://pmc.ncbi.nlm.nih.gov/articles/PMC3048186/). Exact source locations and remaining gaps are attached to individual extraction rows.

## What is ready, and what is not

Factual corrections above are ready to use in working interpretations. Original artifacts are retained for auditability, not endorsed as corrected truth. The separate clinical-event axes can be reviewed now. Final binary outcome assignments are not ready, and no amended AUC is reported.

The difficult issue is not merely finding reports of side effects. Most abstracts do not establish whether a particular injury resulted from recognition of the target on normal tissue. They also do not provide a defensible negative label merely by omitting such a finding. Calling all such cases safe would repeat the original problem in a larger table.

For blood-lineage targets, record intended depletion and its consequences without letting the model's blood flag determine the label. An additional primary pooled analysis reports B-cell aplasia in responders and widespread immunoglobulin replacement after tisagenlecleucel. It overlaps ELIANA, so it is supporting follow-up rather than an independent new cohort. [Pooled safety report](https://pubmed.ncbi.nlm.nih.gov/34353848/).

## Files and use

- `extraction.tsv`: the 25 new representative target/source records, including explicit missing information. Event families are not individual-event incidence data.
- `original_panel_coverage.tsv`: all 30 targets, source links, and review priorities. It contains no new outcome assignments.
- `source_only_review_form.tsv`: a consolidated 30-target form with source links and blank judgments. It omits model scores, historical outcomes, and assistant interpretations. Split into agent/study/event rows during review. This replaces neither the earlier forms nor their audit history.
- `remaining_work.tsv`: one concrete next task per original target.
- `checks.json`: coverage, schema, blank-form and historical-integrity checks; not clinical validation.

Use the existing [reviewer instructions](../Reviewer_Instructions.md) and [coding guide](../clinical_review_pass2/coding_guide.md) with the source-only form. Withhold this report, the assistant extraction, old labels and scores until the reviewer's initial assessment. This is a purposive source set assembled with prior outcome knowledge, not a blinded systematic review.

## Next research step

Complete event-level attribution and observation-adequacy extraction, prioritizing GPRC5D, NECTIN4, CD47, the ADC neuropathy/ocular examples, and managed blood-lineage depletion. These have the greatest potential to change interpretation. Resolve antigen matching for CD44v6 and CLDN18.2. Treat multiple reports on one target as related observations.

Then freeze a coherent amended endpoint and apply the existing score alongside confinement/tau and blood-lineage baselines, without fitting a new rule to rescue performance. Keep possible attribution and uncertain negative observations in explicit sensitivity analyses. If a defensible binary panel cannot be assembled, report the target-level association under its original operational labels and the evidence limitations separately; do not manufacture a replacement ground truth.

This remains an improvement of the same project. The work so far strengthens traceability and identifies a fairer test; it does not yet establish that the clinical association will survive, prove novelty, or guarantee competition/publication outcomes.
