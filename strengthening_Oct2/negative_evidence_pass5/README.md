# Favorable-case evidence and endpoint consistency: pass 5

October 2, 2026. Advisory research analysis, not student submission prose.

**The historical favorable class cannot yet be treated as a set of adequately observed negatives for target-mediated normal-tissue injury.** This is a limitation of the evidence assembled here, not a finding that every treatment is harmful or that better evidence cannot exist.

## Completed work

- Triaged all 19 historical favorable targets against the same question, using prior extraction for coverage and fresh primary-source checks for six targets. Fresh checking was targeted, not an exhaustive literature review.
- Made the observation categories explicit in `endpoint_contract.md`: supported event, explicit nonobservation, indeterminate, or ineligible, at the agent/cohort/event level.
- Separated normal-cell depletion, clinical consequences, supportive care and recovery. Manageability does not erase injury; depletion does not by itself prove severe injury.
- Built a 30-target review template with blank judgments and explicit fields for denominators, monitoring, observation windows, confounding and overlapping cohorts. No scores or historical outcomes appear in that template.

## What the primary-source follow-up adds

The [CD19 pooled report](https://pmc.ncbi.nlm.nih.gov/articles/PMC8344270/) documents on-target B-cell depletion with extensive replacement care. It overlaps the original ELIANA cohort. It supports evaluating consequences under supportive care, not treating normal-cell loss as absent or adding an independent validation cohort.

A [separate BCMA study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9153033/) supports prolonged immune depletion. It selected ongoing responders and is not assumed to represent the historical ide-cel/cilta-cel agents. Infection-event counts must not be mistaken for affected-patient counts or attributed wholly to one mechanism.

The [mesothelin trial](https://pmc.ncbi.nlm.nih.gov/articles/PMC8563385/) supplies an explicit context-specific nonobservation. Regional administration, variable persistence and subsequent combination treatment matter; full negative-observation adequacy is still incomplete.

The [VISION report](https://pmc.ncbi.nlm.nih.gov/articles/PMC8446332/) supplies an explicit zero for severe dry mouth, with a known safety denominator, grading system and treatment-emergent observation window. That is useful evidence for a narrow endpoint. It cannot establish absence of all normal-tissue injury, and median survival follow-up cannot substitute for the toxicity window.

The [rituximab abstract](https://pubmed.ncbi.nlm.nih.gov/9704735/) does not settle target-specific negative-observation adequacy. Detailed immune recovery and infection tables remain a concrete retrieval task.

The [final ALFA-0701 report](https://haematologica.org/article/view/8727) documents treatment cessation for thrombocytopenia and severe bleeding, despite a favorable benefit/risk conclusion. These are distinct event groups. The results undermine using possible marrow recovery as evidence that consequential injury was absent; they do not alone establish the precise target-recognition mechanism.

Exact observations, scope and source locations are in `source_notes.tsv`. Some direct pages were intermittently inaccessible; the table distinguishes inspected abstracts/indexed primary text from full-page access. No secondary review was used to fill missing primary evidence.

## Consequence for the analysis

The triage contains two explicit nonobservations with different scopes, two depletion-related records (one from another agent/cohort), one provisional supported severe event, one severe event with unresolved attribution, and 13 indeterminate cases. These are workflow statuses, not new binary clinical labels. They do not count adequately established safe targets.

This extends uncertainty to cases held fixed in pass 4. Therefore even its broad 0.517–0.924 range was not a bound on every source of clinical-label uncertainty. No replacement AUC is calculated here: choosing a favorable subset or treating unknowns as negatives would recreate the problem.

An important possible bias is now explicit: the historical label rationales accept managed hematopoietic depletion while the model rewards hematopoietic lineage. That creates a risk of definition-driven association. It is not proof that the numerical finding is entirely artificial. The remedy is consistent consequence-based coding and independent review, not automatically changing every blood-lineage label to unfavorable.

## Next concrete dependencies

1. Complete full event/observation records for the 13 indeterminate favorable cases, starting with blood-lineage depletion and ADC mechanism attribution; do not repeat abstract-only coverage as if it resolves the gap.
2. Resolve the named-agent BCMA/CD19 distinctions, CD44v6/CLDN18.2 antigen mapping, and overlapping reports before any study-level expansion.
3. Have an independent reviewer apply the common endpoint with scores withheld. `review_form_score_withheld.tsv` is a target-indexed starting template: duplicate rows into distinct agent/cohort/event records. It retains earlier source links and adds this pass's inspected sources; further sources must be screened consistently, not selectively used to confirm the assistant's triage. Supply the endpoint contract and prior reviewer instructions, but withhold this interpretation and the triage until initial judgments are complete. Source selection itself was not blinded.

The assistant can keep extracting evidence, but cannot supply independence by reviewing its own work again. No external reviewer has been contacted, and no reviewer judgment has been fabricated. Keep the current model and historical outputs intact while this evidence work proceeds.

## Verification and assistance record

`build_and_check.py` verifies exact coverage against the original label table, source-note schema, blank review judgments, and hashes of protected historical files and prior-pass artifacts. `checks.json` records results. These checks verify records and preservation, not clinical truth. No model parameters, original labels or predictions were changed. This pass is an AI-assisted, post-discovery evidence audit with explicit missingness; it is not a preregistration or a systematic review.

Access log: the rituximab primary-paper copy in a USPTO filing was located, but subsequent section retrieval returned 403; no full-text immune/infection counts from that copy were extracted. For CD33, the accessible publisher's final primary report was used instead of secondary summaries. Its patient groups and ascertainment differ from the initial paper; they remain identified by publication. No additional clinical judgments were inferred from access failures.
