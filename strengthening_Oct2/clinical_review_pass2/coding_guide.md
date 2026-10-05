# Clinical evidence coding guide v0.1

Working extraction specification, October 2, 2026. Developed with knowledge of past outcomes and during source review; not a preregistration. The fields below are proposals for review, not adjudicated labels.

## Unit and eligibility

Ultimately use one agent–study/cohort–event record, linked to publication, trial, target, and treatment context. This pilot still contains some event-family summaries; split these before event-rate analysis. Record human versus preclinical evidence. Restrict a clinical primary analysis to human observations, with distinct preclinical supporting evidence. Keep antigen isoform and gene mapping explicit. Record multicomponent agents and combination settings rather than inventing a single-gene mechanism.

Keep severity, protocol-defined DLT, treatment cessation, death, normal-tissue site, and mechanism attribution in separate fields. DLT is a study-defined category and is not synonymous with grade 3 or greater. A serious event need not be high grade, and a high-grade event need not be a protocol DLT. Regulatory approval, efficacy, and program termination do not determine the outcome label.

## Attribution field

| Code | Meaning | What it does not mean |
|---|---|---|
| direct_support | Clinical observations include direct tissue/mechanistic support for normal-tissue target recognition | Universal causality across agents or settings |
| author_probable | Primary authors judge the mechanism probable | Independent adjudication or proof |
| author_possible | Authors propose or cannot exclude the mechanism | Established on-target injury |
| alternative_attribution | Authors attribute the event to another cause | Proof that this target is safe |
| unresolved | Available source does not establish the mechanism | Absence of injury or a negative label |
| no_event_observed | Explicit absence of the specified event in this context | Zero risk or adequate exposure in every patient |

`no_event_observed` requires a named outcome: for example, no DLT, no liver toxicity, or no observed on-target/off-tumor injury. These three statements cannot substitute for one another. Abstract silence is `unresolved`, not `no_event_observed`. Confidence in a source's wording is separate from confidence in its causal inference.

## Endpoint proposals to resolve before amended scoring

The closest match to the biological question is **observed clinically significant normal-tissue injury attributed to target recognition in a specified treatment context**. Proposed positive criteria: human evidence; grade 3 or greater injury OR injury requiring cessation; and direct support or author-probable attribution. Preserve possible attribution in a separate sensitivity set. This is a proposed decision rule, not an applied clinical label. A reviewer should decide whether that composite is coherent and whether cessation reflects injury rather than another reason.

Protocol DLT of any cause should be a separate descriptive endpoint. It can document clinical tolerability but does not specifically validate an expression-based normal-tissue mechanism. Intentional blood-lineage depletion, recovery, supportive care, and clinical consequences require their own fields under the same severity rule; the model's blood flag must not determine how outcomes are coded.

The negative category is the main unresolved issue. An explicit no-event report is necessary but may be insufficient: exposure, follow-up, monitoring, treatment activity, and observation window matter. Do not invent a minimum sample-size cutoff after seeing which genes survive. Ask a reviewer to specify an adequacy criterion, or retain an observational three-way outcome (event / explicitly not observed / indeterminate) and avoid describing it as intrinsic safety.

## Reviewer questions

1. Is the proposed normal-tissue injury endpoint clinically coherent? Which evidence is sufficient for probable on-target attribution?
2. What constitutes adequate observation for a negative study record? Should no-event observations remain a separate descriptive category?
3. How should intentional hematopoietic depletion and managed but severe effects enter the same outcome definition?
4. Which records require antigen-isoform matching or cannot support a gene-level expression comparison?
5. Are records within a study separate cohorts or repeated reporting of the same patients? Which publications overlap?

Return initial judgments without viewing model scores or the assistant's provisional classifications. Log source location, reasons, uncertainty, reviewer identity/date, disagreements, and changes. After reconciliation, freeze a separate amended label version; retain all original artifacts.
