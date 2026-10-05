# Endpoint v0.3 — working observation specification

October 2, 2026. Developed after historical outcomes were seen; not preregistration or adjudication. Supersedes working v0.2 rules only; preserves all historical labels.

Question: in a named human agent/cohort/window, is there clinically significant normal-tissue injury attributed to recognition of the intended antigen? Each event has separate grade >=3, injury-related cessation and death fields. Protocol DLT remains separate. Broader treatment-associated toxicity is recorded even when intended-target attribution is unresolved; it is a different outcome.

For working event_supported state, both qualifying clinical severity and direct_support or author_probable intended-target attribution must be documented. Author_possible, missing mechanism, alternative_attribution, or missing severity stay indeterminate for this endpoint. A fatal immune activation event alone does not establish injury due to normal-tissue target recognition. No review is independent when performed by the assistant exposed to the scores.

explicit_nonobservation requires a named endpoint/grade and an explicit source statement or zero count. No severe dry mouth is restricted to dry mouth; no DLT is restricted to protocol DLT; neither is zero all-organ injury. Indeterminate includes silence, tolerability summaries, mechanism gaps, inadequate monitoring or unresolved identity. Ineligible includes preclinical-only evidence for this human endpoint. These states are event-specific, not universal gene verdicts.

Use identical severity rules for blood and nonblood. Keep intended depletion, loss of function, clinical consequences, required care, recovery and cessation separate. Reversibility or replacement care does not erase injury; depletion alone does not prove severe injury. Conditioning, disease, prior therapies, payload and immune activation are competing causes where relevant.

Adequacy fields: endpoint-evaluable and treated denominators and their units; named endpoint and grading; surveillance; exposure/activity; observation window; censoring and missingness; supportive care; regimen/route/combination/isoform; and cohort overlap. No numerical adequacy cutoff is invented. Outside clinical judgment must resolve adequate-negative eligibility and causal attribution before any amended binary analysis. NR means not reported; not_extracted means unchecked; NA means not applicable; none is zero.

Sampling frame: repair all original 30 outcome targets and screen all other 14 operational rows, including efficacy failures and exclusions. Do not treat efficacy failure as injury absence. New clinical pilot uses the separate deterministic workload frame. All agent/cohort/event records are retained, including non-discordant findings and access gaps. Overlapping reports add information, not independent patients.

Gate: no amended_analysis_labels.tsv or primary performance claim until source/endpoint adjudication and explicit negative eligibility are defensible. If gate cannot pass, retain observational states and clearly named conditional historical analyses. Score direction and rule remain unchanged. Scores/history/assistant judgments are withheld from reviewer initial fields.
