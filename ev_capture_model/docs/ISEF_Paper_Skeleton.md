# Current paper-planning outline — October 2, 2026

Advisory structure and checks only. Dylan writes submission text in his own words and verifies references directly. This replaces the September outline; its complete previous version is preserved in `strengthening_Oct2/before/`. No abstract, research plan, or poster prose is supplied here.

## Scope to settle before writing

- Main evidence: retrospective association in a selected panel of target-agent examples, using normal-cell-type RNA specificity and lineage.
- Supporting theory: response bounds under explicit assumptions, with separate monotone and arbitrary response classes. Not a clinical safety theorem and not a derivation of the score.
- Distinguish original frozen analysis, source-based sensitivity scenarios, and any future independently adjudicated analysis. Do not merge them into one headline.
- Keep the biological project. Relevant expansion is authorized. Prior instructions forbidding all expansion or label-sensitivity analysis no longer apply.

## Introduction planning

- Explain the biological question and why normal expression could contain information about unintended normal-tissue effects.
- Distinguish target, agent, study/regimen, and outcome. Avoid implying that antigen expression guarantees cell killing or that all therapeutic formats share one mechanism.
- Situate tau, SCAN-ACT, BayesTS, and Thomas et al. before claiming a contribution. Confinement is rank equivalent to tau at fixed reference dimension.
- One short origin-story element can connect the EV distribution work to the present question. No new nonbiology pivot is needed.

## Methods planning

- Identify reference version, units, cell-type partition, zero/missing policy, gene/isoform limitations, and the exact lineage mapping used in each historical analysis.
- State the historical score and orientation; distinguish the earlier home-lineage version used in locked predictions from later cleaned flags where relevant.
- Describe discovery-panel selection and literature-based endpoint rules. Do not claim outcome-blind development without dated evidence.
- Explain original 40-row TSV versus 44 operational rows and the dated export. Preserve original hashes; new hashes do not retroactively freeze amended labels.
- Describe separate clinical review and uncertainty categories; do not treat the assistant's review as independent blinded adjudication.
- Explain pairwise AUC, target-level bootstrap, paired differences, and shared-target dependence. Do not treat pairs, cell types, or cells as independent clinical examples.
- Specify monotone versus arbitrary response classes, fixed weights, common response, positive total weight, and the special zero-weight case. See the corrected mathematical note.

## Results planning

- Original analysis: 19 tolerated versus 11 toxic labels; AUC 0.913876; 191 concordant and 18 discordant pairs. Conditional, exploratory result.
- Confinement alone and tau: AUC 0.880383. Blood flag alone: 0.710526. Report paired differences and uncertainty rather than contrasting separate significance tests.
- Targeted source-sensitivity table: reference `audit_followup_Oct2/label_sensitivity.csv`.
- Broader uncertainty diagnostic: deleting seven documented concern rows gives 19 versus 4 and AUC 0.815789. This is not a completed strict-label cohort.
- Modality diagnostic: small and coarse, conditioned on original labels; no claim of causal adjustment or cross-modality validation.
- Autopsy comparison: independent comparator expression, reused clinical labels. Paired difference interval spans zero; no demonstrated superiority or independent clinical replication.
- Protein/reference changes: separate reference robustness from outcome validation; include weak RNA/protein agreement and partitioning limitations.
- New-target validation and locked predictions: retain negative/inconclusive results; do not relabel pre-lock outcomes as prospective resolutions.
- Mathematical results: exact mixture identity; two distinct frontiers; AUC insufficiency; separate versus simultaneous tissue constraints. Tests are implementation checks, not novelty evidence.

## Discussion planning

- Assess whether the observed association adds useful knowledge beyond established specificity and lineage ideas. Do not claim that possible alternatives have been ruled out merely because one control is nonsignificant.
- Explain endpoint ambiguity, selection effects, agent/regimen dependence, expression-to-protein mismatch, and finite reference coverage.
- Explain what survives correction and what remains unknown. Statistical robustness is not a substitute for valid labels.
- Identify independently reviewed agent/study outcomes and fair existing-method comparisons as the most relevant extension. More rows alone do not constitute independent validation.
- Avoid individual patient predictions, clinical recommendations, broad rules about safe targets, or quantified award probabilities.

## Student preparation

- Explain AUC in your own words, including why 0.914 is not 91.4% classification accuracy.
- Reproduce the 30-target count and explain why the operational file contains 44 rows.
- Explain why eight donors or 30 targets remain small independent sample counts despite many cells or pairs.
- Explain the difference between independent expression measurements and independent outcome validation.
- Work through one monotone/arbitrary counterexample without relying on software output.
- Explain the score's equivalence to tau and what contribution remains after acknowledging it.
- Be able to identify which choices and code came from AI and which you independently checked or developed. Update the AI disclosure in your own words with the Adult Sponsor.

Current evidence and remaining tasks: `strengthening_Oct2/Strengthening_Report.md` in the outputs root. Historical detailed numerical records remain in `Aim13_Results_Two_Comparator_Ceiling.md`, with correction notices.
