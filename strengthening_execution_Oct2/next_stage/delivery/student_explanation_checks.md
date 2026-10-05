# Student explanation checks

Before presenting, the student should be able to explain in their own words:

- Why cells are nested in putative participant codes and why an exact code-level holdout is different from an independent clinical cohort.
- Why the pilot-exposed Healthy codes fit the reference while pilot-exposed MM codes were excluded from the primary expansion.
- How `10,000 × target molecules / all-gene molecules` is computed, what integer molecule counts measure, and why RNA detection is not protein abundance or tissue-wide safety.
- What pooled, equal-code, and donor-envelope budgets constrain; why a lower positive fraction under the envelope does not by itself mean the rule is better.
- Why GPRC5D had identical thresholds, why TNFRSF17's aggregate matched-cell gap was robust to the stated worst-case unmatched-cell bounds, and why SLAMF7's evidence is weaker.
- What the 2,697 unmatched selected annotations mean, why they were not imputed as zero, and why the missingness endpoints are bounds rather than confidence intervals.
- Why `Hospital.Code` is treated as a leakage-protection grouping variable but not a globally verified identity; what `Amp.Batch.ID` and `Seq.Batch.ID` describe and why their presence does not separate batch from biology.
- Why the 668 Healthy/Malignant source-label records (659 joined in the locked frame) remain a separate uncertain class rather than being relabeled.
- What the local 5,189-UMI thinning diagnostic measures, how its 100 Monte Carlo draws differ from donor-level uncertainty, and why depth sensitivity does not establish true biological absence.
- Why the clinical evidence table preserves source counts/percentages but leaves all `observation_state` values indeterminate; why a reported zero or `NA` is not an adjudicated negative safety endpoint.
- Why the final contribution decision is “no demonstrated strong extension” despite successful retrieval and analysis, and what evidence would be needed to change that decision.

These are understanding checks, not a script or submission text. The student should revise the project explanation and any school-approved disclosures with their advisor.
