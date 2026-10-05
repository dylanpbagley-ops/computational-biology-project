# Round 6B addendum (written after the DEATH-based weights were computed, BEFORE any target was scored with them)

Date: 2026-09-29. Nothing about any target's score or label has been looked at since faers_precommitment.md was hashed.

What I saw: the measured death-fraction weights have Spearman rho = -0.01 with the mean hand weight of each system's cell types.
Reason (a concept, not a score): the hand weights encode IRREVERSIBLE loss (brain, islets, glomeruli, retina, eye surface all >= 0.7),
whereas P(death | report) encodes ACUTE LETHALITY. Dose-limiting toxicity in the labels includes severe NON-fatal, disabling injury
(ocular keratopathy, neuropathy, skin), which a death fraction under-counts by construction.

Addition (a second measured weight, fixed now): W_severe = P(death OR life-threatening OR disabling | report of the system's PTs), from openFDA
seriousnessdeath / seriousnesslifethreatening / seriousnessdisabling, same PT lists, same rescaling (divide by the maximum), same median default.
Every rule A1-A7 is scored with BOTH weightings. W_death stays the pre-declared primary; W_severe is labelled as added after seeing the weights.
Cumulative count of alternatives tried against the frozen rule after this round: 52 + 13 = 65 (A1-A4, A6, A7 with two weightings = 12, A5 hand-weight control = 1).
Decision rules unchanged (promotion needs separation >= 0.944 and family-wise p < 0.05 across ALL 13 variants; otherwise the frozen rule stays primary).
