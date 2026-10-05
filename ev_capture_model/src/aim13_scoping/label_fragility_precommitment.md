# Label fragility analysis: precommitment (written Sept 29, 2026, before running)

Question: how many of the 30 safe/toxic labels would have to be wrong to erase the frozen rule's result?
No rule, weight or score is changed. Score = top_share + lineage-level blood flag (M2), read from flag_cleanup.json.
This does NOT count toward the 74 alternatives (labels move, the rule does not).

Analyses
F1. Exhaustive flips: for k = 1..6, every set of k targets (of 30) has its label swapped (safe <-> toxic).
    Report per k: number of sets, min / median / max separation, fraction of sets with p < 0.05,
    and the fragility index = smallest k at which ANY set gives p >= 0.05 (the adversarial worst case).
F2. Weak-label set (declared before running; from the Round 7 note): toxic->safe for EGFR, EPCAM, NCAM1, FAP,
    CD47, ERBB2, CEACAM5; safe->toxic for NECTIN4. All 2^8 = 256 subsets applied. Report separation and p for
    every subset, and the all-eight-flipped case.
p-value: Mann-Whitney normal approximation with tie correction and continuity correction (exact enumeration is
infeasible for thousands of relabelings). Validated against the stored permutation p for the frozen labels
(3.5e-5) and against a 20,000-shuffle Monte Carlo on 30 random flip sets; a mismatch of more than a factor of 3
in p on the frozen labels would be reported and the approximation flagged.
Separation = fraction of safe/toxic pairs where the safe target scores higher, ties half.

Predictions (written now)
P1. With 1-2 flips, separation stays >= 0.80 and p < 0.05 in every set.
P2. The fragility index is at least 4 (i.e. no 3 flips erase significance).
P3. The all-eight weak-label flip (F2 worst case) drops separation below 0.75, since it reverses 8 of 30 labels.
Reporting rule: all outcomes reported whether or not they match. No label is changed in the project files.
