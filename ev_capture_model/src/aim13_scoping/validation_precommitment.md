# Validation on newly found labelled targets (precommitment, Sept 29, 2026, BEFORE any expression score was computed for them)

## What this is
Two independent literature searches (run by sub-agents that were not shown expression data) looked for labelled targets outside the
44-gene panel. They found only a handful that meet the frozen labelling rules cleanly, because most clinical toxicity attributions are
hedged by the authors. The labels, their tiers, confidence and sources are frozen in `target_labels_validation.py`
(SHA-256 in `validation_precommitment.md.hashes`); nothing there is edited after scoring.

## Primary validation set (pre-specified now)
safe (5): MAGEA4, SSTR2, CCR4, CTAG1B, IL2RA.  toxic (2): GUCY2C, PSCA.   => 10 safe/toxic pairs.
Hedged toxic (5: TFRC, SLC44A4, MLANA, MAGEA12, EPHA2) are reported descriptively and never pooled into the primary test.
PMEL, SLC6A2, CD276, PRAME are excluded (ambiguous, non-antigen, or too early).

## Rule (frozen, not tuned): score = top_share + 1[blood is the largest lineage family]  (round 4A definition, M2)
Also reported: top_share alone. Data: the same HPA single-cell reference (154 cell types) used for the main panel.

## Test and honest power statement
Separation (AUC) with an exact permutation p. With 5 v 2 targets the smallest attainable p is 1/21 = 0.048, so this set cannot
provide strong evidence either way; it is a directional check that the rule was not overfitted to its 30 development targets.
A pooled 37-target figure is shown for completeness and is NOT independent (30 of the 37 were used in development).

## Predictions (written before scoring)
- V1: the frozen rule orders at least 8 of the 10 primary pairs correctly (separation >= 0.80).
- V2: the two toxic primary targets (GUCY2C, PSCA) score below the median primary safe target.
- V3: the 5 hedged toxic targets score, on average, below the primary safe average (separation of hedged-toxic vs primary-safe >= 0.60).

## Label-consistency sensitivity for the MAIN panel (declared now; a robustness check, not tuning)
The boxed-warning cross-check found: NECTIN4 (labelled safe) carries a boxed warning for fatal SJS/TEN; the ocular boxed warning appears on
F3 and FOLR1 (labelled toxic) and also on a BCMA ADC (TNFRSF17, labelled safe from CAR-T/bispecific evidence). Flips, chosen from the label
text and not from any score: F1 NECTIN4 -> toxic; F2 F3 and FOLR1 -> safe; F3 = F1 + F2; F4 TNFRSF17 -> toxic. Report the frozen rule's
separation under each. Labels in the repository are NOT changed.
