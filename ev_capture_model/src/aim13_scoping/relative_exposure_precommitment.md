# Round 6C precommitment: relative-exposure harm (written after seeing round 6B's results, before running anything below)

Date: 2026-09-29.

## What round 6B showed (already seen)
All 13 measured-severity variants scored BELOW the frozen rule: absolute-threshold worst-organ rules A3/A6/A7 got 0.70-0.72 (paired difference vs frozen about -0.2 with
CIs excluding 0), the sum rule A4 0.79-0.80, the share-weighted A2 0.82-0.84, A1 0.73-0.78, hand-weight control A5 0.80. In every case the measured weights were
indistinguishable from shuffled weights (shuffled-weight p 0.28-0.84). Nothing promoted.
Measured severities are uncorrelated with the hand weights (Spearman ~0.0; W_death vs W_severe 0.90).

## Why absolute thresholds are the wrong comparator (a reasoning step, not a score)
The theorem's harm is measured against the TUMOR comparator: a normal cell is at risk when its antigen level is a meaningful fraction of the level the therapy is
built to hit. An absolute cutoff (1-100 nCPM) says "expressed in a vital organ", but safe and toxic targets alike are expressed in vital organs; what differs is how
that normal level compares with the tumor level. Per-target tumor levels are not available here, so two bracketing assumptions are possible:
 (i) tumor level is the same for every target -> absolute thresholds (round 6B, failed);
 (ii) tumor level tracks the gene's own normal peak, with a floor for tumor-specific genes -> RELATIVE thresholds (this round).

## Relative exposure (fixed)
x_ref(g) = max( peak nCPM of gene g over the 154 cell types, T0 ) with T0 = 10 nCPM (the floor stops a gene with essentially no normal expression, e.g. DLL3,
from looking as exposed as a highly expressed gene).
e_c(g) = mean over r in {0.03, 0.1, 0.3} of 1[ nCPM_c(g) >= r * x_ref(g) ]  (r spans the modality range: T-cell engagers/CAR-T act on very low antigen, ADCs need more).
Harm from a gene, per weight set w over cell types:  MAX form  H = max_c w_c e_c   (worst-organ, weakest link)   or   SUM form  H = sum_c w_c e_c.
Score = -H (higher = predicted safer). This form contains both ingredients of the frozen rule: the home cell type is always exposed (e = 1), so H includes its
expendability weight, and additional cell types above r * x_ref add breadth.

## Variants (9)
Weight sets w_c (cell level): BINARY (0 for the frozen haematopoietic list, 1 otherwise), HAND (expendability_weights, unlisted 0.70), DEATH and SEVERE
(measured system weights of round 6B assigned to each cell type through faers_severity.CELL_SYSTEM; unmapped = median). Each with MAX and SUM = 8 variants.
RE9: system-level, MAX, W_death, exposure = mean(relative RNA exposure, relative IHC exposure), where relative IHC exposure of a system = its level / the gene's highest
level across all IHC rows (only genes with IHC data; n reported).
Robustness only: RE-max/W_death with T0 in {3, 30}.
Diagnostics (not adoptable rules): (a) RE-sum with 28 system-level weights fitted in-sample (upper bound), a permutation null of that optimum, and nested
leave-one-target-out fitting; (b) the same for the MAX form.

## Pre-declared "most realistic rule" (fixed independent of score)
RE-max with W_death (worst-organ, relative-to-tumor-proxy exposure, thresholds integrated over modality range, measured severity). RE9 adds protein evidence for the genes
that have it. Both are recommended as the mechanistic description regardless of score.

## Predictions (fixed)
P1 the BINARY-SUM variant lands within 0.05 of the frozen rule (the formula reproduces the frozen rule's information).
P2 no continuous weight set beats BINARY in the same formula by >= 0.03 (weights are not the bottleneck).
P3 nothing meets promotion: separation >= 0.944 and family-wise (max over the 9 variants) shuffled-weight p < 0.05.
Cumulative count of alternatives tried against the frozen rule after this round: 52 + 13 + 9 = 74.
Reporting rules as before: a lower score is not "worse"; a higher score on 30 targets is not "better" without the promotion test; frozen rule stays primary unless promoted.
