# Round 4C: theory-consistent "vital-mass" measures and absolute-threshold breadth (precommitment, Sept 29, 2026, BEFORE running)

## Why these (reasoning from round 4B, not from looking at these scores)
Round 4B showed confinement measures other than top_share do no better, and total abundance alone gets 0.761. The theorem's harm is a
sum over cell types of (weight x normal expression), with weight 0 for expendable types. Two things follow that 4B did not test:
(1) abundance should be counted only in VITAL types, so that expendable expression (blood, germline, placenta) costs nothing,
which handles blood-lineage targets without a separate flag; (2) a therapy responds to antigen above a threshold, so breadth
should be counted at an absolute level, not relative to each gene's own peak (that relative version was C5).
Expendable lineages (fixed by principle: marrow-regenerated with a clinical replacement, or fertility/pregnancy-only):
blood, blood_other, germline, placental (from `lineage_map.py`). Everything else is vital.

## Schemes (fixed; higher = predicted safer; each alone and "1000*flag + score" with the round 4A flag M2)
| id | score |
|---|---|
| V1 | -log10(1 + sum of nCPM over vital cell types) |
| V2 | -log10(1 + highest nCPM over vital cell types) |
| V3 | share of expression in expendable lineages |
| V4 | mean of panel-rank of top_share and of V1 |
| D1 | -(number of cell types with nCPM >= 1) |
| D2 | -(number of cell types with nCPM >= 10) |
| D3 | -(number of cell types with nCPM >= 50) |
14 new schemes plus baselines C0 and C0+flag; family-wise max-statistic permutation over these 16 (2,000 shuffles).

## Predictions (written before running)
- R1: V1 beats total abundance C6 (0.761) and lands between 0.75 and 0.914.
- R2: V3 scores below 0.80 (germ-line expression such as CA9's spermatid signal will read as expendable).
- R3: D1-D3 land within 0.05 of C5 (0.876).
- R4: no scheme reaches 0.944 with family-wise p < 0.05.

## Decision rule (unchanged from 4B)
Promote to "hypothesis to confirm on new targets" only if separation >= 0.944 AND family-wise p < 0.05. Nothing enters the paper
from these 30 targets alone. Cumulative alternatives tried after this round: 44 (30 + 14), of which tau (C3) was rank-identical
to top_share, so 43 distinct.
