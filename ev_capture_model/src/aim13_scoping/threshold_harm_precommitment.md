# Round 5B: a rule derived from the theorem's threshold form (precommitment, Sept 29, 2026, BEFORE running)

## The gap this addresses (found while running round 5A, not from scoring these rules)
Section 29.3 writes harm as H ∝ Σ_c w_c s_c (weight x share of the gene's expression) and says top_share is its "leading term". But if all weights are equal,
Σ_c s_c = 1 for every gene, so this linear form does NOT reward confinement at all; with binary weights it depends only on how much expression sits in expendable
cell types (that is exactly measure V3, which scored 0.699), not on how few cell types carry the antigen. The confinement effect does follow from the threshold form of the
same harm: a therapy kills a cell when its antigen level exceeds a threshold t, so H(t) = Σ_c w_c π_c q_c(t), where π_c is the cell type's share of the body and q_c(t) the fraction of
its cells above t. More cell types above threshold means more terms, so H grows with breadth. Round 4B/4C's count measures (D1-D3) already used this idea without weights.

## Schemes (fixed; higher = predicted safer; x_c = nCPM of cell type c; a cell type is "above" when x_c >= t)
TH1(t) = -(number of VITAL cell types above t)            [binary weights: vital = not blood, blood_other, germline, placental]
TH2(t) = -(sum over cell types above t of the hand weight w_c)   [expendability_weights.py, default 0.70]
for t in {1, 10, 50} nCPM  ->  6 schemes (cumulative 52; C3/tau duplicate noted earlier). No cell-abundance term (not measured for the 154 types), stated as a limitation.
Family-wise max-statistic permutation over these 6 plus the frozen rule (2,000 shuffles). Validation: same schemes scored on the validation primary set.

## Predictions (before running)
- S1: TH1 and TH2 land between 0.78 and 0.91 on the 30-target panel (near the D-family and V-family results).
- S2: none reaches 0.944 with family-wise p < 0.05.
- S3: TH2 is not better than TH1 by more than 0.03.
## Decision rule (as before)
A scheme is promoted only if separation >= 0.944 and family-wise p < 0.05, and would still need confirmation on new targets.
