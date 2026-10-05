# Round 4A: clean up the blood flag so it does not depend on one arbitrary "home" cell type (precommitment, Sept 29, 2026, BEFORE running)

## Problem being fixed
The section 23.2 rule adds +1 when the single top-share cell type is on the frozen blood list. For a gene spread over dozens
of cell types the top one is close to arbitrary, so a broad gene can collect the bonus by luck (AKAP10 in the discovery
screen). Round 2 showed a *soft* share measure loses the categorical signal, so keep a yes/no flag but define it from the
whole distribution.

## Schemes (fixed; both reported)
Notation as before: s_c = share of the gene's expression in cell type c; top_share = max s_c; h_c = 1 if c is on the frozen
blood list (target_labels.haematopoietic, NOT edited); B = sum_c s_c h_c = blood share.
| id | score | flag meaning |
|---|---|---|
| M1 | top_share + 1[B >= 0.5]                | majority of the gene's expression is in blood lineages |
| M2 | top_share + 1[blood is the largest lineage] | blood has a larger aggregate share than any other lineage family (`lineage_map.py`: 10 families, fixed) |
Baselines re-reported: one-factor 0.880; binary home-cell rule 0.914.

## Predictions (written before running)
- P1: M1 lands between 0.85 and 0.914. Reason: CD22, CD33, CD38 (the blood targets that carry the +0.033 gain) have low top
  shares (0.13-0.27), so a majority-of-expression test may not flag all of them.
- P2: M2 >= M1, because plurality is a weaker requirement than majority.
- P3: neither exceeds 0.914 by more than 0.01.

## Decision rule (a CLEANUP rule, not a "find a better score" rule)
The cleaned flag replaces the home-cell definition if BOTH: separation >= 0.900 AND the paired bootstrap interval of the
difference from the home-cell rule includes 0. Otherwise the home-cell rule stays and its arbitrariness for broad genes is
recorded as a limitation. No claim of "better" is made unless the paired interval excludes 0 on the positive side.
Also reported: the 30 targets' B, home lineage and flag under each definition, the flag counts (safe/toxic, Fisher exact),
a shuffled-flag null (blood flags reassigned across the 154 types, same count), leave-one-target-out minimum.
Cumulative look count after this round: 14 alternative schemes.
