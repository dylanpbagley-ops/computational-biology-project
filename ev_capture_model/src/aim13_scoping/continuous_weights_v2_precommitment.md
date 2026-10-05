# Continuous expendability, second attempt: precommitment (written Sept 29, 2026, BEFORE any of the schemes below was run)

## Why a second attempt
The first continuous scheme (`continuous_weight_sensitivity.py`) scored 0.809 vs binary 0.914 and did not beat
its shuffled-weight null (p = 0.15). Diagnosis (already recorded): it weights only the *home* (top) cell type,
and for broadly expressed targets the home is close to arbitrary (EGFR -> placental tissue). The theorem (section
29.3) does not use the home cell type at all: harm is  H = sum_c w_c * s_c , where s_c is the share of the
antigen's total normal expression sitting in cell type c. So the principled fix is to use ALL the shares, not one.

## Schemes fixed in advance (no other scheme will be tried; all will be reported)
Notation: s_c = share of the gene's expression in cell type c (from the 154-type HPA single-cell consensus, the
same reference as section 23); top_share = max_c s_c; w_c = the hand-assigned weight in expendability_weights.py
(unlisted types = 0.70, as before, not changed); h_c = 1 if cell type c is haematopoietic (target_labels.haematopoietic), else 0.

| id | score (higher = predicted safer) | what it is |
|---|---|---|
| E1 | 1 - sum_c s_c w_c            | theorem harm functional, hand weights, no confinement term |
| E2 | sum_c s_c h_c                | fraction of expression in blood lineages (binary weights, theorem form) |
| E3 | top_share + (1 - sum_c s_c w_c) | current additive rule with the home weight replaced by the expression-weighted weight |
| E4 | top_share + sum_c s_c h_c    | smooth version of the section 23 binary rule |

Baselines re-reported, unchanged: one-factor top_share 0.880; binary two-factor 0.914; old home-weighted continuous 0.809.
Weights, labels, gene list and the cell-type reference are frozen; NO weight will be edited after seeing a result.

## Predictions (written before running)
- P1: E3 and E1 beat the old continuous scheme (0.809), because the home-cell-type artefact is removed.
- P2: E4 lands within 0.03 of binary (0.914) either way; it is the smooth version of the same rule.
- P3: E1 (no confinement term) is weaker than E3 by at least 0.02.
- P4: no hand-weight scheme (E1, E3) beats binary by more than 0.03.
I have no strong prior on whether any scheme reaches 0.914.

## Nulls and uncertainty
- Shuffled-weight null: 1,000 random reassignments of the weight VALUES across the cell types (E1, E3); for E2/E4 the
  haematopoietic flags are reassigned across cell types keeping their count. p = share of shuffles >= observed.
- Percentile bootstrap 95% interval for each scheme's separation, and a paired bootstrap interval for the difference from binary.
- Leave-one-target-out minimum separation for each scheme.
- Coverage: for each scheme with hand weights, the share of expression that falls on cell types with no explicit weight.

## Decision rule (fixed now)
A continuous scheme replaces binary as the PRIMARY rule only if BOTH: (i) separation >= 0.914, and (ii) its
shuffled null p < 0.05. Four schemes are looked at, so any single "win" is described as such (four looks), and a
difference from binary is called real only if the paired interval excludes 0. Otherwise binary stays primary and
the results are reported as a robustness table. Any decrease is reported, not buried.
