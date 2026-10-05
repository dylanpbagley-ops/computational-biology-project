# Continuous expendability, third attempt: EXTERNAL weights (precommitment, written Sept 29, 2026, BEFORE any scheme below was run)

## What is new
Rounds 1 and 2 used weights I assigned by hand. This round uses weights from a published measurement made by other
people for another purpose: the lifespan of each human cell type (Sender & Milo, Nat Med 2021; data file
`ev_capture_data/cell_turnover_sender_milo/Summary.xlsx`). Idea: fast-turnover cell types are replaced by the body,
so their loss is expendable; long-lived, non-renewing ones are not.
`w_c = log10(lifespan in days) / log10(300000)`, clipped to [0,1] (neurons = 1). Module: `external_turnover_weights.py`.

## Fixed in advance
- The mapping of HPA cell types to the 24 lifespan categories (family/lineage rule) is in `external_turnover_weights.py`,
  written from cell biology and NOT edited after any score is computed. 77 of the 154 HPA types map; the rest are
  UNMAPPED and get the median weight of the mapped types (0.391). Where the file has no human lifespan the file's
  rodent value is used (flagged R in the module).
- Honest disclosure: I knew the earlier results and the labels while writing the mapping. It was written from
  lineage only, but it was not blinded in the strict sense.
- Score direction is fixed: higher = predicted safer. A scheme that works "backwards" does not count.

| id | score | what it is |
|---|---|---|
| X1 | 1 - sum_c s_c w_c                     | harm functional, turnover weights, no confinement term |
| X2 | 1 - w_home                            | home cell type only |
| X3 | top_share + (1 - sum_c s_c w_c)       | confinement + expression-weighted turnover harm |
| X4 | top_share + (1 - w_home)              | confinement + home-cell turnover |
| X5 | top_share + 1[home lifespan <= 30 d]  | binary "fast-renewing" flag, same form as the blood flag |
| X6 | top_share + 1[home lifespan <= 365 d] | same, wider cut |

Baselines re-reported: one-factor 0.880; binary two-factor 0.914 (blood flag).

## Predictions (written before running)
- P1: no X scheme reaches 0.914. Lifespan measures how fast a cell is replaced, not whether losing it for the days it
  takes is survivable (colon epithelium turns over in 3.4 d, yet CEACAM5 causes colitis).
- P2: X3 and X4 do not exceed the confinement-only 0.880 by more than 0.02.
- P3: X5 is the weakest of the binary flags, because the T and B cells that make the blood flag work have lifespans of
  63 to 320 days.
- P4: none of the six passes the shuffled-weight null at p < 0.05.
I would be surprised (probability under 15%) if any scheme were adopted.

## Nulls, uncertainty, coverage
Same as round 2: 1,000 shuffles of the weight values (X1-X4) or the flag values (X5, X6) across the 154 cell types;
5,000 bootstrap draws for the interval and for the paired difference from binary; leave-one-target-out; coverage =
mean share of expression on mapped types, and the number of the 30 targets whose home type is mapped.

## Decision rule (unchanged)
Replace binary as the PRIMARY rule only if separation >= 0.914 AND shuffled-null p < 0.05. Six schemes are looked at,
so any single "win" is described as one of six looks (and this is the third round: three rounds; 12 non-binary schemes in
total, counting the first round). Otherwise binary stays primary and the results go in the robustness table.
