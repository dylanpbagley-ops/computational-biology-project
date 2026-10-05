# Round 4B: what else could make the empirical rule more effective? (precommitment, Sept 29, 2026, BEFORE running)

## Diagnosis that motivates these ideas (from reasoning about the method, not from looking at scores)
1. **Shares hide absolute abundance.** top_share is a fraction of a gene's total expression, so a gene that is barely
   expressed anywhere (DLL3 is "essentially absent from adult normal tissue", per the label table) can look as diffuse as
   a gene that is highly expressed everywhere. A therapy responds to actual antigen amounts, not fractions.
2. **top_share depends on how finely cell types are split.** Splitting one lineage into sub-types lowers the top share.
3. **top_share uses one number (the maximum) to summarise a 154-value profile.** Entropy, tau and top-3 use more of it.
4. **The blood flag is the only expendability information**, and the panel has only 8 blood-flagged safe targets.

## Score family (all fixed now; ALL reported; higher = predicted safer)
x = the gene's 154 cell-type nCPM values, s = shares, flag = the round 4A plurality-lineage blood flag (M2).
| id | measure | id | measure |
|---|---|---|---|
| C0 | top_share (baseline)                     | C5 | -n_types (types above 10% of top expression) |
| C1 | top lineage share (`lineage_map.py`)     | C6 | -log10(1 + sum_c x_c)  (total normal expression) |
| C2 | 1 - entropy(s)/ln(154)                   | C7 | -log10(1 + max_c x_c)  (highest single-cell-type expression) |
| C3 | tau (Yanai 2005 tissue-specificity index) | C8 | mean of the panel-rank of C0 and of C6 (confinement + low abundance) |
| C4 | top-3 share                              |    |  |
Each is evaluated (a) alone and (b) as `1000*flag + C` (flag first, then C), giving 18 schemes with the baselines
C0 (0.880) and C0 + flag (0.914) among them, so 16 are new looks. (C8 uses ranks within the 30 panel targets.)

## Predictions (written before running)
- Q1: C6 and C7 (abundance) alone score at least 0.70, because absent-from-adult-tissue targets (DLL3, GPC3) are safe.
- Q2: C2, C3, C4 land within 0.05 of C0 (they summarise the same profile).
- Q3: C1 (coarser, lineage level) is no better than C0.
- Q4: at least one scheme beats 0.914 in raw separation (I put this near 35%; with 16 looks on 30 targets a raw win is
  weak evidence on its own).

## Multiplicity control (fixed now)
A family-wise max-statistic permutation test: shuffle safe/toxic labels 2,000 times, recompute all 18 schemes each time, and
record the MAX separation. Each scheme's family-wise p = share of shuffles whose max >= its observed separation.
This charges every scheme for the fact that the best of 18 was picked. Also: label-permutation p, bootstrap CI,
leave-one-out minimum, and the paired difference vs C0 + flag.

## Decision rule (fixed now)
No scheme replaces the current rule on the strength of these 30 targets alone. A scheme is only PROMOTED to "hypothesis to
confirm on new targets" if separation >= 0.944 (a gain of at least 0.03 over 0.914) AND family-wise p < 0.05. Promoted
schemes enter the paper only after they are confirmed on labelled targets that were not used here. Everything else is
reported as a robustness table. Cumulative alternatives tried after this round: 30 (14 + 16).
