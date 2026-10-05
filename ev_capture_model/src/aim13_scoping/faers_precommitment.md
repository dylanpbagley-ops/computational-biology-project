# Round 6B precommitment: measured organ severity weights and a multi-mechanism "most realistic" rule (written before any count was fetched)

Date: 2026-09-29. Frozen before running: faers_severity.py (28 systems, PT lists, 147 cell types mapped + 7 unmapped, IHC tissue mapping),
target_labels_frozen.tsv, wholebody_v2.json, expendability_weights.py, lineage_map.py, flag_cleanup.py. Hashes in faers_precommitment.md.hashes.

## What is being built
1. MEASURED WEIGHTS. For each of the 28 organ systems: severity = (# FAERS reports recording death among reports that list any PT of that system)
   / (# reports listing any PT of that system), from openFDA drug/event counts (all drugs, all years to the fetch date).
   Weight w_sys = severity / max severity. Unmapped cell types/tissues get the MEDIAN weight; sensitivity runs use MIN and MAX.
   The weights are computed with no reference to any target score or label.
2. RULES (7, all scored the same way; higher score = predicted safer):
   A1  top_share + (1 - w_meas(home cell type))
   A2  top_share + (1 - sum_c s_c * w_meas(c))                      s_c = share of the gene's expression in cell type c
   A3  -H_RNA,  H_RNA = max_sys w_sys * E_sys^RNA                  WORST-ORGAN aggregation, absolute expression, thresholds integrated
   A4  -sum_sys w_sys * E_sys^RNA                                   same as A3 with SUM aggregation (tests the aggregation choice)
   A5  A3 with the existing per-cell-type hand weights instead of measured system weights (control; unlisted cell type = 0.70)
   A6  -H_IHC,  H_IHC = max_sys w_sys * E_sys^IHC                   PROTEIN evidence, organ level
   A7  -max_sys w_sys * (E_sys^RNA + E_sys^IHC)/2                   RNA and protein together
   Exposure. RNA: e_c(g) = mean over t in {1,3,10,30,100} of 1[nCPM_c(g) >= t]; E_sys^RNA = max over the system's cell types (unmapped types form one
   extra system with the median weight). IHC: level Not detected/Low/Medium/High = 0/1/3, 2/3, 1; rows with reliability 'Uncertain' or level not in that set
   are dropped; E_sys^IHC = max over rows mapped to the system (unmapped tissues form one extra system with the median weight).
   Genes with no IHC rows are excluded from A6/A7 only (n reported).
   Robustness only (not counted as a new rule): A3 with rank weights (rank of severity / 28).

## Pre-declared realism ranking (fixed now, independent of any score)
Mechanisms modelled: (i) absolute level not share, (ii) threshold uncertainty integrated, (iii) worst-organ aggregation, (iv) measured severity,
(v) protein evidence. A7 has all five, A3 and A6 have four, A4 three, A5 three (hand severity), A2/A1 one to two, the frozen rule none of them
(it is share-based, uses a binary lineage flag). The rule I recommend as "most realistic" is therefore A7, and it is recommended REGARDLESS of its score.

## Outputs and decision rules
Report for every rule: separation on 19 safe v 11 toxic (exact/label-permutation p), paired bootstrap difference vs the frozen rule (0.914), leave-one-out
minimum, separation among the 22 targets not flagged blood (11 v 11), validation-set separation (4 safe v 2 toxic, scorable only), and the weight-shuffled null
(permute the 28 system weights, 1000 times; for A5 permute cell-type hand weights): does the measured weighting beat arbitrary weights of the same structure?
Family-wise null = maximum over the 7 rules of each shuffled run.
Weights themselves: report the table, Spearman correlation of measured severity with the mean hand weight of each system's cell types, and flag any
result that looks absurd (e.g. a clearly mild reaction above a clearly severe one).
 * PROMOTE a rule over the frozen rule only if separation >= 0.944 AND family-wise p < 0.05; then it would still need confirmation on new targets.
 * Otherwise the frozen rule stays primary; the most realistic rule is reported next to it with its score and paired CI, and described as statistically
   tied if the CI includes 0.
 * A rule that scores lower is NOT called worse (see the noise simulation); a rule that scores higher on 30 targets is NOT called better without the promotion test.
Cumulative count of alternatives tried against the frozen rule after this round: 52 + 7 = 59 (57 distinct after tau=top_share and A5~E-type overlap is not deducted; count is conservative).
