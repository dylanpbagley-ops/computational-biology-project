# Adjudicating "more realistic" rules without leaning on 30 safe/toxic labels (precommitment, Sept 29, 2026, BEFORE running)

## Why
A truer rule can lose on a small panel to noise (simulation: `rule_comparison_noise_simulation.py`). So realism should be judged by evidence that is
different from, and richer than, one safe/toxic bit per target. Test 1 asks whether a model of harm puts the harm where it actually fell.

## Test 1: tissue location ("naming")
For each of the 13 toxic targets (11 panel + GUCY2C and PSCA from the validation set) the damaged tissue is fixed in `damage_tissue_map.py`.
A model assigns each of the 154 cell types a harm contribution. Score = the best (smallest) rank of any documented-damaged cell type
(1 = the model puts the most harm exactly there); hit@3 = best rank <= 3.
| id | harm contribution of cell type c |
|---|---|
| N1 | s_c (share of the gene's expression) |
| N2 | s_c x 1[c is vital]  (vital = not blood, blood_other, germline, placental) |
| N3 | s_c x hand weight (expendability_weights.py, default 0.70) |
| N4 | s_c x lifespan weight (external_turnover_weights.py, unmapped = median) |
Predictions (before running): N2 improves on N1 for CA9 and EGFR (germline/placental removed) and worsens CD47 (its damaged tissue, erythrocytes, is
"expendable"); N4 is the worst (colonocytes get weight 0.10, so CEACAM5's harm shrinks); N3 within one target of N1 either way.
Comparisons vs N1: exact sign test on per-target best-rank changes (ties dropped).

## Test 2: candidate "truer" rule built from the vital-share idea (2 more looks; cumulative 46)
R1 = 1 - (largest share of the gene's expression in any VITAL cell type);  R2 = R1 + 1[M2 blood flag].  R1 says how much of the gene's expression sits in
the single most exposed vital cell type; it needs no separate flag because expendable cell types are simply not counted.
Decision rule for treating a rule as "truer" and adopting it (ALL must hold): (a) N2 beats N1 on naming with sign-test p < 0.05; (b) R1 or R2 has separation >= 0.88
on the 30-target panel with a paired-bootstrap interval against the frozen rule that includes 0; (c) on the validation primary set it is at least as good as the frozen rule.
Otherwise the frozen rule stays and the result is reported as information.
