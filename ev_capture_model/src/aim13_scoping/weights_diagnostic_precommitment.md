# Round 6A precommitment: are the expendability weights the bottleneck? (written before running)

Date: 2026-09-29. Inputs frozen: target_labels_frozen.tsv, wholebody_v2.json, lineage_map.py, flag_cleanup.py.

Question (Dylan): "the expendability weights might be messing things up." Test: how much could ANY weights on the
10 lineage families add to the frozen rule (top_share + blood flag, separation 0.914 on 19 safe v 11 toxic)?

Model class F (family offsets): score(g) = top_share(g) + a[k(g)], k(g) = the lineage family with the largest
aggregate share of gene g's expression (as in flag_cleanup.flags, blood = frozen haematopoietic list),
a_k in [0, 1.5], 10 free numbers. The frozen rule is the point a_blood = 1, all others 0.
Model class G (family composition): score(g) = top_share(g) + sum_k a_k * share_k(g), a_k in [0, 1.5].

Procedure:
 1. Oracle (in-sample): maximise separation (AUC) over a. This is an UPPER BOUND on what any family-level weighting
    can do on this panel, and it is optimistic (10 numbers fitted to 30 labels).
 2. Null for the oracle: the same optimisation on 200 label permutations, to see how large an in-sample gain
    the optimiser finds from noise alone.
 3. Honest estimate: nested leave-one-target-out. Fit a on 29 targets, score the held-out one, pool the 30
    held-out scores, compute separation.
 4. Stratified view: separation among non-blood-largest targets only (8 safe v 11 toxic, 88 pairs), where all 16
    misordered pairs of the frozen rule live.

Interpretation rules (fixed now):
 * If oracle gain over 0.914 is < 0.03, family-level weights cannot materially help: weights are not the bottleneck.
 * If nested-CV separation of learned weights <= 0.914, learned weights do not beat the binary out of sample.
 * Learned weights are diagnostic only. Nothing is adopted from this round.
