# Mathematical scope correction — advisory technical note

October 2, 2026. This note supersedes the overstatements in the original section 29. It is an explanation and verification record, not student submission prose. The two distributions and every weight must be specified; normal expression alone does not determine a tumor-versus-normal frontier.

Let p_T and p_c be probability distributions on a common ordered measurement grid, q_c nonnegative compartment prevalences summing to one, and w_c fixed nonnegative harm weights. All compartments share the same response f(x) in [0,1]. Define T = E_T[f], H = sum_c q_c w_c E_c[f], and W = sum_c q_c w_c.

If W>0, define p_W = sum_c q_c w_c p_c / W. Then H/W = E_W[f]. This identity is exact under these assumptions. It is not a measured clinical severity model, and the weights are not inferred from the expression score. If W=0, H=0 for every response and the normalized mixture is undefined; f=1 attains T=1 without modeled harm.

## Two different response classes

For **nondecreasing f**, the feasible set in (H/W,T) coordinates is the convex hull of the antigen-threshold ROC points, including f=0 and f=1. On n ordered bins, write f as a nonnegative weighted sum of upper-tail indicator functions, with the remaining weight on f=0. Conversely every such convex combination is a nondecreasing response. The upper boundary is the best T at the specified budget within this response class.

For **arbitrary f**, the feasible set is the linear image of the cube [0,1]^n: the convex hull of all deterministic bin-subset responses. Its upper boundary is obtained by likelihood-ratio ordering, with randomization at a boundary if needed. It is not in general the convex hull of only the upper likelihood-ratio ROC curve; arbitrary responses can also occupy the lower part of the full feasible set. Bins with positive tumor mass and zero weighted-normal mass are handled without dividing by zero. This is classical Neyman-Pearson/fractional-knapsack geometry, not a new statistical principle.

The arbitrary-response upper boundary also bounds the monotone class, but need not be attainable within it. Equality of the upper boundaries follows when likelihood-ratio ordering is compatible with the antigen order, with appropriate treatment of ties. Physical realizability of every mathematical response is not assumed or demonstrated.

Verified example: p_T=(0.45,0.10,0.45), p_W=(0.10,0.80,0.10), budget 0.20. Arbitrary responses attain 0.90; nondecreasing responses attain 23/45, approximately 0.5111. The full ROC shape matters; AUC alone is insufficient, including for likelihood-ratio ROC curves (see audit_followup_Oct2/recheck_results.json).

## Multiple tissue constraints

Selecting a single nonzero weight recovers that compartment's reference. However, the minimum of separately optimized compartment ceilings is only an upper bound on the optimum under all compartment budgets simultaneously. Each separate problem can have a different optimizing response.

Counterexample: p_T=(0.5,0.5), normal references (1,0) and (0,1), with a budget of 0.5 in each. Each separate arbitrary-response optimum is 0.75, but simultaneous constraints require f_1<=0.5 and f_2<=0.5, yielding at most 0.5, attained at f=(0.5,0.5). Therefore one-hot identities or AUC comparisons do not verify joint attainability.

## Relationship to the empirical score

The identity sum w_c s_c = 1 - sum (1-w_c)s_c for expression shares s_c summing to one is algebraically correct. It does not establish H proportional to that quantity. A sufficient additional assumption would be q_c E_c[f] = a s_c with the same constant a across compartments; that assumption is not established by the project's measurements. The identity does not imply that concentration in fewer types minimizes harm or derive the blood bonus.

For fixed n and nonzero total expression, top_share = max(x)/sum(x). The standard tau specificity index equals (n - 1/top_share)/(n-1), so it has exactly the same ranking. This is already recognized in the project's section 37, and must be central to novelty framing. All-zero profiles require a no-data policy, not an automatic safety verdict. The breadth lemma concerns expression exceeding a relative cutoff under fixed partitioning; it does not infer cell prevalence, protein abundance, tissue accessibility, or clinical injury.

## Verification

`frontiers.py` implements the two finite-grid upper boundaries separately. `test_frontiers.py` includes the counterexample, zero-mass handling, input validation, mixtures, tau, random monotone-response containment, and independent enumeration of arbitrary-response linear-program vertices. `frontier_results.json` records both counterexamples. The legacy weighted verifier tests arbitrary-response upper bounds and likelihood-ratio utility, not the full monotone theorem or clinical relevance. Its old special-case output names are retained for historical compatibility and must be read with this correction.

Sources: [likelihood-ratio ROC theory](https://link.springer.com/article/10.1007/s13171-021-00244-5); [specificity metric benchmark and tau formula](https://pmc.ncbi.nlm.nih.gov/articles/PMC5444245/). The counterexamples and finite-grid proofs above are directly checkable mathematical arguments, not empirical clinical claims.
