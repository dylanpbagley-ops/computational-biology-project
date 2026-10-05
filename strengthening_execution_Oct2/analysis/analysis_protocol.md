# Selected narrow B analysis protocol

October 2, 2026, after pilot results were inspected. Post-pilot development specification; it is not an untouched test. Primary specification in feasibility/candidate_b/pilot_specification.md is preserved with its prior hash. No endpoint or target chosen to rescue historical AUC.

Retain exactly the pilot data, groups, three genes, normalization, deterministic threshold class, budget .10, policies and held-out split. Extend reporting with donor-specific held-out summaries, 2000 donor bootstraps of equal-donor threshold coverage and normal mismatch (seed 20261002), and secondary raw-count analysis at all four predeclared budgets. Bootstrap independently resamples the 8 MM and 9 NBM donors; never resample cells as biological subjects. Percentile intervals describe this small empirical donor distribution, not universal coverage or a rigorous uncertainty bound. Do not bootstrap the 72 dependent folds as independent subjects.

Primary interpretation: donor heterogeneity limits a shared RNA threshold on these measured plasma cells. Report whether constraints observed in training transfer to withheld donors. Negative outcomes/poor malignant coverage are preserved. No universal guarantee on unseen donors, no RNA-to-protein conversion, no efficacy/harm interpretation. Figures link to exact gene/policy/budget/output versions. Depth/raw-count comparisons are diagnostics of measurement dependence, not a second mechanistic uncertainty model.

External check feasibility: inspect local myeloma ADT array names and original source gating. It is independently measured protein but cannot independently validate malignant versus healthy-plasma-cell RNA discrimination; source gating itself uses target markers. Report this mismatch rather than forcing a replication. No new inaccessible data search.

Original analysis: reproduce pairwise AUC on literal TARGETS and unchanged top_share + M2_flag. Compare specificity, tau and flag on identical original rows. This is historical reproduction only. No amended labels are created while clinical gate remains closed; no new clinical p-values or selective subsets.

Tests: concrete threshold examples/exhaustive finite cases, input validation, donor weights, source/observational missingness, score-free initial fields, held-out disjointness, protected input hashes. Run every new standalone entry point. Preserve protocol/source hashes, seed, runtime and outputs in delivery/reproduction.md.
