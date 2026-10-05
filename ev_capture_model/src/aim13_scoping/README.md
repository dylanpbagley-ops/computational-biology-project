# aim13_scoping — scripts behind docs/Aim13_Therapeutic_Index_Ceiling_PROPOSAL.md

Exploratory / scoping scripts from the September 6, 2026 session. They are the
evidence for every number in that proposal document. They are NOT yet
production aims: no tests, crude percentile gating, no integration with the
existing `src/` module conventions. Promote them properly (module docstring,
tests, shared helpers) if the proposal is adopted.

Run order:
1. `verify_cdf_identity.py` — proves escape_fraction == F(n*) against real data.
2. `cohort_gating_and_batch_scan.py` — writes cohort_scan_results.json (41 files).
3. `confound_quantification.py` — the gating/batch confound numbers.
4. `pilot_tumor_vs_bystander_auc.py` — writes pilot_roc_results.json (41 files x 37 antigens).
5. `rank_antigens.py` — the two rankings and their Spearman correlation.
6. `ground_truth_check.py` — validated-vs-broad separation test.
7. `gating_robustness.py` — 40/50/60/70th-percentile gate sweep.

All read from `ev_capture_data/raw_cytof/` and import the FCS parser from
`aim12_single_cell_cytof_reordering.py`. Nothing in `src/` was modified.

Note: the scripts write intermediate JSON to ~/scratch by default; paths are
inline at the top of each file.
