# Reproduction and verification

Root: C:/Users/dylan.bagley/Desktop/outputs. New outputs isolated in strengthening_execution_Oct2. input_hashes.json protects 63 original label/prediction/score/preceding-pass/input files. No Git repository/commits/worktree.

Runtime: bundled Python 3.12.14, NumPy 2.3.5. SciPy absent despite plan's older runtime record; no optimizer used. Matplotlib 3.11.2 and dependencies installed into execution/dependencies only. Initial network attempt failed; elevated installation succeeded. Sandbox cannot read those installed directories, so elevated figure execution was necessary and succeeded. Do not install into legacy package folders.

Run from ROOT:

```powershell
& 'C:/Users/dylan.bagley/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' strengthening_execution_Oct2/evidence/build_batch.py
& 'C:/Users/dylan.bagley/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' strengthening_execution_Oct2/analysis/run_b_pilot.py
& 'C:/Users/dylan.bagley/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' strengthening_execution_Oct2/analysis/run_extension.py
& 'C:/Users/dylan.bagley/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' -m unittest discover -s strengthening_execution_Oct2/analysis/tests
& 'C:/Users/dylan.bagley/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' strengthening_execution_Oct2/analysis/make_figures.py
```

startup.py successfully ran once and refuses to replace initial hashes; do not reset them. Other scripts verify protected input hashes before/after and replace only new execution outputs. Protocol lock rejects changes: use an explicitly versioned new analysis for departures. Figure command on this host requires access to installed packages.

Frozen records: evidence/endpoint_v03.sha256; feasibility/candidate_b/pilot_specification.sha256 (430812ec57a6c738aac1672f656a160b9078be252b3ca232a1c170681bf2b9df); analysis/frozen_specification.sha256 (10e5a8098af0228b067e45ca14d1a653196c29c8daee4e599a33c90262694400). Pilot reuses previously inspected data; full protocol explicitly post-pilot. Neither timestamp creates pristine validation. Seed 20261002,2000 donor draws per gene; RNG stream consumed sequentially across genes. No clinical fitting/tuning.

Fresh verification:14 tests passed. Seven numerical tests use independent expected values/exhaustive finite examples; initial import failure preceded implementation, then tests passed. Seven evidence/output checks cover44-row screening, missing versus zero, blank reviewer fields, trial overlap, donor split exclusion, observed-envelope constraints and protected hashes. Batch/pilot/extension/figures standalone entry points passed after loader change. Direct pairwise historical AUC reproduced 19/11 and 0.913875598086; direction/ties tested. All63 protected inputs unchanged. No unrelated legacy whole-suite reruns.

Numerical record: analysis/results/*.tsv/JSON. b_summary contains exact empirical threshold enumeration; b_held_out and b_held_out_donors preserve dependencies; b_bootstrap is small-donor sampling diagnostics, not guaranteed bounds; b_raw_sensitivity is measurement-scale diagnostic. extension_results records runtime/protocol/seed; pilot_run.log and extension_run.log retain final output. output_hashes.json records new sources/reports/results/figures (dependencies/cache/itself excluded).

Figure traces: b_observed_donor_tradeoff uses b_summary at budget.10 and matched donor cell values; b_held_out_donors uses corresponding donor-level table (one dot per individual, averaged over partner folds); historical_baselines uses historical_reproduction.json; diagnostic_a_source_types uses candidate_a/comparability.tsv (workload, not outcomes). Four300dpi PNG/editable SVG pairs generated; main donor plots visually inspected.

Self-review checked finite-threshold optimality, common-response/oracle scope, denominators/missingness, disjoint sample IDs, protocol timing and claim traceability. It is author self-review, weaker than independent review; no agents spawned. No manuscript/submission prose written. Original outline retained; new advisory outline supplies proposed structure.

Remaining limits: full clinical priority extraction incomplete; Thomas full methods not retrieved; several direct PMC pages blocked while indexed primary sections were available; no reviewer identified; no comparable independent RNA replication; prediction timestamp audit not done. Raw-scale changes and withheld-donor failures retained. Passing software checks do not supply source adjudication or scientific novelty.
