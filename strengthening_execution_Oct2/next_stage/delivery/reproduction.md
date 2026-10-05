# Reproduction and verification

Run these commands from `C:\Users\dylan.bagley\Desktop\outputs` in PowerShell using the bundled Python 3.12.14 runtime:

```powershell
& 'C:/Users/dylan.bagley/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' -m unittest discover -s strengthening_execution_Oct2/next_stage/analysis/tests -v
& 'C:/Users/dylan.bagley/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' strengthening_execution_Oct2/next_stage/analysis/metadata_diagnostics.py
& 'C:/Users/dylan.bagley/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' strengthening_execution_Oct2/next_stage/analysis/reanalyze_external_full_from_selected_cells.py
& 'C:/Users/dylan.bagley/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' strengthening_execution_Oct2/next_stage/evidence/normalize_clinical_observations.py
& 'C:/Users/dylan.bagley/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' strengthening_execution_Oct2/next_stage/evidence/enrich_clinical_table2.py
& 'C:/Users/dylan.bagley/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' strengthening_execution_Oct2/next_stage/delivery/plot_external_policy_results.py
& 'C:/Users/dylan.bagley/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' strengthening_execution_Oct2/next_stage/delivery/plot_metadata_depth_final.py
& 'C:/Users/dylan.bagley/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' strengthening_execution_Oct2/next_stage/delivery/build_output_hash_manifest.py
```

The clinical commands are sequential. The second reads the hash-locked v3 18-row intermediate and writes the canonical 34-row table. The new next-stage analysis suite was independently rerun by the parent after the root status snapshot test was added: **36 passed**. Parent also reran the clinical stages and confirmed 18 then 34 rows, valid `observation_state`, and final canonical TSV SHA-256 `2ad7972cea12becfda4e8638b896c90cdb96c40e631d37a78807c2903fbf5ad6`. Separately, the original-stage tests under `strengthening_execution_Oct2/analysis/tests` were run as read-only preservation checks: **14 passed**. The original analysis generators were not rerun.

`metadata_diagnostics.py` verifies the 63 protected inputs and starting-state snapshot, the locked protocol hash, source matrix totals, and local diagnostics. The full deterministic reanalysis verifies the annotation workbook, selection/protocol/file manifests, every one of the 456 downloaded count-table hashes and sizes, the first-parse selected-cell-table seal, the exact keyed joins and sorted code order before refreshing the policy outputs. It does not download more data. The preserved first in-flight policy outputs remain under `analysis/results/inflight_policy_outputs_before_deterministic_regeneration/`; they are explicitly superseded, not deleted.

The policy figure script reads only `external_full_newcode_thresholds.tsv` and `external_full_newcode_missingness_bounds.tsv`; the depth figure script reads `metadata_depth_retention.tsv` and `metadata_thinning_paired.tsv`. Both write a value TSV alongside SVG and 300-dpi PNG. Matplotlib 3.11.2 already exists under `strengthening_execution_Oct2/dependencies`, but the default sandboxed Python process sees only an incomplete namespace; figure generation used approved read access to that existing package. No package was installed. The older `metadata_depth_tail.svg` remains preserved; the readable, numeric-axis final plot is `delivery/figures/metadata_depth_thinning_final.svg` and its PNG companion.

Key lock and data hashes:

- v3 expansion protocol: `62eb79d20342c5fa220b75375c3b847108873830486664583ee2ed4fe31a4461`.
- Marker-blind selection lock: `c67c4e8e2fccf1ca0b3a36e61a5064b49cc446b6e0a411b191fe3f845ddb07c0`.
- Actual 456-file count-table manifest: `9275d41c22cc9359aa2c21b1b1b72e3fb2bad080534f229f4c5b146ae214f15f`.
- First-parse selected-cell TSV seal: `1043299440b74ddad9ecd49ad50d249a09ce2ad2564791493f9b029d8af9f1c8`.
- Canonical clinical observations TSV: `2ad7972cea12becfda4e8638b896c90cdb96c40e631d37a78807c2903fbf5ad6`.

Protocol chronology is not described as pristine preregistration. The preserved v1 artifact is `799a246abb8fe939e4937fd5bd9f97eaf04e6a6222b8aee56d720b3544c8478c`; v2 is `aa7b16216786b308bd8fb49f375fb481793fc3169bd3e4719bc9966892a7029b`; the v2 serialization addendum is `0a1e667e83883cfd644be48ac24aee885fb2b4dd808b6d6da5d677711b20b833`; and v3 is the post-pilot, pre-expansion lock above. An earlier execution note's `5f5c961d…f793cd2` lock was not found in preserved files and remains an unverified chronology note.

`delivery/output_hashes.json` is rebuilt last with the command above. It hashes regular next-stage deliverables, including reports, logs, scripts, manifests, and final figures; it excludes this self-referential manifest, downloaded data, bundled dependencies, and caches. The separate 456-file source manifest records raw download integrity.

Known source limitations remain: 2,697 selected annotation IDs have no exact count-table match; one bounded partition check found `W...` count IDs versus workbook `WMC...` IDs and no verified alias; source `Hospital.Code` is not a globally proven person ID; GSE193531 overlap is unknown; and the GSE305971 paper's full Methods/Figure route and runnable z-score code were not accessible for exact method replication. No biological values were imputed for unmatched cells.
