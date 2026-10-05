# Diagnostic protocol — source metadata and technical alternatives

Date: October 2, 2026. This protocol is frozen before running the new donor-level diagnostics. It is exploratory and will not amend any biological identity labels or clinical outcome labels.

## Question and scope

For the existing TNFRSF17, GPRC5D, and SLAMF7 common-RNA-threshold result in Boiarsky et al. GSE193531, do observable library-depth or cell-quality differences track donor-level threshold tails, and do the main donor-tail patterns persist after a valid common-depth sensitivity? This analysis can identify measured associations and sensitivity. It cannot assign a causal biological or technical explanation.

## Locked samples, units, and inputs

- Use the existing `run_b_pilot.load()` eligibility rule exactly: `disease_stage == 'NBM' and normal_or_neoplastic == 'normal'`, and `disease_stage == 'MM' and normal_or_neoplastic == 'neoplastic'`.
- Retain all nine NBM donors and all eight MM individuals; no donor is removed for cell count, depth, marker level, or diagnostic result. MGUS and SMM are outside this myeloma comparison. The MM-8 individual also has an SMM-1 specimen, which is not in the selected comparison.
- Biological unit is the individual/sample, never a cell or threshold fold. NBM and MM are analyzed as separate strata. No pooled depth-versus-marker association is computed.
- Match cell metadata to extracted markers by the exact `index`/`cell_ids` key. Duplicate/conflicting metadata keys, missing keys, or unrecognized selected sample IDs are fatal errors.
- Required fields are `n_counts`, `n_genes`, `frac_mito`, `sample_ID`, `disease_stage`, and `normal_or_neoplastic`. Blank or absent covariates remain unknown. No imputation, surrogate batch assignment, or causal multivariable model is allowed.
- All protected inputs must match `input_hashes.json`. Existing threshold, label, result, or protocol files are read-only.

## Primary donor-level summaries

For each eligible donor and each of the three markers, report eligible cell count; median and interquartile range of `n_counts`, `n_genes`, and `frac_mito`; raw marker mean, median, detection fraction, and 90th percentile; and normalized marker mean, median, detection fraction, and 90th percentile. Normalized marker abundance remains `raw UMI / n_counts * 10,000`, as used by the existing analysis.

Use the existing frozen 10% pooled-policy normalized threshold for each marker, read from `analysis/results/b_summary.tsv`, to calculate the per-donor fraction at or above that threshold. Do not refit thresholds. Within each disease stratum separately, calculate descriptive Spearman rank correlations (average ranks for ties) between donor median `n_counts` and donor threshold-tail fraction, and between median `n_counts` and marker detection fraction. Report donor points and n (9 NBM or 8 MM); no inferential p-values, pooled correlation, regression, or causal language. The corresponding existing raw-scale thresholds are used for raw-scale tail summaries only.

## Common-depth sensitivity gate and algorithm

Only run thinning if the following provenance checks pass on the exact selected cells: every target UMI value and `n_counts` is finite, nonnegative/integer-valued, and `0 <= target UMI <= n_counts`; the 29,387-cell raw matrix header matches the metadata cell IDs without duplicates; the matrix contains integer UMI counts; and the per-cell sum over all matrix rows equals metadata `n_counts` exactly. If any check fails, create a refusal record, test that refusal condition, and omit thinning.

If the gate passes, set common depth `m` to the floor of the linear-interpolated 10th percentile of the 17 donor-specific median `n_counts` values, calculated from all selected cells before inspecting marker outcomes. No marker values determine `m`. In every donor, retain cells with total `n_counts >= m`; report retained/eligible cell counts for every donor. If any donor would have zero retained cells, refuse the entire thinning analysis rather than drop that donor.

For each marker and retained cell, draw the post-thinning marker count from the exact hypergeometric distribution with good items equal to that marker's observed UMI count, bad items equal to `n_counts - marker UMI`, and draws `m`. Use one NumPy `Generator` with seed `20261003`, in fixed order `TNFRSF17`, `GPRC5D`, `SLAMF7`, then sorted sample ID and original cell order. Convert thinned marker counts to per-10,000 abundance using `10,000 * thinned_count / m`. Report each donor's thinned detection fraction, median and 90th percentile, and its fraction above the pre-existing normalized 10% pooled threshold. Do not interpret this as calibrated full-library dropout or protein absence.

## Outputs, tests, and claims

Write donor rows, stratum-specific association rows, raw-matrix verification, thinning decision/results, and a single depth/tail figure under this `next_stage` directory. Include fixed seed, selected IDs, input hashes, software versions, exact commands, and all failed checks in reproduction. Tests must cover keyed joins, duplicate/conflicting IDs, unknown samples, preserved missing values, donor grouping, integer/count validation, hypergeometric boundaries and deterministic output, and the refusal path when totals are not verified.

Allow: measured associations or their absence, sample-level variability, observable technical/context associations, and whether a count-valid common-depth sensitivity changes the specific threshold-tail summaries. Prohibit: causal attribution, general population inference, universal RNA or protein selectivity claims, causal correction of donor values, or describing an absent metadata field as a shared batch or as proof of no effect.
