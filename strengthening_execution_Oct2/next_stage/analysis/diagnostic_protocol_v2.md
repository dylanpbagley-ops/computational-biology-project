# Diagnostic protocol v2 — source metadata and technical alternatives

Date: October 2, 2026. This version supersedes the hashed v1 draft before any donor-level association or thinning result was computed. The initial v1 file remains preserved for the record. The amendment adds a paired same-cell pre/post comparison and repeated thinning to separate depth-based cell selection from Monte Carlo read sampling.

## Question and scope

For TNFRSF17, GPRC5D, and SLAMF7 in Boiarsky et al. GSE193531, do observable library-depth or cell-quality differences track donor-level threshold tails, and do those tails persist in a common-depth sensitivity after comparing the same retained cells before and after thinning? The analysis can identify measured associations and sensitivity; it cannot assign a causal biological or technical explanation.

## Locked samples, units, and inputs

- Use the existing `run_b_pilot.load()` eligibility rule exactly: `disease_stage == 'NBM' and normal_or_neoplastic == 'normal'`, and `disease_stage == 'MM' and normal_or_neoplastic == 'neoplastic'`.
- Retain all nine NBM donors and all eight MM individuals; no donor is removed for cell count, depth, marker level, or diagnostic result. MGUS and SMM are outside this myeloma comparison. MM-8 has a second SMM-1 specimen from the same person, but that precursor specimen is excluded by the locked disease-stage rule.
- Biological unit is the individual/sample, never a cell or threshold fold. NBM and MM are analyzed as separate strata. No pooled depth-versus-marker association is computed.
- Match cell metadata to extracted markers by exact `index`/`cell_ids` key. Duplicate/conflicting metadata keys, missing keys, or unrecognized selected sample IDs are fatal errors.
- Required fields are `n_counts`, `n_genes`, `frac_mito`, `sample_ID`, `disease_stage`, and `normal_or_neoplastic`. Blank or absent covariates remain `NR`. No imputation, surrogate batch assignment, or causal multivariable model is allowed.
- All protected inputs must match the pre-execution hashes in `input_hashes.json`. Existing threshold, label, result, or protocol files are read-only.

## Primary donor-level summaries

For each eligible donor and marker, report eligible cell count; median and interquartile range of `n_counts`, `n_genes`, and `frac_mito`; raw marker mean, median, detection fraction, and 90th percentile; and normalized marker mean, median, detection fraction, and 90th percentile. Normalized abundance remains `raw UMI / n_counts * 10,000`, as in the existing analysis.

Use the existing frozen 10% pooled-policy normalized threshold for each marker from `analysis/results/b_summary.tsv` to calculate each donor's fraction at or above threshold. Do not refit thresholds. Use the existing raw-scale thresholds for raw-scale tail summaries only. Within NBM and MM separately, calculate descriptive Spearman rank correlations (average ranks for ties) between donor median `n_counts` and donor threshold-tail fraction, and between median `n_counts` and marker detection fraction. Report each donor point and n (9 NBM, 8 MM); no inferential p-values, pooled correlation, regression, or causal language.

## Common-depth sensitivity gate and cell-selection rule

Only run thinning if all checks pass on exact selected cells: marker counts and `n_counts` are finite, nonnegative, integer-valued with `0 <= marker count <= n_counts`; the 29,387-cell source matrix has unique cell columns matching every metadata ID; all matrix counts are integer UMI counts; and per-cell sums across every matrix row equal metadata `n_counts` exactly. If any check fails, record a refusal, test that refusal, and omit thinning.

Set common depth `m` to the floor of the linear-interpolated 10th percentile of the 17 donor-specific median `n_counts` values. Calculate it using depth only, before inspecting any marker values. Within every donor retain the exact cells with `n_counts >= m`; report retained and eligible counts for each donor. If any donor has no retained cell, refuse the full thinning comparison rather than dropping the donor.

## Paired repeated thinning

For each marker and each of 100 independent Monte Carlo replicate streams, draw each retained cell's post-thinning marker count from a hypergeometric distribution with good items equal to its observed marker UMI count, bad items equal to `n_counts - marker UMI`, and draws `m`. Use NumPy's `default_rng(20261003)` and a fixed loop order: replicate 1–100, marker order `TNFRSF17`, `GPRC5D`, `SLAMF7`, sorted sample ID, then original cell order. Marker streams are marginal simulations; they do not represent a measured joint-marker distribution.

For each retained cell, calculate the pre-thinning per-10,000 value as `10,000 * observed marker count / n_counts`, and each post-thinning value as `10,000 * thinned marker count / m`. Compare pre- and post-thinning donor detection fraction, median, 90th percentile, and fraction at or above the frozen normalized 10% pooled threshold on these exact same retained cells. Summarize post-thinning metrics by their mean and 2.5th/97.5th percentiles across the 100 Monte Carlo replicates. These intervals describe simulation variability conditional on the observed cells and counts; they are not patient-level intervals, confidence intervals, or evidence of population transportability.

## Outputs, tests, and claims

Write donor rows, stratum-specific associations, raw-matrix verification, depth-retention counts, paired pre/post summaries, 100-draw simulation summaries, and one depth/tail figure under `next_stage`. Record selected IDs, input hashes, software versions, fixed seed, protocol hash, exact commands, and every failed check in reproduction. Tests cover keyed joins, duplicate/conflicting IDs, unknown samples, preserved unknown covariates, donor grouping, integer/count validation, matrix-total verification/refusal, hypergeometric boundaries, depth selection, and reproducible Monte Carlo output.

Allowed: measured association or lack of association, sample-level variability, observable technical/context associations, and whether a count-valid common-depth sensitivity changes the specific threshold-tail summaries. Prohibited: causal attribution, general population inference, universal RNA or protein selectivity claims, causal correction of donor values, or interpreting missing metadata as a common batch or proof of no effect.
