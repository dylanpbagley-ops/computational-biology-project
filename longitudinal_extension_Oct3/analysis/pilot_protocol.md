# Locked descriptive pilot protocol

Protocol lock must be timestamped in pilot_protocol.lock.json before downloading or parsing any expression matrix.

## Question and claim boundary

For three documented t(11;14) participants, describe the change in same-cell transcript-detection distributions for the fixed target panel TNFRSF17, GPRC5D, SLAMF7 in paired CD138-selected pre/post-venetoclax bone-marrow specimens. This is a descriptive specimen-level RNA question. It does not establish malignant identity for each barcode, antigen protein availability, functional escape, a causal venetoclax effect, tracked cell lineage, response prediction, or generalization beyond these three participants.

## Locked specimens

Use exactly these pairs, all derived from the RNA columns in the primary Blood Table 1 and cross-checked against GEO sample records:

- PM01: GSM5971019 / P1393 pre → GSM5971023 / P1555 post.
- PM07: GSM5971021 / P1389 pre → GSM5971025 / P1479 post.
- PM10: GSM5971020 / P1217 pre → GSM5971024 / P1438 post.

Do not use PM09 (post RNA sample P1698 is not in the checked GEO sample records) or substitute PM01's P1426 DNA specimen for its P1393 RNA specimen. Do not add GSE226335/GSE301286 records. The six sample files were selected by source identity and order, not by target values. Keep the chemistry metadata conflict (sample-level v2 versus series-level v3) visible in every interpretation.

## Input and cell scope

Use the public filtered_feature_bc_matrix.h5 files for the six GSM records only. Before counting targets, verify HDF5 structure, unique barcodes, unique feature identifiers, exact gene symbols, feature types, integer nonnegative raw UMI values, matrix dimensions, and sample-level file correspondence. Sum total UMI depth over Gene Expression features only; exclude antibody/capture features from the depth denominator. If feature types are absent, use all features only when metadata confirms a Gene Expression-only filtered matrix. If any target feature is absent, duplicated ambiguously, noninteger, or not represented as raw counts, stop before target outcomes are calculated. A sparse missing coordinate is a zero only when that measured feature is present in the same matrix; missing features/metadata are unavailable, never zero. Preserve the unresolved sample-record v2 versus series-design v3 metadata discrepancy.

Include every barcode in each deposited filtered matrix. No additional QC threshold, state filter, or after-result cell exclusion is allowed. The GEO source is primary CD138+ MM cells, but no per-barcode malignant/clone labels are verified. Report all cells as CD138-selected filtered-matrix barcodes and no inferred malignant-state results. Each specimen is a separate biopsy; barcodes are never matched across specimens.

## Primary endpoints

Use strictly positive observed integer counts for detection. For each specimen and each of the three unordered pairs (TNFRSF17/GPRC5D, TNFRSF17/SLAMF7, GPRC5D/SLAMF7), report the four joint cell counts:

- both detected (n11);
- first detected only (n10);
- second detected only (n01);
- neither detected (n00).

Also report per-marker marginal detected/total fractions, OR coverage ((n11+n10+n01)/N), double-undetected fraction (n00/N), and both directed conditional complementarity fractions (B detected among A-undetected and A detected among B-undetected). A zero conditional denominator is undefined/NA, not zero. Report each specimen separately, plus within-participant post-minus-pre changes. Participants are the independent units; cells are subsamples. Do not calculate cell-as-person p-values or claim a population-level treatment effect.

No pair-selection algorithm, ranking, prediction baseline, or method-superiority analysis is part of this pilot. The target family and all pair comparisons are fixed before expression retrieval; expected OR coverage is only a descriptive quantity.

## Depth sensitivity

Primary detection uses the deposited raw UMI counts with no normalization threshold beyond count > 0. Report each specimen's cell count and per-barcode total UMI distribution. As a sensitivity check, compute each specimen's median total UMI over Gene Expression features and set integer common pseudo-depth D = floor(minimum of the six specimen medians). For each barcode, draw min(D, observed_total_UMI) molecules without replacement from that barcode's complete RNA count vector using a multivariate hypergeometric draw over (target A, target B, target C, all other Gene Expression counts combined). This preserves the joint target distribution during thinning. Retain every barcode; report the fraction below D and 100 deterministic sensitivity replicates using seeds 20261003 through 20261102. This is a depth sensitivity only, not causal correction.

If D < 1 or NumPy lacks a reliable multivariate hypergeometric primitive, do not thin each marker independently; mark depth sensitivity unavailable and report unresolved depth confounding.

## Missingness, exposure, and decision criteria

List missing genes, specimens, barcodes, and any failed joins. Do not convert absent modalities, failed annotations, or unmeasured features to zero. Already-viewed prior-stage GSE305971 expression remains development data and is not a validation cohort. The GSE199359 matrix is newly exposed once downloaded.

A result is descriptive only. Consistent pair-level distributions across participants would support a narrow feasibility observation, not a method or biological validation. Inconsistent directions, very small cell denominators, matrix/comparability failures, or changes erased by joint depth thinning count against an interpretable extension. No threshold, target, population, or endpoint may be changed after seeing results. No full analysis follows from this three-participant pilot.

