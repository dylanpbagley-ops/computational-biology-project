# Contribution gate decision: stop after the descriptive pilot

## Decision

The paired-data gate passed for a three-participant descriptive pilot, but the strong-contribution gate failed. Do not begin a full analysis of this cohort, make a prediction model, or claim a substantial extension from these results. Retain the reproducible pilot as feasibility evidence and a negative contribution-gate result.

G1 passed only narrowly: the accessible main text for the closest sources did not report this exact fixed-panel joint-detection contrast, but the broader biology of single-cell target co-expression and longitudinal target changes is established. The relevant supplements were not fully verified, so novelty remains conditional. G2 passed for sample identity/order and file access: primary Blood Table 1 links PM01 P1393→P1555, PM07 P1389→P1479, and PM10 P1217→P1438 in pre/post-venetoclax RNA columns; all six GEO files were downloaded after lock. This supports a small pilot, not generalization.

## Pilot evidence

All 44,981 filtered-matrix barcodes were retained. They are CD138-selected specimen-level RNA barcodes; the public records provide no per-barcode malignant/clone labels. The data contain nonnegative int32 counts for 33,538 Gene Expression features, with one exact feature for each locked target. We report transcript detection only. No barcode is joined between biopsies.

| Participant | Pre→post GSM | Cells pre→post | Median Gene Expression UMIs pre→post | Fraction post cells below D=1,318 |
|---|---|---:|---:|---:|
| PM01 | GSM5971019→GSM5971023 | 6,207→5,028 | 12,293→13,135.5 | 4.93% |
| PM07 | GSM5971021→GSM5971025 | 5,787→4,462 | 5,761→5,553 | 2.60% |
| PM10 | GSM5971020→GSM5971024 | 9,792→13,705 | 8,591.5→1,318 | 49.96% |

The common pseudo-depth is the floor of the minimum of the six specimen medians. The six sample records say 10x CG00052 v2 while the series description says v3; this discrepancy remains unresolved. All 100 replicates jointly sampled the three target counts and all other RNA counts without replacement. Replicate bands describe within-specimen thinning sensitivity, not participant uncertainty.

For each pair, none of the two directed conditional changes had the same direction across all three participants. For example, SLAMF7 detection among TNFRSF17-undetected barcodes changed by +0.056, +0.047, and −0.158 pre-to-post in PM01, PM07, and PM10. After joint thinning, the replicate medians were −0.003, +0.012, and −0.062; PM01's direction reversed and the PM10 magnitude shrank. For SLAMF7 detection among GPRC5D-undetected barcodes, the raw changes were −0.014, +0.038, and −0.403, with depth-sensitivity medians −0.006, +0.012, and −0.063. The fixed-panel distributions differ across specimens, but this pilot does not identify a consistent participant-level pattern beyond those observations.

The very low PM10 post depth makes its large raw shifts especially sensitive to sampling. Joint thinning narrows these changes but does not make the data comparable in tumor composition, sequencing chemistry, or cell state. It is a sensitivity check, not a correction for those confounders.

## Why full execution fails

- The independent unit count is three participants from one study. The 44,981 barcodes are subsamples, not independent patients; no independent validation cohort was identified.
- Participant-level directions vary, and one specimen has a marked depth imbalance. No reproducible pattern or precision/validation justification supports a broad biological conclusion.
- Per-barcode malignant/clone/state labels, matched surface protein, functional antigen measurements, lineage tracking, and exact chemistry confirmation are unavailable in the checked records; the sample-level source does identify a CD138-positive-selected fraction.
- The pre/post contrast is venetoclax-relative, not antigen-directed therapy. There is no untreated comparison and no basis for a causal treatment-effect claim.
- The near prior art already establishes broad same-cell target co-expression and longitudinal target expression. An exact three-pair table in another treatment context is too narrow to establish a substantial scientific contribution on its own.
- No selection method is tested; therefore a baseline or method-superiority claim is not appropriate. A prediction holdout would be impossible with three participants and would not validate this exposed descriptive contrast.

These are limits of this candidate's evidence, not proof that all serial myeloma RNA data are unavailable.

## Closest primary work and comparison

The 2023 study “Single-Cell Discovery and Multiomic Characterization of Therapeutic Targets in Multiple Myeloma” reports same-cell co-expression/mutual exclusion, flow checks, and longitudinal target expression across six people; it includes TNFRSF17, GPRC5D, and SLAMF7. The 2024 study “Integrated epigenetic and transcriptional single-cell analysis of t(11;14) multiple myeloma and its BCL2 dependency” profiles pre/post-venetoclax single-cell RNA/ATAC and analyzes cell-state changes in this disease context. The accessible main text does not state the exact joint-category endpoint used here, but its supplement was not fully inspected. These sources bound the narrow G1 rationale and prevent a generic novelty claim. [Cancer Research DOI 10.1158/0008-5472.CAN-22-1769](https://doi.org/10.1158/0008-5472.CAN-22-1769); [Blood DOI 10.1182/blood.2023020276](https://doi.org/10.1182/blood.2023020276). Computational target selection and antigen escape are also established, so neither is claimed as new. [“Single-cell mapping of combinatorial target antigens for CAR switches using logic gates,” DOI 10.1038/s41587-023-01686-y](https://doi.org/10.1038/s41587-023-01686-y); [“Mechanisms of antigen escape from BCMA- or GPRC5D-targeted immunotherapies in multiple myeloma,” DOI 10.1038/s41591-023-02491-5](https://doi.org/10.1038/s41591-023-02491-5); [“Multimodal antigenic escape to GPRC5D-targeted T cell engagers in multiple myeloma,” DOI 10.1038/s41591-025-04175-8](https://doi.org/10.1038/s41591-025-04175-8).

## Recommendation

Stop work on this cohort after this reproducible pilot and do not present it as a validated method, biological mechanism, or competitive-strength extension. If a new direction is considered, first complete only the metadata tests in `adjacent_gap_data_cards.md`; do not download more data or start a new analysis until one card passes its stated gate and the student/advisor chooses the scope.
