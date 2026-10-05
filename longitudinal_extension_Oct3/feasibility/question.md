# Narrow question and scope

## Candidate question

Among the three t(11;14) participants with documented matched pre- and post-venetoclax single-cell RNA specimens in GSE199359, how do joint transcript-detection states for the predeclared targets TNFRSF17, GPRC5D, and SLAMF7 differ between the CD138-selected specimens? This is a small, descriptive treatment-context pilot of sampled-cell populations. It does not track cells or lineages, establish surface-protein availability, establish functional antigen escape, or estimate a treatment effect.

## Increment over the closest work

The closest myeloma target-discovery study (Yao et al., Cancer Research 2023, DOI 10.1158/0008-5472.CAN-22-1769) already measured co-expression/mutual exclusion among candidate surface targets in single cells and flow cytometry, and showed longitudinal expression changes in six people across disease progression, including TNFRSF17 and SLAMF7. Therefore neither generic same-cell target co-expression nor longitudinal target expression is new here. A separate study (Leblay et al., Blood 2024, DOI 10.1182/blood.2023020276) already profiled pre/post-venetoclax single-cell RNA/ATAC and resistance-associated cell-state changes in t(11;14) disease, but its reported focus is B-cell/plasma-cell phenotype and venetoclax biology, not joint transcript-detection states for the three predeclared therapeutic targets.

The only potentially unanswered contrast retained for a minimal pilot is whether the joint detection distribution of this fixed three-target panel differs between documented pre/post-venetoclax CD138-selected specimens in that distinct t(11;14) context. It is not a claim of a new target-selection algorithm, prediction method, therapy escape, or clinical benefit. The novelty gate is a narrow conditional pass for descriptive feasibility only; the closest papers substantially bound interpretation. Any claim that the exact three-target joint distribution was not shown in supplements remains provisional because the relevant supplementary tables/figures were not opened in this pass.

## Population and assay scope

- Independent unit: participant (three paired participants); specimens are paired within participant; cells are subsamples.
- Assay: deposited single-cell RNA. H5 data were verified as nonnegative int32 Gene Expression UMI counts; detection uses count >0.
- Selection: source paper describes CD138-positive selected cells. This is not a per-cell malignant-state annotation. If no per-cell tumor/clone labels are present in the deposited records, any pilot must be explicitly limited to CD138-selected specimen-level RNA and cannot call every cell malignant.
- Time: “pre” and “post” are relative to venetoclax in the 2024 source table, not to CAR/TCE therapy. The specimens are separate biopsies; no lineage tracking.
- Endpoint if the pilot proceeds: within each specimen, four joint categories for every fixed target pair; marginals; OR coverage; B detection among A-undetected cells, with undefined denominators reported as unavailable. No cell-level inferential p-values.

## Decision

G1 passed narrowly for the treatment-context descriptive pilot only; the prediction/selection branch failed. G2 passed for the three verified pre/post RNA pairs and six processed files. The locked pilot completed on 44,981 CD138-selected barcodes with all fixed pair categories and joint depth sensitivity (D=1,318). Participant-level changes were heterogeneous and the PM10 post specimen was much shallower; therefore the full strong-contribution gate failed. Do not start full analysis or claim a substantial extension from these three participants.

