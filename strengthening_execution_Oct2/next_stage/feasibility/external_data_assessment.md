# External data and independence assessment

The usable external lead was the MARS single-cell atlas in GSE305971, part of the GSE305967 superseries. The GEO annotation workbook contained 416,184 per-cell records; its selected locked frame comprised first-biopsy (`Biopsy.Sequence=1`), bone-marrow (`Tissue=BM`), MARS (`Method=MARS`) source labels: Healthy/Normal_PC reference cells, NDMM/RRMM/Malignant cells, and a separate Healthy/Malignant ambiguous-label stratum. Public count supplements were downloaded by exact `(GSM, Amp.Batch.ID, filename)` identity; 456 files totaling 306,858,728 compressed bytes were verified. The entire data directory remained within the authorized 500 MiB total cap. This was processed integer-molecule data, not a raw sequencing archive.

The published study is a distinct source and assay from GSE193531, but the public records available here do not provide a person-level crosswalk sufficient to verify patient non-overlap. `Hospital.Code` was used as a conservative putative-person key and kept all same-code batches together, not treated as globally verified identity. Thus the completed evaluation is an internally held-out, same-study code comparison. It can assess how a policy behaves among codes not exposed to the pilot; it cannot be called independent cohort replication or a direct transfer of numeric RNA thresholds between assays.

The GEO record and publication differ slightly in their reported study totals (GEO page: 322 MM plus 20 healthy individuals; article abstract: 341 total). This discrepancy is preserved rather than reconciled by inference. No patient-level deduplication against GSE193531 was attempted. The published study already conducts target prioritization; the exact per-patient z-score implementation could not be reconstructed from available Methods/code access, so the selected comparison is limited to policy outputs on the new assay.

Primary files and reproducibility records:

- Annotation workbook: `data/GSE305971_cells_annotation_for_geo_ver_20250306.xlsx`, SHA-256 `2fc3f9e19ab82097f6a50689de7bf4083aeb81aba442f263ea81afe8b150b94d`.
- Exact public count-file list and hashes: `analysis/results/external_full_actual_download_manifest.tsv`, SHA-256 `9275d41c22cc9359aa2c21b1b1b72e3fb2bad080534f229f4c5b146ae214f15f`.
- Metadata-only selection: `analysis/results/external_full_selection_lock.json`, SHA-256 `c67c4e8e2fccf1ca0b3a36e61a5064b49cc446b6e0a411b191fe3f845ddb07c0`.
- Exact join/missingness: `analysis/results/external_full_file_join_audit.tsv`, `external_full_unmatched_selected_annotations.tsv`, and `external_full_newcode_missingness_bounds.tsv`.

The first count sample GSM9189168 verified access and cell-key joins but was not sufficient by itself for evaluation; it contained 57,874 gene rows, 384 count columns, 105 exact metadata matches, and only two Malignant cells. The bounded full selected run superseded that feasibility-only sample. Two count-table partitions used a `W...` cell-ID namespace while workbook IDs used `WMC...`; absent a verified crosswalk, these selected annotations remain unmatched rather than being guessed into the analysis.
