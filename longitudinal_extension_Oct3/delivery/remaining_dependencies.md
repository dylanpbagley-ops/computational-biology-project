# Remaining dependencies and limits

- **Clinical review:** The existing project still has no independent clinical reviewer or adjudication; this stage does not revisit or change those prior clinical records.
- **Cell identity:** GEO sample records identify primary CD138-positive MM cell preparations, but no per-barcode malignant-state or clone mask was verified. Current tables are restricted to CD138-selected matrix barcodes.
- **Protein/function:** GSE199359 matrices provide RNA counts only. There is no matched surface-protein assay in these six files. RNA detection does not establish protein availability or functional escape.
- **Chemistry:** individual GSE records say 10x CG00052 v2 while the series design says v3. This remains unresolved and limits cross-specimen comparability.
- **Treatment timeline:** Table 1 supports pre/post-venetoclax order; precise interval and other intervening treatments are not recorded in the accessible pair metadata.
- **Prior-art boundary:** the accessible primary main text for the closest studies was reviewed, but some supplementary material was blocked or not verified. The exact fixed-panel endpoint's absence from every supplement is not established.
- **Independent units/validation:** only three participants from one study are included. No independent holdout, prospective cohort, or protein-validated extension is present.
- **Future GSE301286 lead:** participant/timepoint/modality matching and exact talquetamab order are unresolved. Its paired WGS subset must not be treated as paired RNA/protein.
- **Future GSE213691 lead:** whether the flow data are matched at same-specimen level to longitudinal single-cell RNA is unresolved.
- **School process:** Dylan/advisor should decide whether either metadata-only card merits a new scoped protocol. No new cohort download, external communication, publication, registration, paid access, or laboratory activity is part of this stage.
