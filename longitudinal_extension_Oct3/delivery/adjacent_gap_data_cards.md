# Adjacent gap/data cards (metadata tests only)

## Card 1 — Talquetamab progression with sample-resolved single-cell RNA

**Question if feasible:** In verified same-participant pre/progression specimens around GPRC5D-directed talquetamab, do tumor-state-labeled single-cell RNA profiles show a repeatable change in the joint target-detection distribution? This would be an antigen-directed therapy context, distinct from the venetoclax pilot, but RNA remains a surrogate for surface protein and cannot establish functional escape.

**Inspected primary records:** the 2026 GPRC5D escape paper (DOI 10.1038/s41591-025-04175-8; [primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC13004696/)) cites RNA accessions GSE199359 and GSE301286, ATAC accessions GSE199268 and GSE301287, and controlled WGS access through EGA. The GSE301286 GEO sample titles include four P-coded records; title strings do not resolve patient identity, modality per specimen, order, or therapy exposure. The paper reports paired WGS in a subset; that does not establish paired RNA or protein. A primary figure/supplement crosswalk was not verified because the direct supplement route was blocked in the bounded pass.

**Independent unit and feasibility gate:** participant, with separate biopsies nested within participant. First verify from a primary caption/table and GEO sample records at least three exact RNA specimen pairs, stable participant IDs, pre/progression order, intervening talquetamab exposure, actual RNA file links, and per-cell tumor/state labels. Keep WGS, ATAC, RNA, and protein as separate modalities unless a source explicitly links the same specimens. Stop if the only verified pairs are WGS or if identities/order remain unresolved.

**Novelty gate:** GPRC5D-targeted antigen escape is already established in this paper and elsewhere. The only candidate must be a specific result beyond republishing escape, such as a fixed-panel paired RNA joint-distribution change that is not already shown in the article/supplement and has a defensible assay comparison. Require the paper-to-GSM crosswalk and one exact supplement/figure overlap check before any analysis; if either remains inaccessible, stop rather than repeat blocked-route attempts.

**Next test:** one targeted metadata-only reconciliation of all GSE301286 GSMs to the paper's MM-03/MM-10/MM-20 modality and timepoint labels, then a single novelty check against accessible primary figures/supplement. No matrix download before both gates pass. If passed, return for a new scoped decision; do not assume this is a solved or novel extension.

## Card 2 — Cross-check single-cell RNA targets against surface protein

**Question if feasible:** In myeloma specimens with same-sample single-cell RNA and flow/surface-protein measurements, does the fixed-panel RNA joint-detection distribution correspond to measured antigen co-expression, and does that relationship persist across verified longitudinal samples?

**Inspected primary lead:** the 2023 myeloma target-mapping paper (DOI 10.1158/0008-5472.CAN-22-1769; [primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC10102848/)) reports flow-cytometry validation and longitudinal specimens among six participants; GSE213691 is cited as a data lead in the source inventory. The exact patient/timepoint-to-flow crosswalk and whether the flow measurements are matched to the same RNA specimens were not established in this pass. This is a metadata question, not evidence that a matched cohort exists.

**Independent unit and feasibility gate:** participant for serial comparison; specimen pair for same-sample cross-modality matching. Verify the GEO sample records and primary supplement for stable patient IDs, dates/order, sample-level RNA and flow/protein pairing, target-panel coverage, and eligible pair count. Require at least three ordered participant pairs with comparable assays for a descriptive feasibility pilot; do not interpret flow from a different specimen as matched protein.

**Novelty gate:** this source already reports single-cell target co-expression, longitudinal expression, and flow validation. Matched RNA/flow availability alone is not new. Proceed only if one precise additional contrast survives a single review of the exact primary supplement/figure details; otherwise treat this as prior art and stop.

**Next test:** one metadata and supplement check for GSE213691, if the primary table/supplement is accessible, then report a matched/unmatched table and the exact novelty overlap. If fewer than three verified paired participants, no same-specimen protein, or no remaining incremental contrast exists, stop. Do not retry blocked routes repeatedly. Any extension would need a new protocol and advisor discussion before analysis.
