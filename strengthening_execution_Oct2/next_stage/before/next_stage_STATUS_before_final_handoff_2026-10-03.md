# Next-stage execution — October 2, 2026

Plan: `docs/superpowers/plans/2026-10-02-donor-reference-extension-gates.md`.

Status: orientation and preservation gate complete; new bounded retrieval and diagnostics underway.

Authorization: public primary-source research/public-data retrieval and reversible workspace analysis only. No outside messages, clinical label changes, publication, registration, accounts, lab access, or student submission prose. Deadlines are late January (school fair) and mid-March (state); weekly hours and an independent reviewer are not identified.

Starting preservation: every entry in the existing 63-file protected-input manifest was checked before new diagnostics. A new hash snapshot and copies of root/compact status pointers are under `before/`. No Git repository or local `AGENTS.md` was found. Existing files remain read-only; all new work is here.

Current verified metadata, prior to numerical diagnostics: local matched cell metadata has 29,387 cells and fields `index`, `sample_ID`, `disease_stage`, `n_genes`, `frac_mito`, `n_counts`, `cluster`, and `normal_or_neoplastic`. Group-level counts have not yet been reconciled to the loader's exact filter. Local metadata does not directly encode age, sex, preparation, or sequencing batch.

Next: confirm exact eligible cells/count units; freeze diagnostic protocol; retrieve sample-level source metadata and bounded independent data candidates; run locked tests/diagnostics; complete the specific BCMA and GPRC5D clinical extraction; document contribution gate and path decision.

## Resume checkpoint — October 3, 2026

- Reopened only the bounded donor-metadata, external-data, prior-art and selected clinical-source questions in the finalized execution plan. The existing within-study B analysis is retained as executed; its headline is a measurement limitation, not a validated clinical or novel biological result.
- The previously generated `analysis/results/external_annotation_by_biopsy.tsv` is not a valid biopsy summary: it groups on `Biopsy.Sequence` alone, collapsing the same sequence number across distinct `Hospital.Code` values. Preserve it as a superseded artifact; do not use its rows or any counts derived from them.
- The downloaded GEO count file `GSM9189168_AB14564.txt.gz` contains 384 cell columns. Exact join to the workbook yields 105 annotated cells, including only 2 labelled `Malignant`; this file verifies a working retrieval and partial join, but cannot support donor-policy evaluation.
- The workbook includes `Hospital.Code`, `Biopsy.Sequence`, `Amp.Batch.ID`, `Seq.Batch.ID`, `Disease`, `Tissue`, `Method`, and `Populations`. Primary documentation inspected so far does not define the two identity fields as a validated person key; treat donor grouping and overlap with GSE193531 as unknown pending provenance.
- Next gate: construct an auditable per-batch/per-Hospital.Code+Biopsy.Sequence inventory; use only GEO source metadata to choose a deterministic healthy bone-marrow plasma-cell and malignant bone-marrow subset; verify patient grouping/overlap and exact count-file joins before any marker-value analysis. The current sample pilot is not evidence that the broader dataset is unusable.
- Clinical repair sources were freshly checked for ide-cel (PMID 33626253) and talquetamab (PMID 36507686); combined grade-count strings and the agent-specific denominator unit need normalization into separate integer event/grade rows before they can be used. Keep percentage-only reports as percentages and do not infer exact counts.

