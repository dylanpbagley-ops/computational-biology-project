# Longitudinal target extension — final stage status

Last updated: 2026-10-03T23:13:55-07:00 (America/Phoenix)

## Outcome

The authorized feasibility pass and fixed three-participant descriptive pilot are complete. The sample identity/order and file-access gates passed narrowly for PM01, PM07, and PM10 pre/post-venetoclax RNA in GSE199359. The strong-contribution gate failed. This is not a substantial scientific extension, validated method, causal treatment result, or evidence of functional antigen escape. Stop analysis on this cohort.

The saved protocol SHA-256 is f32b335d2bc88a1571f4045452033d422049304b501738c032c966bb8c45a9ac, locked at 2026-10-03 22:18:08 America/Phoenix (2026-10-04 05:18:08 UTC) before any matrix download or expression parsing. Six public matrices total 100,468,485 bytes; file sizes and SHA-256 values are recorded in data/download_manifest.tsv. All six matrices passed structural checks. The run retained 44,981 filtered-matrix barcodes and produced 18 specimen-by-marker-pair summaries plus 1,800 joint-thinning sensitivity rows.

## Scientific limits

Primary Blood Table 1 links the RNA pairs PM01 P1393→P1555, PM07 P1389→P1479, and PM10 P1217→P1438. PM01's P1426 is a different DNA specimen and was not substituted. The analyzed sample material is CD138-positive selected, but no per-barcode malignant, clone, or state labels are available. RNA counts support transcript detection only; there is no matched surface protein, functional antigen assay, lineage tracking, or untreated comparator.

The paired directions vary across the three participants. PM10 post median Gene Expression depth is D=1,318 and 49.96% of its barcodes fall below D; the other post specimens have 2.60% and 4.93% below D. The reason for PM10's depth difference is unresolved. Joint hypergeometric thinning is a within-specimen sensitivity analysis, not an adjustment for tumor composition, chemistry, cell state, or a source of participant-level uncertainty. The GEO sample records say 10x CG00052 v2 while the series description says v3; the discrepancy remains unresolved.

The exact fixed-panel joint endpoint was not found in the accessible main-text review of the closest studies, but supplement overlap was not fully verified, so no novelty claim is made. Generic same-cell target co-expression, longitudinal expression, target selection, and antigen escape are established prior art. See feasibility/prior_art.tsv, feasibility/gap_gate.md, and delivery/path_decision.md.

## Preserved inputs and checks

preservation/verify_preservation.py passed with 301 file references checked, both pre-existing status suffix reconciliations, and zero failures. It confirms the 63 protected inputs and 155 prior outputs match their recorded hashes, with the two documented status suffixes reconciled; the prior output inventory SHA-256 is b63865ce9249efbaf19e3f275b5e5d773f2ffdb9c51730c8f376b51b7d8c4305. No previous project file was edited.

The parent independently reran the final pilot and independently recomputed the primary raw-H5 outputs and barcode table. It checked the joint thinning against exact combinatorial probabilities and verified paired sensitivity summaries from 9,000 replicate-delta records. The eight focused tests and figure generator passed. Details are in logs/assistance_log.tsv; runnable commands are in delivery/reproduction.md.

## Final files

- Evidence-by-claim table: delivery/claim_evidence_matrix.tsv
- Final decision and bounded pilot interpretation: delivery/path_decision.md
- Reproduction instructions: delivery/reproduction.md
- At most two conditional metadata-only leads: delivery/adjacent_gap_data_cards.md
- Stage inventory with SHA-256 and byte size for every other stage file: delivery/output_hashes.json (the inventory does not hash itself)

## Next actions

1. Review the pilot and stop decision with the student/advisor; keep any submitted prose student-authored.
2. If the team still wants to continue, choose only one conditional metadata/novelty card and run its bounded check if the primary source records are accessible. Do not start a full analysis unless both the data/identity and novelty gates pass and a new scoped decision is recorded.

The final decision is stored in delivery/path_decision.md for delivery grouping, rather than feasibility/path_decision.md as the plan's expected location; this location choice is recorded in the reproduction note and decision log. No full follow-up analysis or further download is recommended. Revisit an adjacent card only after its one bounded metadata reconciliation and exact prior-art overlap test pass.