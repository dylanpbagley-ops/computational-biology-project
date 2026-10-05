# Cycle 3 status — parent-review correction complete

## Outcome

Verdict: stop_current_question for the selected candidate direction. P01’s neural peptide-row access remains limited. The broader correction search surfaced P05 as a separate public processed resource; it does not establish a distinct unreported contrast or support a pilot.

Q11: captured 2026-10-05T02:46:08Z; immunopeptidom* AND (interferon OR IFNG); 51 total, 50 returned at offset 0, next offset 50, more available. Screened only title ranks 1–35.

Q12: captured 2026-10-05T02:53:17Z; (peptidome OR ligandome OR "HLA ligands") AND interferon; 1,859 total, 50 returned at offset 0, next offset 50, more available. Screened only title ranks 1–5. This is the captured repeat of the same exact query response; no further query or pagination.

Together, title_screen.tsv has 40 unique title-only records. P04, P05, and P06 are the three metadata/abstract candidates. Only P05 received closer article and file inspection.

P05’s Supplementary Data 2 contains 784 selected unique peptide strings and per-condition detection fractions from four preparations of one ECN90 line; two primary islet preparations in the article are untreated-only. Supplementary Data 7 contains 292 unique strings in ten published clusters with binary condition flags and no replicate/sample IDs. A fixed exact-string join against 74 JCI Table 5 peptide sequences found 0 matches; Data 2 is not exhaustive, so this is not evidence of absence or safety. The broad IFN-α and inhibitor-cluster comparisons are already reported.

Two official workbooks were downloaded after HTTP 200 size checks: 85,288 and 37,674 bytes (122,962 total). Per-file and aggregate caps passed. See P05_data_inspection.md, D02, and data_manifest.json.

Correction clock: start 2026-10-05T02:44:25Z; final status timestamp 2026-10-05T03:19:24Z; elapsed 34m59s, including 14m59s beyond the initial 20-minute limit for the parent-directed public-file access and final record/integrity checks. The correction start timestamp is captured; parent dispatch time was not.

Final integrity checks passed after the evidence update: cycle3 6 sources, 12 search rows, 40 unique titles, 2 cards, 2 downloads/122,962 bytes, 6 protected inputs, 13 output hashes; cycle2 validator passed with 14 sources, 12 searches, 6 cards, 2 gates, 25,465,927 download bytes, 8 protected inputs, and 18 hashes.
