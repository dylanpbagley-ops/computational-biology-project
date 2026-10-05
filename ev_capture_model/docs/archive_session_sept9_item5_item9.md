# Archive: item 5 and item 9 closed out — September 9, 2026

Full detail behind two lines in `EV-Capture-CURRENT.md`'s "What changed" section. Kept here per
`PROCESS.md` so `CURRENT.md` doesn't creep back up in size.

## Item 5 — a third real calibrated agent on the ROC (§11a)

Dylan found and confirmed the correct free full-text link for Bluhm et al. 2018 (BCMA CAR-T,
Molecular Therapy, PMID 30078440) via the browser pane, which got past the reCAPTCHA gate that
had blocked every earlier attempt (PubMed MCP full-text, WebFetch, ScienceDirect, Molecular
Therapy's own site all 403'd or returned empty full text — Elsevier-licensed, not in the open PMC
corpus). Read in full: BCMA density quantified with **QuantiBRITE PE calibration beads** — the
same absolute-molecule standard as CLN-978/blinatumomab — reporting an activation threshold "in
the range of 100 BCMA molecules" (JeKo-1, the lowest-density responding line). This is a genuine
methodological match to the project's existing rigor bar, unlike the B7H3 paper Dylan had
correctly flagged as "too approximate" in an earlier round.

Dylan's instruction: "if you think that data is accurate enough to be beneficial then let's add
to ROC."

**What actually got built.** BCMA isn't in the CD19/DLBCL CyTOF panel Section 11's table uses, so
this isn't a drop-in row. Wrote `bcma_operating_point.py`, applying the identical test Section 11
already ran for CD19 to BCMA in the myeloma CITE-seq data instead (CD269/BCMA raw ADT counts,
11 samples, same population gating as `myeloma_ceiling.py`). Honest result: **the same absolute-
calibration gap Section 11 found for CD19, replicated independently in a second disease and a
second technology (CITE-seq vs. mass cytometry).** Raw CD269 counts across the cohort top out at
357 (one sample, one cell); per-sample tumor-compartment medians mostly single digits — nowhere
near the 100-molecule threshold, and the raw counts don't clearly separate tumor from bystander T
cells in some samples, itself a sign they aren't tracking real molecule numbers.

Written up as new **Section 11a** in the results document, immediately after Section 11. Full
test suite re-run after the edit: 195/195 passing.

**Why this counts as a genuine positive result despite no operating point landing:** the same gap
appearing independently in two diseases, two technologies, two labs is stronger evidence the gap
is structural (how these single-cell assays are generated) rather than a fluke of one panel.

## Item 9 — an independent healthy bone-marrow reference for myeloma (§33)

Prior rounds (documented in `Work_Completed_Sept9.md`'s "Round 2" and "Round 3" sections) had
found four real candidates and narrowed to two live ones: Da Vià et al. 2025 (*Blood*, 21 donors,
54,419 healthy plasma cells, but EGA controlled-access) and Boiarsky et al. 2022 (*Nature
Communications*, GEO GSE193531, 9 NBM donors, freely downloadable). Dylan found the Da Vià PDF
himself and asked which was the right call.

**Da Vià verification.** Read the PDF in full. Confirmed correct article (DOI match). Corrected
two things from the round-3 summary: the real healthy-donor plasma-cell count is 54,419, not the
234,789 combined-across-all-groups figure previously cited; and the paper's own data-availability
section explicitly confirms EGA accession EGAC50000000232, access "directed to and fulfilled by"
the corresponding author on request. Closed as not implementable in this session's timeframe —
filed, not pursued further.

**Boiarsky judged worth adding, then implemented.** Dylan's instruction: "Is Boiarsky a good
addition to the project or are we better off without it — if we're better off adding it then go
ahead." Reasoning before touching code: every existing myeloma comparison (§15, §19) is within-
patient (tumor vs. normal from the same biopsy) — no independent, externally-collected healthy
check exists on that side, unlike CD19/DLBCL which already has one (§16, PBMC). That's a real,
if modest, gap worth closing even though Boiarsky is RNA-only, not CITE-seq/protein.

**Downloading it was the hard part.** This session's device shell can't reach any NCBI/GEO host
(`curl -v` confirmed a hard proxy-allowlist block on both `ftp.ncbi.nlm.nih.gov` and
`www.ncbi.nlm.nih.gov`). Tried fetching via the Claude Browser pane (separate, unrestricted
network access) and streaming the ~52MB gzip matrix in-browser via `DecompressionStream`, then
either (a) triggering a blob-URL download — files never appeared in the real Windows Downloads
folder despite being granted access to check — or (b) base64-encoding the decompressed content
for transfer through tool-call text output — 52MB decompresses to tens of millions of characters,
impractical to round-trip. Also tried a targeted in-browser line-by-line scan for just the 11
needed gene rows, using `window`-object state to persist a `ReadableStreamDefaultReader` across
separate tool calls (local `const`/`let` variables do NOT persist between calls in this browser
tool — only `window.*` properties do) — this was still running, having scanned 1,084 of 22,274
gene rows after reading 128MB of decompressed stream, when **Dylan downloaded both GEO files
himself in the browser and placed them in the outputs folder**, which is what actually unblocked
this.

**Once the files were on the device, the rest took minutes, not hours.** Verified gzip integrity
(clean), inspected the real structure directly with `zcat`/`awk`/`gzip -t` instead of guessing:
22,274 gene rows (a filtered/HVG set, not the full transcriptome), gene-symbol row names, 29,387
total cells matching the metadata file's row count exactly. All 11 target genes (TNFRSF17/BCMA,
CD38, SDC1/CD138, SLAMF7, GPRC5D, CD19, CD3D, CD3E, CD14, NCAM1, MS4A1) found on a direct `awk`
grep in under a minute — confirming the earlier in-browser scan hadn't hit a real problem, just
bad luck with row order (genes aren't alphabetically sorted; the found rows range from line 1546
to line 16962).

Wrote `extract_boiarsky_markers.py` (one-time extraction of the 11 gene rows + full cell metadata
into a 3.4MB JSON, so no script need re-touch the 52MB raw matrix again; verified matrix column
order exactly matches metadata row order via an explicit assert) and
`boiarsky_healthy_reference.py` (the actual analysis).

**Found an unplanned scope limit before analyzing anything:** canonical T/myeloid/B/NK markers
sit under 12% positivity across all 29,387 cells (most under 4%), while plasma-cell markers
CD38/SDC1 sit at 60-67%. **This file is already CD138-sorted/plasma-cell-enriched, not whole bone
marrow** — consistent with the study's design, but it rules out repeating Section 16's exact
"restriction against other lineages" AUC test (too few other-lineage cells left to test against).

**What the data could still test, and what it found.** Restricted to the 9,329 genuinely healthy
NBM cells: BCMA (TNFRSF17) detected in **82.87%** — more often than the plasma-cell identity
markers CD38 (73.78%) and CD138 (57.44%) themselves. Checked whether the 15% of cells negative
for both CD38 and CD138 were a real separate population: they're not (median sequencing depth
roughly half the rest — 5,627 vs. 10,288 median UMI counts — the signature of dropout, not a
different cell type; BCMA still detected in 72.6% of them). For contrast, GPRC5D (a known lower-
abundance myeloma target) comes back sharply lower (5.67%), reassurance the pipeline distinguishes
genuinely sparse from genuinely common expression within the same sorted population.

Written up as new **Section 33** ("Fifth addendum"), following the document's existing
conventions (limits stated before the result; explicit "closes the artifact question, not the
safety question" framing). Full test suite re-run: 195/195. `verify_documented_numbers.py`:
22/22 unchanged (§33's numbers aren't wired into that checker, same precedent as §11a — both are
narrative findings from standalone scripts, not numbers persisted to a checked JSON).

## Close-out actions taken the same session

- `Work_Completed_Sept9.md` and `Full_Assessment_Impact_ISEF_Weaknesses_NextSteps_Sept8.md` moved
  from a `Claude outputs/` chat-delivery inbox into their proper home, `ev_capture_model/docs/`,
  matching the Claude project's canonical copies (kept in sync going forward).
- The Da Vià PDF moved from the outputs-folder root into `data_sources/pdfs_and_figures/`.
- A one-off helper script (`append_section33.py`, its job done) and the two now-stale `Claude
  outputs/` duplicate files were moved into a new `_to_delete/` folder — this device session has
  no delete permission, so Dylan needs to remove that folder (and the resulting empty `Claude
  outputs/` folder) himself.
- `WHERE-EVERYTHING-IS.md` updated: new §33 row, new Boiarsky data/script rows, doc-location
  changes, and a note on what actually worked for large browser-fetched files (asking Dylan to
  download and place the file himself, after several in-session approaches failed or proved
  impractical for a ~50MB file).
