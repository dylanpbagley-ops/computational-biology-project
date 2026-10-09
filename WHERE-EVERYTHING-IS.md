# WHERE EVERYTHING IS — read this first

**Purpose.** A map of this project, written September 8, 2026 (updated September 9) so that a new
session never has to hunt. If you are starting fresh: read `EV-Capture-CURRENT.md`, then this
file, then section 29 of the results document. That is enough to work.

**One-line orientation.** This project proves there is a hard ceiling on what any
antigen-targeted cancer therapy can achieve, set by how the target is distributed across normal
cell types — and shows that ceiling predicts which real drugs poisoned patients.

---

## The four documents that matter

| file | what it is | read it when |
|---|---|---|
| **`EV-Capture-CURRENT.md`** | the living handoff. Small on purpose. | always, first |
| **`ev_capture_model/docs/Aim13_Results_Two_Comparator_Ceiling.md`** | **THE definitive results record, sections 1–33.** Every number, every method, every recorded error. ~117 KB. | whenever you need a number or a method |
| **`PROCESS.md`** | how the handoff system works and how to maintain it | before editing `EV-Capture-CURRENT.md` |
| **`WHERE-EVERYTHING-IS.md`** | this file | when you can't find something |

**Rule:** the results document is authoritative. If `CURRENT.md` and the results document
disagree, the results document wins, and `CURRENT.md` needs fixing.

---

## Results document — what is in which section

| § | topic |
|---|---|
| 1–8 | the single-compartment ROC ceiling: derivation, the n* identity, numerical verification |
| 9–12 | placing real agents on the ROC; figure corrections |
| **11a** | **a third real calibrated agent (Bluhm et al., BCMA), and the same calibration gap found again independently in a second disease/technology** |
| 13–15 | the third normal compartment; gating robustness; **myeloma replication** |
| 16–17 | healthy-donor CITE-seq reference; the **two-antigen (multi-target) theorem** |
| 18–22 | the outcome flip; circularity removed; **a wrong diagnosis of mine, recorded**; the outcome test at its ceiling |
| **23** | **the whole-body extension — 19 safe vs 10 toxic, 0.868, p = 0.00031; two-factor rule 0.905; solid tumours alone 0.792** |
| 24 | the pancreatic thread end to end, including an honest negative |
| 25 | shuffled-label null for the two-antigen work |
| 26 | the protein cross-check, **the wrong control I ran**, and the corrected resolution conclusion |
| 27–28 | independent solid-tumour replication (rapid autopsy) + the two-antigen work there |
| **29** | **the weighted multi-compartment ceiling theorem — the unification. Lead with this.** |
| 30 | genome-wide discovery screen + its 57-fold positive control |
| 31 | why Deep Visual Proteomics cannot answer the protein question |
| 32 | the protein question answered as far as any data allows — **a limitation, presented as one** |
| **33** | **an independent healthy bone-marrow reference for the myeloma side (Boiarsky et al.), closing the asymmetry with the CD19/PBMC healthy reference in §16** |
---

## Code — `ev_capture_model/`

```
ev_capture_model/
  src/                  the model itself (binding_model.py, population_model.py, aims 1-12)
  src/aim13_scoping/    ALL the ceiling work (sections 13-33). 56 scripts.
  src/figures_writeup/  every figure script + figure_style.py (shared visual identity)
  src/pancreatic_calibration/   CLDN18.2 and mesothelin dose-response fits
  src/cptac_msln_outcomes/      CPTAC PDAC survival data (n = 135)
  tests/                341 tests. Run: python -m pytest -q  (~2.5 min)
  docs/                 all write-ups and archives
```

### The scripts that matter most, in `src/aim13_scoping/`

| script | what it does |
|---|---|
| **`target_labels.py`** | **THE FROZEN LABEL FILE.** 40 targets, clinical safety history, one citation each. No analysis script may write to it. Contains the word-boundary guard that failed once. |
| **`verify_documented_numbers.py`** | **run this after ANY change.** Re-asserts all 22 headline numbers against the result JSONs and fails loudly on drift. (Sections 11a and 33 are standalone narrative findings, not wired into this checker — same precedent for both.) |
| `weighted_ceiling_theorem.py` | §29 — the theorem and its 16,400-check verification |
| `wholebody_v2.py` | §23 — the main whole-body analysis, all six tests |
| `pancreatic_thread.py` | §24 — calibrations, CPTAC survival, the 14-target ranking |
| `two_antigen_null_control.py` | §25 — shuffled-label null |
| `protein_cross_check.py`, `matched_resolution_check.py` | §26 — protein references; the second is the CORRECT control |
| `resolution_sensitivity.py` | §26.2 — **the WRONG control, kept deliberately** |
| `adc_autopsy_ceiling.py`, `adc_two_antigen.py`, `trial_registry_control.py` | §27–28 — the solid-tumour replication |
| `target_discovery.py` | §30 — the discovery screen and its positive control |
| `hpa_summary_reference.py` | §31 — the DVP detection-limit finding |
| `protein_celltype_citeseq.py`, `adt_gene_map.py` | §32 — protein vs RNA confinement at matched resolution |
| **`bcma_operating_point.py`** | **§11a — places Bluhm et al.'s BCMA n*=100 on the myeloma CITE-seq data; finds the same absolute-calibration gap independently** |
| **`extract_boiarsky_markers.py`** | **§33 — one-time extraction of 11 target-gene rows + cell metadata from the raw GSE193531 matrix into a compact JSON (`extracted_markers.json`); no script needs to re-touch the 52MB raw matrix again** |
| **`boiarsky_healthy_reference.py`** | **§33 — the healthy-NBM BCMA restriction check itself, run on the extracted JSON** |
| `prospective_predictions.py` | the **locked** forward predictions, dated 2026-09-08 |
| `build_confinement_table.py`, `targetscan.py` | the released tool: genome-wide table + CLI |

Every script writes a matching `.json` next to itself. Those JSONs are what the verifier reads.
---

## Data — `ev_capture_data/` (about 4.5 GB total)

| folder | size | what it is | used by |
|---|---|---|---|
| `raw_cytof/` | 1.7 GB | **41-patient DLBCL CyTOF**, 6,865,361 cells. The Aim 12 cohort. | Aims 12, 19–21 |
| `raw_cytof_rln/` | 2.0 GB | **36 reactive lymph node (RLN) files**, 8,143,673 cells — the normal-B comparator | §19–21 |
| `gated_cache/` | 595 MB | 77 gated `.npz` files (B / T / other compartments) built from the two above | §13–21 |
| `gated_cache_v1_2compartment/` | 407 MB | **superseded** two-compartment version. Kept for reproducibility; do not use. | — |
| `healthy_pbmc_citeseq/` | 146 MB | **8 healthy donors, 161,764 cells, 228 surface proteins by antibody**, 31 blood cell types. A genuine protein cell-type reference. | §16, §32 |
| **`healthy_bone_marrow_boiarsky/`** | **55 MB** | **Boiarsky et al. 2022 (GSE193531): 29,387 cells, 9 of them genuinely healthy (NBM) donors, CD138-sorted/plasma-cell-enriched. Raw matrix + metadata (both gzip), plus the extracted 11-gene JSON.** | **§33** |
| `normal_tissue_atlas/` | 44 MB | **the six HPA reference files** — see below | §23, §26, §30–32 |
| `adc_autopsy/` | 9.1 MB | **rapid-autopsy solid-tumour cohort**: 30 patients, 749 metastases, 96 normals, 72 ADC targets, plus a 623-trial registry | §27, §28, §30 |
| `myeloma_adt/` | 4.8 MB | 11 myeloma samples, 73,755 cells, 63 proteins | §15, §11a |
| `cd19_populations/`, `raw/`, `raw_svatlas/` | 71 MB | EV-era source data (Aims 1–9) | Aims 1–9 |
| `*_digitized.csv` (root) | small | figure-digitised CA19-9 / EGFR dose-response points | Aims 3–4 |

### The six HPA files in `normal_tissue_atlas/` — which is which

| file | what it measures | resolution | result |
|---|---|---|---|
| `rna_single_cell_type.tsv.zip` | RNA | **154 cell types** ← the primary reference | **0.868** |
| `normal_ihc_data.tsv.zip` | protein, antibody staining | ~96 tissue × cell type | **0.785** |
| `rna_tissue_consensus.tsv.zip` | RNA | 51 organs | 0.726 |
| `ms_tissue.tsv.zip` | protein, mass spec | 20 organs | 0.650 (fails) |
| `ms_tissue_sample_data.tsv.zip` | protein, mass spec | 60 samples / 20 organs | 0.650 (fails) |
| `proteinatlas.tsv.zip` | **per-gene annotation, 119 columns** — protein class, subcellular location, cancer + cell-line specificity, prognostics for 21 cancers | one row per gene | powers §30–31 |

**Do not download more HPA files.** `proteinatlas.tsv.zip` turned out to contain everything the
Deep Visual Proteomics and TCGA downloads were wanted for, and §31 explains why the DVP data
cannot answer its question regardless.

**Downloading in general:** the container and the device shell are both blocked from
proteinatlas.org, nature.com and NCBI/GEO (both `ftp.ncbi.nlm.nih.gov` and `www.ncbi.nlm.nih.gov`
return a hard proxy-allowlist block — confirmed directly with `curl -v`). Downloads go through
Claude in Chrome, or — when that also hits a bot-check (PMC/reCAPTCHA, Cloudflare) — the Claude
Browser pane, which has its own separate, unrestricted network path. Getting a *large* file (tens
of MB) from either browser context onto the actual Windows filesystem has repeatedly been the
hard part, not finding it: browser-side blob-download tricks and base64 round-trips through tool
calls have both failed or been impractical for files over a few MB (see the Boiarsky data,
below). **The one approach that has reliably worked for a large browser-fetched file: have Dylan
download it himself in the browser and drop it in the outputs folder** — this is what actually
unblocked GSE193531 (§33) after in-session approaches stalled. PowerShell `Move-Item` commands
and direct download URLs have both failed repeatedly in practice too — **give Dylan a click path**
(which page, which link text) instead, or ask him to download and place the file directly.
---

## Write-up materials

`Poster, Paper, and Write Up/Figures/` — every figure as **300 DPI PNG + editable-text SVG**:
- `fig0`–`fig8` the original set (graphical abstract, superselectivity heatmap, CA19-9 flip, stem-cell crossover, Aim 11 reordering, Aim 12 CyTOF example, distribution skew, conceptual schematic, evidence-strength comparison), plus poster variants of fig0/1/5/8.
- **`fig9_wholebody_restriction`** — the whole-body result and the solid-tumour-only check.
- **`fig10_pancreatic_thread`** — calibration, the negative survival result, the 14-target ranking.
- `figure_style.py` in `src/figures_writeup/` holds the shared visual identity. Use it for any new figure.

## Documents worth knowing about, in `ev_capture_model/docs/`

| file | what it is |
|---|---|
| `Aim13_Results_Two_Comparator_Ceiling.md` | **the results record, §1–33** |
| `Label_Review_Sheet.md` | 40 labels with citations and tick boxes, for an outside reviewer to sign |
| `OSF_Preregistration_DRAFT.md` | pre-registration incl. the H5/H6/H7 outcomes as they actually came out |
| `Project_Status_Weaknesses_and_Expansion.md` | first status assessment (Sept 7) |
| **`Full_Assessment_Impact_ISEF_Weaknesses_NextSteps_Sept8.md`** | **the Sept 8 full weakness/impact assessment (22-item weakness table, ISEF odds, phased next-steps plan). Moved here Sept 9 from a `Claude outputs/` inbox copy — this is now the canonical location, matching the Claude project.** |
| **`Work_Completed_Sept9.md`** | **the Sept 9 tracking doc — what's been done from the Phase 0-2 plan above, what's blocked, walkthroughs for submitting the locked predictions and the pre-registration. Read this before re-doing anything in the assessment doc above. Moved here Sept 9 for the same reason.** |
| `archive_sessions_sept3-8_2026.md` | **all session narrative and every recorded error, through Sept 8** |
| `ev_capture_model_narrative_archive_*.md` | four older archives (pivot/Aim 12, pancreatic search, figures, novelty diagnosis) |
| `CAR_T_Pivot_Validation_and_Findings.md` | the precisely-scoped account of Claim 2 |
| `Prediction_Registry.md` | four falsifiable predictions stated before checking |
| `AI_Use_Disclosure_Log.md` | ISEF requires this. Keep it current — **still stops at Aug 26-27, does not yet cover the Sept 6-9 work. Open gap, needs Dylan's own review.** |
| `Limitations_and_Calibration_Gap.md` | the calibration gap, stated honestly |

Also in the **Claude project** (visible across Claude surfaces, not on disk):
`Status_ISEF_Ratings_and_Expansion_v2.md` — the ISEF odds, the weakness ranking, and the
expansion angles with their evidence. (The assessment and tracking docs above are now saved in
**both** places — the Claude project, and `ev_capture_model/docs/` on this computer — kept in
sync going forward.)

---

## How to check nothing is broken

```bash
cd ev_capture_model && python -m pytest -q                       # 341 tests, ~2.5 min
cd src/aim13_scoping && python verify_documented_numbers.py      # 22 headline numbers vs the JSONs
```

Both should pass before and after any change. The verifier exists because this project has
already been bitten once by prose drifting away from code (§20).

## Loose items in the root folder, so nobody wonders

The root folder was cleaned up September 9, 2026 — see "What changed in the outputs folder,
September 9" below for the full account. As of that cleanup:

| item | what it is | keep? |
|---|---|---|
| `ChatGPT_Transfer_Package.zip` + folder (~137 MB combined) | a one-off export made to hand the project to another tool | **still flagged safe to delete** once Dylan is sure he doesn't want it — nothing new needed from it this session either |
| `GSE234261_BCMA_CAR-T_data/` | the BCMA CAR-T generalisation dataset (10 patients, directionally consistent, not significant at 5 vs 5) | keep — cited in the write-up as a supporting check |
| `archive/` | handoff versions v3–v15 | keep, never read by default |
| `data_sources/` | source PDFs and figure crops behind the digitised curves — **now also holds the Da Vià et al. 2025 PDF** (`blood_bld-2024-025643-main.pdf`, moved here from root Sept 9 — the article behind the item-9 search's EGA-blocked candidate, see §33) | keep |
| `ev_capture_data/tool_outputs/confinement_web.json` | the genome-wide confinement table the published tool embeds | keep |
| `ev_capture_data/gated_cache_v1_2compartment/` (407 MB) | superseded two-compartment gating | keep for reproducibility, do not use |
| `superseded_matched_resolution_test.py.bak` | renamed out of the way — pytest auto-collects `*_test.py` and this file's helper takes arguments, which caused a spurious collection error | keep as a `.bak`, do not restore the old name |
| **`_to_delete/`** | **a folder created Sept 9 holding items Claude could not delete itself (no delete permission on this device) — see the delete list below** | **Dylan: review and delete this whole folder** |

### What changed in the outputs folder, September 9

Root-level clutter is now cleaned up:
- `blood_bld-2024-025643-main.pdf` (Da Vià et al.) moved to `data_sources/pdfs_and_figures/`.
- `Full_Assessment_Impact_ISEF_Weaknesses_NextSteps_Sept8.md` and `Work_Completed_Sept9.md`,
  which only existed as a `Claude outputs/` chat-delivery inbox copy (stale — missing the Sept 9
  status-pointer line the project version has), are now saved in their proper home,
  `ev_capture_model/docs/`, alongside every other write-up document, matching the Claude project.
- The stale `Claude outputs/` copies, and a one-off helper script left over from editing the
  results document (`append_section33.py` — its job is done; it's not a reusable pipeline script
  the way `extract_boiarsky_markers.py`/`boiarsky_healthy_reference.py` are), were moved into
  `_to_delete/` since this device session cannot delete files directly.
- **Files Dylan should delete himself** (Claude cannot — no delete permission on this device this
  session): everything in `_to_delete/` (safe, all superseded/disposable), and the empty
  `Claude outputs/` folder left behind after its two files were moved out. The `ChatGPT_Transfer_
  Package.zip` + folder remains a separate, larger, optional deletion — see the table above.

---

## Two standing constraints

1. **Never use the built-in Cowork browser** — it closes the Claude desktop app. Use Claude in Chrome.
2. **ISEF bars generative AI from writing the submitted abstract, research plan or poster text.** Dylan writes the prose; Claude does structure, figures, analysis, feedback and drafts he rewrites in his own words.
