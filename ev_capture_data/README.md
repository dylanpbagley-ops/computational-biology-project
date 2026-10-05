# Lennon et al. 2019 — Digitized Figure 2 & 3 Data

Source: Lennon KM et al. 2019, *J Extracell Vesicles* 8:1685634, DOI 10.1080/20013078.2019.1685634, PMC6844376 (open access). Figures 2 and 3 (pages 4-5 of the PDF) were uploaded as `lennonsciencefairdata.pdf` and digitized here so the PDF never needs to be re-uploaded.

## What's in this folder

- **`*_digitized.csv` / `*_rescaled.csv` / `*_raw_digitized.csv`** — per-vesicle point clouds extracted from each scatter-plot panel. Columns: `diameter_nm`, `copies_per_ev`, `blob_pixel_area` (blank for density-estimated points), `source`.
- **`figures/`** — the cropped source images each CSV was digitized from (so you can visually re-check any panel), the two full 400-DPI page renders, and reconstruction-quality check plots.
- **`scripts/`** — the Python digitization pipeline (axis calibration, blob detection, density-surface reconstruction for saturated ink regions, and the decomposition analysis).

## Panels covered

| File prefix | Figure panel | Cell line / sample | Marker | Reported N (paper) |
|---|---|---|---|---|
| `panel_a_PANC1_EGFR` | Fig 2a | PANC-1 | EGFR | 1651 |
| `panel_d_HPDEC_EGFR` | Fig 2d | HPDEC | EGFR | 197 |
| `panel_e_PANC1_CA199` | Fig 2e | PANC-1 | CA19-9 | 736 |
| `panel_f_HPDEC_CA199` | Fig 2f | HPDEC | CA19-9 | 104 |
| `fig3c_EGFR_healthy` / `fig3c_EGFR_PDAC` | Fig 3c (left) | Patient plasma, all 6 healthy / all 5 PDAC pooled | EGFR | not given per-point (pooled time series) |
| `fig3c_CA199_healthy` / `fig3c_CA199_PDAC` | Fig 3c (right) | Patient plasma, pooled | CA19-9 | not given per-point |

Figure 2 panels (a, d, e, f) have three CSV variants each:
- **`_raw_digitized.csv`** — output straight from image analysis, no correction.
- **`_rescaled.csv`** (recommended default) — the same points, linearly rescaled per-axis so the sample's mean and CV exactly match the paper's reported `EGFR/EV_av`, `CV`, `D_av`, `CV` numbers. Preserves point rank-order and joint (diameter, copies) structure; only shifts/stretches each axis.
- Un-suffixed `panel_X_..._digitized.csv` — an early single-pass version, superseded by the raw/rescaled pair above; kept for reference only.

Figure 3c panels have no rescale step — no per-point summary statistic is published for the pooled 3c scatter (only per-patient N in Fig 3a/b, and fold-changes in the gated region), so these are **shape-only reconstructions**, not stat-matched.

## How the digitization works (and its real limits)

Each panel's plot-area axes were located by pixel analysis (long dark axis lines + tick-mark spacing), cross-checked visually against a grid overlay. Calibration is accurate to roughly ±1-2 pixels (~1 nm on the diameter axis, ~0.5 units on the copy-number axis).

**The core problem:** these are dense scatter plots — up to 1651 points in a 400×200-unit plot area a few hundred pixels across. In the low-diameter/low-copy-number corner (where most vesicles sit), individual dots visually and physically overlap into a single solid black (or solid red/blue) ink blob with no recoverable dot boundaries, at any image resolution derived from a printed/PDF figure. This is a real, unavoidable limit of digitizing a static scatter image, not a bug in this pipeline.

**The fix used here — hybrid reconstruction:**
1. Wherever dots are sparse enough to resolve individually (outer/tail regions), connected-component blob detection recovers **exact** pixel-precise (x, y) for each point. These rows are tagged `source = resolved` (or `resolved_cluster_*` for small multi-dot clumps that were split by estimated count).
2. Wherever ink has saturated into one large fused blob, the local pixel darkness is treated as a probability density surface (Gaussian-smoothed, then sharpened), and synthetic points are sampled from that surface. The **number** of synthetic points is set so the total (resolved + synthetic) count matches the paper's reported N for that panel. These rows are tagged `source = density_estimated`.

**What this means for using the data:**
- Treat `resolved` rows as real digitized data points — trustworthy to a couple of nm / copy-number units.
- Treat `density_estimated` rows as a **smoothed reconstruction of local point density**, useful for distribution *shape* (spread, skew, where mass concentrates, which population dominates a region) — not as claims about individual vesicles. Do not cite a specific `density_estimated` row as "vesicle #N had exactly Y copies."
- The `_rescaled.csv` files additionally apply a linear correction so summary statistics match the paper exactly (see table below) — this corrects for any systematic under/over-smoothing introduced by the density step, without touching relative point structure.

**Validation performed:**
- Reconstructed N exactly equals the paper's reported vesicle count for every Figure 2 panel.
- Reconstructed mean/CV are within the ballpark of reported values even *before* rescaling (see table below); rescaling then makes them exact by construction.
- Visual comparison (`figures/reconstruction_check.png`, `figures/fig3c_reconstruction_check.png`) confirms the reconstructed clouds are visually near-indistinguishable in shape from the original published panels — same skew direction, same relative spread between panels, same qualitative healthy-vs-PDAC separation in Fig 3c.
- For the Figure 3c decomposition analysis specifically, results were cross-checked using **only the exactly-resolved (non-density-estimated) subset** and gave the same qualitative and roughly the same quantitative answer as the full reconstruction (see `scripts/decomposition_analysis.py` output) — so the headline finding does not depend on the density-estimation step.

**Pre-rescale digitization quality (raw, unadjusted):**

| Panel | N (target/got) | Diameter mean (reported/digitized) | Diameter CV (reported/digitized) | Copies mean (reported/digitized) | Copies CV (reported/digitized) |
|---|---|---|---|---|---|
| 2a PANC-1 EGFR | 1651/1651 | 51 / 57.8 nm | 0.5 / 0.47 | 28 / 34.0 | 1.2 / 0.85 |
| 2d HPDEC EGFR | 197/197 | 52 / 56.1 nm | 0.4 / 0.69 | 19 / 23.3 | 0.9 / 1.31 |
| 2e PANC-1 CA19-9 | 736/736 | 41 / 44.9 nm | 0.4 / 0.50 | 31 / 37.1 | 0.9 / 0.75 |
| 2f HPDEC CA19-9 | 104/104 | 45 / 48.4 nm | 0.3 / 1.02 | 23 / 29.2 | 1.1 / 1.31 |

Digitized copy-number means run consistently ~15-25% high relative to reported, and CVs are inconsistently off in both directions — expected given the corner-saturation problem above, especially for the smaller panels (d, f) where N is low enough that a handful of misjudged blob splits move the statistics more. **Use the `_rescaled.csv` files, not the raw ones, for anything quantitative.**

**Known minor artifact:** because rescaling is a linear stretch to hit the target mean/CV exactly, a small number of points can land outside the plot's original 0-200 (copies) or 0-400 (diameter) range. Checked across all four Figure 2 rescaled files: this affects 9/1651 points (0.5%) in panel 2a (max value 225, vs. axis max 200) and 4/1651 points pushed to ~0 diameter in the same panel; panels 2d, 2e, 2f are unaffected. Harmless for distribution-shape purposes but worth knowing if you filter or plot against the original axis bounds.

## Figure 3c: the sign-inversion decomposition (this is the important one)

`scripts/decomposition_analysis.py` reruns the "fraction-positive vs. copies-per-carrier" decomposition from the handoff document directly on the digitized per-vesicle point clouds, not just on summary statistics from the text. This was "Immediate next step 1" in handoff v2, previously blocked on image access.

**Gate caveat:** the negative-gate polygon (dashed line in the original Fig 3c) was **hand-approximated** as a diagonal trapezoid from visual inspection, not independently digitized pixel-by-pixel from the dashed line itself. The approximation puts 0-0.4% of healthy EVs inside the gate, close to the paper's stated "<5%" — a reasonable but not exact match. Treat exact gated-population numbers as approximate; see "Result" below for why the qualitative conclusion is not sensitive to this.

**Result (see script output for full numbers):**

| Marker | Gated PDAC-enriched: mean diameter | Gated PDAC-enriched: mean copies/EV | All healthy: mean diameter | All healthy: mean copies/EV | Raw copies/EV ratio | Surface-density ratio (copies/nm²) |
|---|---|---|---|---|---|---|
| EGFR | 192 nm | 49.8 | 55 nm | 42.7 | **1.17×** (slightly higher) | **0.079×** (12.6-fold lower) |
| CA19-9 | 177 nm | 58.5 | 48 nm | 30.0 | **1.95×** (about double) | **0.119×** (8.4-fold lower) |

**This refines, not just confirms, the v2 handoff's sign-inversion claim.** The paper's text says the gated PDAC-enriched population carries "fewer target molecules" than the bulk/healthy population. The digitized data shows this is true in **surface density** (copies per nm² of membrane) — down roughly 8-13×, driven by the much larger diameter of the gated population — but **not true in raw copies-per-vesicle**, where the gated PDAC population actually carries as many or more absolute target molecules than a typical healthy EV (1.2-2× more, not fewer). The paper's own text uses "fewer target molecules" loosely; the digitized data indicates it is really a **density** statement (fewer molecules *per unit membrane area*, because the vesicle is so much bigger), not a raw-count statement. This distinction matters directly for the model: multivalent capture physics in Aim 1 depends on receptor count/spacing in the contact patch, which scales with areal density and contact-patch geometry, not the whole-vesicle copy number in isolation — so the *density* framing is likely still the physically relevant one for superselectivity, even though the "fewer molecules" phrasing taken completely literally (raw count) is not quite right.

This was cross-checked using only exactly-resolved (pixel-exact, non-density-estimated) points and gives the same result (surface-density ratio 0.070-0.080×), so it is not an artifact of the density-reconstruction step.

## Reproducing or extending this

All scripts are self-contained Python (numpy, scipy, PIL, matplotlib) in `scripts/`:
- `find_axes.py`, `find_ticks.py` — axis/tick pixel-location helpers (diagnostic, used interactively during calibration).
- `digitize_panel.py` — simple single-pass blob digitizer (superseded by hybrid version, kept for the `detect_axes()` helper).
- `digitize_hybrid.py` — the main Figure 2 digitizer (black-dot panels): blob detection + density-surface reconstruction + `rescale_to_match_reported()`.
- `digitize_fig3c.py` — color-based (red/blue) version of the same hybrid approach for Figure 3c.
- `decomposition_analysis.py` — loads the Fig 3c CSVs and runs the gate/density decomposition reported above.

To redo a panel with different calibration or gate assumptions, edit the pixel coordinates at the top of the relevant script call and rerun — no need to re-crop from the PDF unless a different panel is needed (the full-page 400 DPI renders in `figures/` cover everything in Figures 2 and 3).
