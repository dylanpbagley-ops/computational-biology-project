> **October 2 correction notice:** This is a historical analysis record. For current scientific scope and unresolved issues use `strengthening_Oct2/Strengthening_Report.md` and `strengthening_Oct2/mathematical_scope.md` in the outputs root. Original scores and labels are preserved, not newly adjudicated. The autopsy comparison reuses the project outcome labels; it is not independent clinical replication. Nominal discovery-panel p-values do not account for label, selection, or full research-history uncertainty. Pre-correction copies are retained in `strengthening_Oct2/before/`.

# Aim 13-15 results: the therapeutic-index ceiling, with a real normal-B-cell comparator

**Date: September 6, 2026.** This is the results document for the proposal in
`Aim13_Therapeutic_Index_Ceiling_PROPOSAL.md`. Everything below was computed this session
from real data. Scripts and result JSON: `ev_capture_model/src/aim13_scoping/`.

**What changed since the proposal:** the proposal's single largest open gap was that the
comparator was bystander T cells rather than normal B cells. That gap is now closed. All 36
B-panel reactive lymph node (RLN) control files were downloaded from the Mendeley release
(2.0 GB, 8,143,673 real cells), verified, and analysed. The result is more interesting than
the proposal anticipated, and it changes the headline claim.

---

## 1. Data and verification

- **41 DLBCL tumor samples** (existing, 6,865,361 cells) + **36 B-panel RLN control files**
  (new: 31 individual donors + 5 pooled batch-control aliquots, 8,143,673 cells).
- All 36 new files parse with the project's existing FCS3.0 reader. **0 parse failures.**
- **Panel identity confirmed channel by channel.** 30 of 36 RLN files carry the identical
  65-stain panel as the tumor files. The 6 that differ (all batch CytB4) are missing only
  two *derived QC columns* (`bc_separation_dist`, `mahalanobis_dist`), not antibodies — and
  the CytB4 *tumor* files are missing exactly the same two, so tumor and control are
  consistent within batch. Every antigen channel used here sits at an identical parameter
  index in all batches and both sample types (verified explicitly for 19 channels).
- The dataset's own description states each batch run included an aliquot from a pool of 13
  RLN samples "to serve as an internal staining/acquisition batch control." Those 5 pooled
  files were downloaded and used as controls rather than discarded.

Gating is percentile-based throughout (intact by DNA, live by 195Pt, bead-excluded, CD45+;
then B = CD19+ or CD20+, CD3−; T = CD3+, CD19−CD20−). Percentile gates are invariant to the
~10x cross-batch signal-scale difference documented in the proposal. Cached gated
populations for all 77 files live in `ev_capture_data/gated_cache/` (407 MB).

---

## 2. Negative controls, run first

These decide how to read every number that follows.

**Normal-vs-normal, same batch, independent healthy donors (95 donor pairs, 37 antigens):**
median |AUC − 0.5| = **0.089**, 90th percentile **0.220**, max **0.440**.

Healthy people's B cells differ from each other substantially. This is the noise floor: an
apparent tumor-vs-normal window has to clear it to mean anything. The antigens with the
largest healthy donor-to-donor spread are CD23 (0.166), CXCR5 (0.154), CD20 (0.144), BTLA
(0.138), CD19 (0.136), HLA-ABC (0.130).

**Cross-batch normal-vs-normal (40 pairs):** median |AUC − 0.5| = **0.109**, i.e. only
**+0.020** above the within-batch figure. Batch effects contribute far less than genuine
donor biology. This is the empirical confirmation of the proposal's claim that AUC is
effectively immune to the batch-scaling confound — asserted there, measured here.

---

## 3. The retrodiction: the four approved-drug targets have no window against normal B cells

Median AUC across the 41 patients, tumor B cells versus each comparator:

| antigen | vs bystander T cells | vs normal B cells | approved drug in B-cell lymphoma |
|---|---|---|---|
| CD19 | 0.976 | **0.576** | axi-cel, tisa-cel, liso-cel, blinatumomab |
| CD20 | 0.971 | **0.554** | rituximab, mosunetuzumab, glofitamab, epcoritamab |
| CD79B | 0.819 | **0.489** | polatuzumab vedotin |
| CD22 | 0.913 | **0.359** | inotuzumab ozogamicin, CD22 CAR-T |

Against everything that is not a B cell, these four are near-perfect. Against normal B
cells they collapse to chance or below — CD19's 0.576 sits *inside* the healthy
donor-to-donor variation measured above.

This is exactly right, and it was not fitted to anything. Every approved CD19/CD20/CD22/CD79b
therapy causes B-cell aplasia and hypogammaglobulinemia, because these targets genuinely
cannot distinguish a lymphoma B cell from a normal one. The model recovers a universally
known clinical fact from single-cell data alone, with no clinical input.

CD22 at 0.359 is a second, sharper retrodiction: tumor cells carry *less* CD22 than normal B
cells, which is the documented reason CD22-directed therapy is antigen-density-limited.

---

## 4. The joint window, and the headline claim

A target must spare **both** normal compartments, so the achievable index is bounded by the
worse of the two. Taking the minimum:

- **Best joint window of any antigen in the whole 37-marker panel: 0.616 (CD80).**
- **Best joint window among the four approved-drug targets: 0.576 (CD19).**
- Definitively excluded by data — below chance on at least one comparator: CD22 (0.359),
  CD72 (0.376), CD23 (0.259), IgD (0.218), CD21 (0.189), CD27 (0.172), CD43 (0.215),
  CD48 (0.312), TRAIL-R1 (0.245), FasR (0.357).

> **No antigen in a 40-marker B-cell panel offers a therapeutic window against both normal
> compartments that clearly exceeds healthy donor-to-donor variation. The antigens that
> became approved drugs are not the ones with the best window. They are the ones whose
> unavoidable off-tumor compartment — the normal B-cell lineage — happens to be
> expendable.**

That reframes target selection in B-cell malignancy: it was never solved by finding a
selective antigen. It was solved by choosing a dispensable lineage. And it predicts, from
the same principle, why the identical playbook fails in solid tumors — there is no
dispensable lineage there. This is a statement the field's own framing (find a
tumor-specific antigen; engineer a sharper threshold) does not make.

---

## 5. Per-patient windows: the window is patient-specific, not target-specific

- CD19's joint window ranges **0.082 to 0.977** across the 41 patients (median 0.576).
- **CD19 is the best available target in only 4 of 41 patients.**
- The best available antigen differs across patients — 17 different antigens win across 41
  patients — with a median best-available joint window of **0.860** versus CD19's 0.576.
- **33 of 41 patients** have at least one antigen whose joint window exceeds the healthy
  donor-variation 90th percentile; 40 of 41 exceed the median.

The actionable reading is that antigen selectivity in DLBCL is a per-patient property, and a
single fixed target is a compromise for most patients. That is a real, testable,
previously-unstated claim.

---

## 6. The honest scope: this is an exclusion tool, not a discovery tool

**Read this before citing anything above as a target recommendation.**

The comparator set here is lymphoid only — bystander T cells and normal lymph node B cells.
The true off-tumor compartment is every tissue in the body. That has one clean mathematical
consequence and one blunt practical one.

**The clean consequence:** adding normal compartments can only *lower* a ceiling, never
raise it. The joint window is a minimum over compartments, and a minimum over a larger set
is no larger. So every number above is a valid **upper bound**, and an antigen that already
fails here is **definitively excluded** — that conclusion is safe and cannot be rescued by
better data. CD22, CD72, CD21, CD23, CD27, CD43, CD48 are excluded on this basis.

**The blunt practical one:** the converse does not hold, and section 5's per-patient
"winners" must not be read as recommendations. HLA-ABC wins for 7 patients and is MHC class
I, present on every nucleated cell in the body; PD-L1, CD80, CD86, CD24 and CD40 are all
broadly expressed on cell types simply not present in a lymph node sample. They score well
here only because the compartments that would sink them were not measured. Their ceilings
will fall, possibly to nothing, when more normal tissue is added.

So the correct claim is: **a necessary condition, rigorously computed, that a large fraction
of candidate targets already fail.** Section 4's headline survives this caveat completely,
because it is a statement about targets failing, not passing.

### Other limitations, unchanged or newly identified

1. **Monotonicity.** The ceiling theorem assumes kill probability increases monotonically
   with antigen count. The project's own CLDN18.2 pancreatic-line finding is a real measured
   non-monotonicity — state it as the boundary, do not bury it.
2. **The likelihood-ratio corollary is real and has bite.** Where the empirical likelihood
   ratio is non-monotone, a likelihood-ratio rule beats *every* density threshold. Measured
   on real data: CD27 by +0.79 Youden's J, CD48 by +0.38, CD38 by +0.10, while CD19, CD20,
   CD22 and CD79B gain exactly 0.000 (thresholding is already optimal for those). This is a
   criterion for when a NOT-gate is required instead of a CAR, and it connects to Aim 5.
3. **Gating is percentile-based and mine, not an expert's.** Robust across a 40th-70th
   percentile sweep (separation held at 1.000 in the proposal's ground-truth test), but a
   clustering-based alternative should be compared.
4. **Normal B cells are a mixed population** (naive, memory, germinal-centre) and donor
   subset composition varies. Part of the 0.089 noise floor is composition, not marker
   biology. Subsetting normal B cells is a real refinement.
5. **Tumor purity.** The tumor "B" gate does not separate malignant from residual normal B
   cells inside the biopsy. Light-chain restriction (IgKappa/IgLambda are both in the panel)
   would fix this; sample O6357 shows a 416:1 kappa:lambda ratio, so the signal is clearly
   there, but surface-Ig-negative DLBCL will limit coverage.
6. **Novelty is bounded, and the boundary is documented.** See section 7.

---

## 7. Novelty, as checked

**Prior art that exists.** Single-cell-RNA target prioritization pipelines: SCAN-ACT
(*Genome Medicine* 2025, doi 10.1186/s13073-025-01514-9), Dannenfelser et al. (*Cell Systems*
11:215-228, 2020, doi 10.1016/j.cels.2020.08.002), and a 2025 bioRxiv logic-gate paper
(doi 10.1101/2025.03.19.644074). Measured antigen-density activation thresholds. Affinity-tuned
CAR panels. The Lim lab's ultrasensitive density-sensing circuit (*Science* 371:1166, 2021,
doi 10.1126/science.abc1855).

**What is not there.** The 2025 logic-gate paper binarizes expression at "greater than 0",
reports a single (efficacy, safety) operating point, and states outright that it "did not
distinguish between tumor-targeting and tumor-killing efficacy" — precisely the link the
n\* identity supplies. A 2025 *Frontiers in Immunology* review of CAR-T mathematical modeling
confirms directly that no reviewed model treats antigen density as a distribution yielding an
escape fraction, derives any selectivity bound, uses ROC or likelihood-ratio reasoning, or
offers an antigen-ranking framework. PubMed returns **zero** results for Neyman-Pearson or
likelihood-ratio combined with immunotherapy target discrimination, and only three papers with
"therapeutic window" and "antigen density" both in the abstract.

**Relationship to the Lim lab's circuit, stated plainly.** Their ultrasensitive synNotch→CAR
circuit is a sharper *monotone* response to antigen density. A sharper monotone rule moves the
operating point along the ROC toward its corner; it cannot move the curve. So the theorem
predicts and bounds their result rather than competing with it. Verified numerically: across
**642 monotone kill functions** on 8 real antigens — the project's own `kill_probability` at
45 parameter settings plus Hill functions with cooperativity up to n = 50, i.e. exactly that
class of ultrasensitive circuit — the largest excess above the ROC convex hull was
**+2.3e-15**, floating-point zero. (6 further evaluations overflowed at extreme k_intra and
were not evaluable; the analytic proof via the layer-cake representation covers them.)

---

## 8. What to do next

1. **Light-chain purification of the tumor gate** (item 5 above) — cheap, uses data in hand,
   removes the most obvious methodological objection.
2. **Re-run Aim 12 on gated tumor cells** with the batch analysis, and report honestly how
   much of the 59.3% reordering survives.
3. **Express the ceiling as the achievable operating point of real agents**, using Aim 10's
   real calibrated thresholds (81 / 843 / 8,690 molecules/cell) to place actual drugs on
   each antigen's ROC. This is what turns a bound into a statement about real therapies.
4. **Add non-lymphoid normal compartments** if any public single-cell protein atlas allows
   it, to tighten the bound. Every addition lowers ceilings and strengthens the exclusion
   claims.
5. **Pre-register the ground-truth label list** before any further ranking claim, ideally
   with the criteria written down and checked by someone else.

---

# Addendum, same session: light-chain validation, the Aim 12 re-run, and a documentation error found in the project's headline number

## 9. Light-chain purification and the same-tube internal control

The panel carries IgKappa and IgLambda, so clonality can be called by the standard clinical
criterion and the minority light chain used as a **same-patient, same-tube, same-run**
normal control — a better comparator than the RLN files.

**A real bug, caught and fixed, recorded so it is not repeated.** The first implementation
set the kappa and lambda positivity cutoffs at the same *percentile* of each channel's own
B-cell distribution. That forces |kappa+| == |lambda+| exactly, by construction: every one
of the 72 samples returned K:L = 1.00 and the method was vacuous. The cutoff must come from
a genuinely marker-negative population. T cells in the same file are surface-immunoglobulin
negative, so their own IgKap/IgLam distribution is the correct internal negative control.

**Calibration check, which the method has to pass:** with the corrected cutoff, the 31 RLN
donors read K:L median **1.50**, range **1.00-2.82** — the textbook physiological polyclonal
range. The method measures real light-chain biology.

**Tumor calls** (restriction threshold = the more conservative of the clinical 3:1 / 1:2
convention and the measured healthy-donor interval, plus a requirement of ≥5% surface-Ig+
B cells so sIg-negative tumors are not called on a handful of events): 15 kappa-restricted,
4 lambda-restricted, 4 surface-Ig-negative, 18 not restricted. **18 patients** have both a
clone (≥500 cells) and a residual internal normal population (≥150 cells).

**The validation.** Comparing AUC(clone vs same-tube internal normal) against
AUC(clone vs batch-matched RLN normal), across 35 antigens (excluding IgKap/IgLam, which
define the gate and are therefore circular — they return exactly 1.000 and 0.000, the
required sanity check):

- Pearson r = **0.625**, Spearman = **0.685**
- median |difference| = **0.097**
- mean signed difference (RLN minus internal) = **+0.046**

The two comparators agree in rank and magnitude. The cross-sample RLN comparator **slightly
inflates** windows relative to the perfectly-matched internal control, which is the
conservative direction for this project's exclusion claims: real windows are, if anything,
*smaller* than section 4 reports.

**And the central retrodiction gets stronger under the better control:**

| antigen | same-tube internal control | RLN comparator |
|---|---|---|
| CD19 | **0.473** | 0.576 |
| CD20 | 0.575 | 0.554 |
| CD22 | 0.563 | 0.359 |
| CD79B | 0.717 | 0.790 |

Against the patient's own residual normal B cells, CD19 sits *below chance* — tumor cells
are marginally harder to distinguish from the patient's own normal B cells than a coin flip.
The strongest available control gives the strongest version of the result.

Largest comparator disagreements (lymph-node normal B cells differ from intra-tumoral
residual normal B cells): CD27 (+0.299), CD23 (−0.275), FasR (+0.236), CD49D (+0.232),
IgD (−0.220). These are subset-composition effects and should be stated as such.

---

## 10. Aim 12 re-run on gated cells — and a documentation error in the project's headline number

### 10a. What is and is not reproducible

`aim12_single_cell_cytof_reordering.py` contains the FCS parser plus a set of **hardcoded
result constants transcribed from a browser session**. It contains no function that computes
the reordering result, and the test suite does not check those constants. **The project's
single headline number has therefore never been reproducible from this repository.** That is
not a claim it is wrong — it is a statement that nothing here could confirm it. This session
checked what is confirmable:

- **The worked example reproduces exactly.** Recomputed escape fractions at
  CLN-978_cytolysis / kill_threshold 0.5, from the local files: O6357 19.12% (stored 19.1),
  Q5971 33.67% (33.7), O8155 36.05% (36.0), Q7998 44.68% (44.7), P7676 51.15% (51.1). Every
  value within 0.05 percentage points.
- **The stored 41-patient naive-by-mean ranking is bit-identical** to the ranking recomputed
  from the local files.

So the parser, the local data, and the escape-fraction computation are all correct and
agree with the browser session. The discrepancy is isolated entirely to the counting rule.

### 10b. The counting rule, resolved

The documented text describes the criterion as:

> "Reordered with a MEANINGFUL spread (≥2 percentage points of real escape-fraction
> difference **between the reordered patients**): 16/27 (59.3%)"
> ... "Largest single real spread observed: 81.8 percentage points (one real patient's
> predicted escape fraction 81.8 points higher than another's)"

Implementing that criterion literally — a grid point counts when some **reordered pair**
differs by ≥2pp in escape fraction — gives **13/27 = 48.1%**, largest reordered-pair spread
**53.1 pp**.

Testing alternative conventions found the one actually used. A grid point counts when the
ranking flips at all **and the spread between the most- and least-escaped patients in the
whole cohort is ≥2pp** — a condition on the cohort's full range, not on the pair that
reordered. That convention reproduces the stored numbers **exactly: 16/27 = 59.3%, largest
spread 81.8 pp.**

**So the arithmetic is right and the prose is wrong.** 81.8 pp is the cohort's full
escape-fraction range at that grid point, not the gap between two patients who swapped
places. The ≥2pp filter, as actually applied, is nearly always satisfied whenever the model
is not fully saturated, so it removes far less than the description implies.

**What to do:** report **48.1% (13/27), largest reordered-pair spread 53.1 pp** as the
figure matching the stated criterion, or keep 59.3% and describe the criterion accurately.
Do not keep the current pairing of number and description. Figure 8's four-tier panel
(8.65 / 5.56 / 48.10 / 59.26) and any caption citing 81.8 pp need updating. This is very
much better found now than by a judge in March.

### 10c. The effect survives gating

Recomputed under one fixed convention (any reordered pair, ≥2pp escape spread), across three
populations:

| population | n | grid points reordered | largest reordered-pair spread |
|---|---|---|---|
| ungated (what Aim 12 used) | 41 | 13/27 (48.1%) | 53.1 pp |
| B-gated | 41 | 14/27 (51.9%) | 42.3 pp |
| light-chain clone | 19 | 11/27 (40.7%) | 15.0 pp |

The reordering effect is **not an artifact of ungated whole-biopsy data** — it survives, and
is marginally stronger, after gating. It weakens at the clone level (40.7%, and the spread
drops sharply to 15.0 pp), which is the honest reading: part of the apparent magnitude came
from mixing cell populations, but a real effect remains in pure tumor cells.

### 10d. The continuum replacement for the grid

Using escape = F(n\*), the arbitrary 27-point grid can be replaced by the exact criterion:
do two patients' CDFs **cross** anywhere inside the real n\* range (12.3 to 28,868
molecules/cell)? If they cross, mean-based ranking is provably wrong for some achievable CAR
sensitivity; if not, mean-based ranking is safe.

| population | patient pairs whose CDFs cross |
|---|---|
| ungated | 241/820 (29.4%) |
| B-gated | 196/820 (23.9%) |
| clone | 18/171 (10.5%) |

This has no tunable convention, no grid, and no "meaningful spread" threshold to argue
about. It should replace the grid metric in the write-up.

---

## 11. Placing real agents on the ROC — an honest negative result

Applying the three real measured activation thresholds (n\* = 81, 843, 8,690 molecules/cell)
to the tumor and normal-B distributions:

| agent | n\* | CD19 tumor killed | CD19 normal B killed | window |
|---|---|---|---|---|
| CLN-978 | 81 | 45.1% | 34.7% | +10.4 pp |
| blinatumomab (cytolysis) | 843 | 0.1% | 0.0% | +0.0 pp |
| blinatumomab (IFN-γ) | 8,690 | 0.0% | 0.0% | +0.0 pp |

Only the most sensitive agent lands inside the data's dynamic range at all; the other two
fall off the top of the CyTOF signal distribution entirely. **This is the absolute-calibration
gap, now quantified rather than described:** the CyTOF ion-count scale (medians ~5-400)
cannot be equated with molecules/cell, so operating points are not currently computable.

That is a useful negative result, because it sharpens which claims are available. The
**ceiling** (an AUC) is calibration-free and stands. The **operating point** of a specific
drug is not computable without absolute calibration and must not be claimed. CD22's window at
n\* = 81 is **−20.6 pp** (more normal B cells killed than tumor), directionally consistent
with its AUC of 0.359, so the qualitative ordering survives even where the absolute numbers
do not.

---

## 11a. A third real calibrated agent, and the same gap found again independently (September 9, 2026)

Section 11 used the two real, absolute-molecule-calibrated agents this project had on file at
the time: CLN-978 and blinatumomab, both CD19-directed, both fit to a continuous dose-response
titration. A third real agent was added this session: **Bluhm et al. 2018** (*Molecular
Therapy*, doi:[10.1016/j.ymthe.2018.06.012](https://doi.org/10.1016/j.ymthe.2018.06.012), PMID
30078440), a BCMA CAR-T study that quantified antigen density with **QuantiBRITE PE
calibration beads** — the same absolute-molecule standard as CLN-978/blinatumomab — and
reported an activation threshold "in the range of 100 BCMA molecules" (JeKo-1, the
lowest-density line that still triggered lysis, carried ≤~100 molecules/cell; NCI-H929
carried ~25,000). One honest methods difference, stated plainly: CLN-978/blinatumomab were fit
to a continuous titration curve, giving a precise k_intra; Bluhm's ~100-molecule number is read
off a handful of discretely bead-quantified cell lines correlated with lysis/IFN-γ output, not
a fitted curve. It is real, peer-reviewed, absolute-molecule data — just a coarser derivation
than the other two.

BCMA is not measured in the CD19/DLBCL CyTOF panel section 11 uses, so it cannot be added to
that table directly. It **is** measured in this project's own myeloma CITE-seq cohort (section
16 onward) as CD269, so the natural test is the same one section 11 already ran, replicated
independently in the second disease: does the myeloma ADT scale support placing n\*=100 as an
operating point?

Raw (non-CLR) CD269 antibody-derived-tag counts across all 11 myeloma samples, tumor
(malignant plasma cell) compartment vs. the three normal compartments:

| sample | n tumor cells | tumor median | tumor p90 | tumor max | T median | myeloid median | NK median |
|---|---|---|---|---|---|---|---|
| 01 | 420 | 10.0 | 19.0 | 32.0 | 14.0 | 6.0 | 6.0 |
| 03 | 692 | 3.0 | 8.0 | 18.0 | 3.0 | 3.0 | 3.0 |
| 05 | 762 | 2.0 | 5.0 | 9.0 | 2.0 | 1.0 | 1.0 |
| 07 | 567 | 3.0 | 6.0 | 11.0 | 3.0 | 2.0 | 2.0 |
| 12 | 576 | 11.0 | 16.0 | 27.0 | 7.0 | 7.0 | 6.0 |
| 13 | 869 | 6.0 | 9.0 | 41.0 | 4.0 | 3.0 | 2.0 |
| 15 | 271 | 5.0 | 10.0 | 54.0 | 4.0 | 5.0 | 3.0 |
| 17 | 737 | 3.0 | 7.0 | 357.0 | 2.0 | 3.0 | 2.0 |
| 19 | 136 | 2.0 | 3.0 | 10.0 | 0.0 | 1.0 | 0.0 |
| 22 | 212 | 1.0 | 4.0 | 7.0 | 2.0 | 1.0 | 1.0 |
| 24 | 759 | 3.0 | 8.0 | 20.0 | 2.0 | 4.0 | 2.0 |

The highest single-cell raw CD269 count observed anywhere in the 11-sample cohort is 357; the
median of per-sample tumor-compartment 90th-percentile counts is 8.0. Bluhm's n\*=100 sits
above nearly every value in this table, tumor or normal. **This is the same absolute-calibration
gap section 11 reported for CD19, found again independently, in a second disease, a second
technology (CITE-seq vs. mass cytometry), and a second calibrated agent.** Raw CITE-seq
antibody-derived-tag counts are UMI/read counts, not molecule counts, and cannot be equated
with a QuantiBRITE-derived n\* without an independent bridging calibration this project does
not have the wet-lab access to run. Notably, in several samples (e.g. sample 01: tumor median
10 vs. T-cell median 14) the tumor compartment does not even clearly exceed the normal
compartments on this raw scale — a sign that the raw ADT scale is not tracking real molecule
counts in a biologically meaningful way here, reinforcing that this is a genuine unit mismatch
rather than a borderline-but-usable approximation.

**The consistent finding across two independent placements, in two diseases and two
technologies, is itself the reportable result:** the ceiling (AUC, calibration-free) is the
quantity this project's data can support; the operating point of a specific real drug or CAR
construct is not computable from either dataset without absolute molecule calibration that
neither dataset has. That is a sharper, more defensible statement for having been checked
twice rather than once, and it is the honest reason no third calibrated ROC point is plotted
in this document despite a third real, peer-reviewed, bead-calibrated threshold now being on
file. Code: `bcma_operating_point.py`.

---

## 12. Figure corrections made, and one that needs Dylan's own decision

**Corrected this session** (factual accuracy only — see the ISEF note below):

- `aim12_single_cell_cytof_reordering.py` — the `REAL_SINGLE_CELL_REORDERING_RESULT` dict now
  carries an explicit block recording which counting convention produced 16/27 and 0.818, the
  strict pairwise alternative (13/27, 0.531), and the convention-free CDF-crossing metric. The
  original values are unchanged, so every figure still renders identically; the correction is
  additive, not a silent rewrite.
- **Figure 8** caption: removed the sentence describing 81.8 pp as a between-patient spread;
  it now states which criterion is in use and gives the stricter figure (48.1%, 13/27).
- **Figure 8 poster** banner: "7-11x more target reordering" was derived from the disputed
  convention (it becomes roughly 5-9x under the strict one), so it now reads "several times
  more" — true under either.
- **Figure 6, Figure 5, Figure 0** captions: now cite both criteria explicitly rather than
  presenting 59.3% as unqualified.

**Needs Dylan's decision — Figure 5's hero example is gating-dependent.** Figure 5 and its
poster version are built on the comparison "P7676 (mean 109.1) displaces Q7915 (mean 111.8)
despite Q7915 having the higher average density." Both means are **ungated whole-biopsy**
values. After gating to B cells the two reverse: **Q7915 = 227.0, P7676 = 139.0**. The figure
is still an accurate depiction of the Aim 12 analysis as actually run, but the sentence "even
though Q7915 has a higher average antigen density" is only true of ungated data, and a judge
who asks "were those gated?" will find the example inverts.

Three options, in order of how much work they cost:
1. Add a clause to the caption stating the means are whole-biopsy, ungated.
2. Pick a replacement example whose direction survives gating, from the B-gated re-run.
3. Rebuild Figure 5 on gated data entirely, which is the strongest but requires re-deriving
   the top-5 lists.

**ISEF note.** Society for Science bars generative AI from writing the submitted abstract,
research plan, or poster text. The edits above are factual corrections to numbers and
criteria, not new prose to submit. Every caption that will appear on the poster should be
rewritten by Dylan in his own words before submission.

---

# Second addendum (September 7, 2026): third compartment, gating robustness, and replication in a second disease

## 13. A third normal compartment (step 6)

The gated cache was rebuilt for all 77 files to emit a third population: QC-passing CD45+
leukocytes that are neither B nor T — the myeloid and NK bystanders that were previously
discarded. A target must now spare **three** compartments, so the ceiling is the minimum over
T cells, myeloid/NK, and batch-matched normal B cells.

**The theorem's monotonicity requirement is satisfied exactly.** Computing the two- and
three-compartment ceilings from the *same* AUC estimates: **13 of 37 ceilings fell, 0 rose.**
(Comparing against the earlier saved run instead shows small apparent rises; that is
resampling jitter, measured here at median +0.0003, max 0.0087, which is why the check is done
within a single run.)

**Adding the compartment removed the spurious winners.** The largest falls are exactly the
antigens flagged in section 6 as scoring well only because the compartments that would sink
them were unmeasured:

| antigen | 2-compartment | 3-compartment | change |
|---|---|---|---|
| CD38 | 0.538 | 0.399 | −0.138 |
| CD86 | 0.609 | 0.481 | −0.127 |
| HLA-ABC | 0.546 | 0.459 | −0.087 |
| CD49D | 0.538 | 0.458 | −0.080 |
| CD25 | 0.506 | 0.431 | −0.075 |
| PD-L1 | 0.605 | 0.535 | −0.069 |
| CD80 | 0.616 | 0.548 | −0.068 |

**And the framework now converges on the real drug target.** With three compartments the
highest ceiling in the entire 37-antigen panel is **CD19 at 0.580** — the antigen that actually
became the drug. CD80, the previous top scorer, falls to 0.548. This is a meaningful
strengthening of section 6's caveat: the earlier per-patient "winners" were an artifact of a
narrow compartment set, and adding compartments corrects it in the direction the theory
predicts.

The binding compartment is now recorded per antigen: for CD19, CD20, CD22 and CD79B it is
normal B cells; for CD38, HLA-ABC, CD44, CD86 and PD-L1 it is the new myeloid/NK compartment.

## 14. Gating robustness (step 10)

The whole tumor-versus-normal ceiling was recomputed with two structurally different, fully
automatic gating methods on arcsinh-transformed signal: **Otsu's method** (cutoff maximising
between-group variance, no percentile assumption) and a **two-component Gaussian mixture**
fitted by EM. All three methods were run on the *same* 12 tumors and 6 RLN controls so the
comparison is exact.

| | percentile | Otsu | GMM |
|---|---|---|---|
| CD19 | 0.398 | 0.417 | 0.381 |
| CD20 | 0.543 | 0.660 | 0.610 |
| CD22 | 0.291 | 0.282 | 0.280 |
| CD79B | 0.467 | 0.533 | 0.491 |

Rank agreement: percentile vs Otsu **0.922**, percentile vs GMM **0.956**, Otsu vs GMM
**0.984**.

**The conclusions do not depend on the gating method.** An earlier version of this check
appeared to show large instability in CD19 (0.580 vs 0.415), but that comparison was confounded
— it compared the full 41-tumor cohort against a 12-tumor subset. On matched files the three
methods agree closely. The real sensitivity is **cohort size**, not gating: percentile gating on
12 tumors versus all 41 gives a rank correlation of 0.857.

Under all three methods, on matched files, the four approved-drug targets sit at or below
chance (CD19 0.38–0.42, CD22 0.28–0.29), which strengthens rather than weakens the headline.

## 15. Replication in a second disease: multiple myeloma (step 9)

The framework was re-run end to end on **GSE234261** (Rade, Grieb, Weiss et al., *Nature Cancer*
2024) — a BCMA CAR-T multiple myeloma cohort already held by this project but never analysed
this way. It is a genuinely independent test: **different cancer, different lab, different
technology** (CITE-seq antibody-derived tags rather than mass cytometry), **different targets**.
11 samples, 73,755 cells, **63 surface proteins per cell**, CLR-normalised.

Populations were gated from the panel's own lineage markers: malignant plasma cells
(CD138+ and/or BCMA+, CD3E−CD14−CD19−) against four normal compartments — T cells, myeloid,
NK cells, and normal B cells.

**Results.**

| protein | ceiling | rank of 63 | clinical status | limiting compartment |
|---|---|---|---|---|
| CD138 | 0.698 | 1 | indatuximab ravtansine — **discontinued** | T cells |
| **CD269 (BCMA)** | **0.608** | **2** | **ide-cel, cilta-cel, teclistamab — approved** | T cells |
| FITC | 0.528 | 3 | panel control fluorophore | — |
| CD319 (SLAMF7) | 0.524 | 4 | elotuzumab — approved | normal B cells |
| PE | 0.495 | 14 | panel control fluorophore | — |
| CD307 (FcRH5) | 0.451 | 31 | cevostamab — in trials | — |
| CD38 | 0.287 | 49 | daratumumab, isatuximab — approved | **NK cells** |

Plasma-cell-lineage targets (BCMA, SLAMF7, FcRH5, CD138) rank at mean position **9.5 of 63**,
with a 0.888 probability of outranking a random other protein.

**Four things worth stating precisely.**

1. **The negative controls behave.** FITC and PE are fluorophores in the panel with no
   biological target; they land at 0.528 and 0.495, i.e. chance. Their *ranks* look high only
   because most of the 63 proteins are compressed near 0.5 — which is itself the honest
   reading: **only CD138 (0.698) and BCMA (0.608) clearly separate from the control values.**
   Everything else, SLAMF7 included, is within noise of a fluorophore.
2. **BCMA is the top viable target, out of 63, picked by the method.** It is the flagship
   myeloma CAR-T and bispecific target. Its ceiling (0.608) is comparable to CD19's in lymphoma
   (0.580) — both barely above chance, in two diseases, two technologies, two labs.
3. **The "approved despite a bad window" pattern replicates.** CD38 is an approved myeloma drug
   target with a ceiling of 0.287, limited by NK cells — and daratumumab's well-documented
   on-target effect is exactly NK-cell depletion. As with CD19 and B-cell aplasia in lymphoma,
   the drug was approved not because the target is selective but because the collateral is
   tolerable. **That is the section-4 thesis, reproduced in an independent disease.**
4. **CD138 is a real, named illustration of the framework's boundary.** It tops the ranking on
   bone-marrow data and failed in humans: indatuximab ravtansine was discontinued because
   CD138 (syndecan-1) is heavily expressed on normal **epithelium** — skin, lung, kidney — a
   compartment a bone-marrow sample cannot contain. This converts section 6's caveat from a
   hypothetical into a documented case with a drug name attached: the ceiling is a necessary
   condition against measured compartments, never a sufficient one.

**Not attempted here:** an outcome analysis. The cohort has CR/non-CR labels at n=5 vs 5, which
an earlier session already tested (Mann-Whitney p=0.22 on means, best p=0.0952 across a
threshold sweep). No reframing rescues that sample size; it is reported as directionally
consistent and nothing more.

---

# Third addendum (September 7, 2026): a healthy human reference, and the two-antigen extension

## 16. Testing the thesis directly in healthy people (step 8)

The central claim — *the antigens that became drugs are not tumor-selective, they are restricted
to a lineage the body can do without* — is a claim about **normal** biology. So it can be tested
entirely inside a healthy dataset, with no cross-technology comparison at all.

**Data.** GSE164378 (Hao et al., *Cell* 2021): **161,764 PBMCs from 8 healthy donors, 228 surface
proteins per cell**, with the authors' own expert cell-type annotations (Mono 49,010; CD4 T
41,001; CD8 T 25,469; NK 18,664; B 13,800; other T 6,789; DC 3,589). The compartments are theirs,
not gates of mine. The panel covers the targets of **both** diseases studied here, plus four
isotype-control antibodies and eleven antigens measured with two independent clones.

**Measure.** For each protein, *restriction* = the largest AUC separating any one healthy lineage
from all the others. 1.0 means it marks exactly one lineage; 0.5 means it is spread evenly.

**Controls first.**
- The four isotype antibodies bind nothing and score **0.562–0.598**, not 0.5. **The effective
  noise floor is therefore 0.598**, and any restriction below that is indistinguishable from an
  antibody with no target. This is an essential calibration and is used throughout.
- The eleven duplicate clone pairs agree to a **median of 0.034**. One pair is a real outlier —
  CD275's two clones differ by 0.315 — and is excluded from interpretation.

**Result.** The four approved lymphoma targets are, in healthy humans, among the most
B-lineage-restricted proteins measurable:

| target | restriction | home lineage | rank of 224 |
|---|---|---|---|
| CD19 | **1.000** | B | **1** |
| CD22 | 0.990 | B | 8 |
| CD20 | 0.986 | B | 12 |
| CD79b | 0.886 | B | 56 |

**CD19 is the single most lineage-restricted protein of 224 measured, in healthy people.** Taken
with section 3 — where the same four score at chance against normal B cells — this is the thesis
stated as two measurements rather than an argument: *not tumor-selective, maximally
lineage-selective.* The drugs work because the lineage is expendable.

**Honest limit of this reference.** It is blood, so it can only score targets whose home lineage
is present in blood. Plasma cells live in bone marrow, so every myeloma target scores at or below
the isotype floor here (BCMA 0.580, SLAMF7 0.612, FcRH5 0.577, CD138 0.594) and those numbers are
**uninformative, not evidence against them**. Scoring myeloma targets needs a healthy *marrow*
reference.

### 16b. A forward-looking prediction, and two hits

Filtering all 224 antibodies to those whose home lineage is B (dispensable) and whose restriction
clears the weakest approved target (CD79b, 0.886) gives ten proteins. Four are the known approved
targets. The other six are predictions. Two of them are, on checking the literature afterwards,
exactly the targets the field independently chose:

- **CD268 = BAFF-R (0.990).** Qin et al., *Science Translational Medicine* 11:eaaw9414 (2019),
  doi 10.1126/scitranslmed.aaw9414 — *"CAR T cells targeting BAFF-R can overcome CD19 antigen loss
  in B cell malignancies."* Developed specifically as the CD19-relapse salvage target.
- **CD72 (0.973).** Izgutdina et al., *J Immunother Cancer* 13(12) (2025),
  doi 10.1136/jitc-2025-012013 — affinity-matured CD72 nanobody CAR-T, concluding that
  *"for B-ALL in particular, CD72 may be a preferable second-line immunotherapy target over CD22."*

The method was given healthy-donor protein measurements and nothing else. It ranked CD19 first,
recovered all four approved targets, and surfaced the field's two leading next-generation B-cell
targets. The remaining predictions (CD21, CD196/CCR6, CD275, HLA-DR) are untested; HLA-DR should
be discounted because it is also high on monocytes and dendritic cells.

## 17. The two-antigen extension (the multi-target theorem)

**This section was rewritten on September 7, 2026 after two bugs of mine were caught by the
myeloma replication. The first version overstated the value of combining antigens by roughly
tenfold. Both bugs and their effect are documented below rather than quietly corrected,
because the correction is the more interesting result.**

### The hierarchy

With one antigen, cells carrying the same amount are indistinguishable, which is what forces
the single-antigen ceiling. With two, each cell is a point on a plane, and populations that
overlap completely along either axis alone can still be separated. Four nested families:

| | rule | contains |
|---|---|---|
| 1. SINGLE | kill where one antigen is high **or** where it is low | — |
| 2. AND gate | kill where x > a **and** y > b | (1) |
| 3. NOT gate | kill where x > a **and** y < b | (1) |
| 4. LIKELIHOOD RATIO | kill where P_tumor(x,y)/P_normal(x,y) > c | all of the above |

Each family contains the one above, and family 4 is optimal over *all* rules by the
Neyman–Pearson lemma, so **SINGLE ≤ max(AND, NOT) ≤ LR** is forced for Youden's J
(= sensitivity + specificity − 1).

### Two bugs, both found by replication

**Bug 1 — the single-antigen baseline only tested "kill where the marker is HIGH."** Many real
markers are *low* on tumor (CD3E in myeloma; IgD in lymphoma), for which the correct
one-antigen rule is "kill where it is low." Omitting that made a NOT gate look like a
two-antigen gain when it was really a one-antigen effect with the opposite sign. The baseline
now tries both directions on both antigens.

**Bug 2 — the likelihood-ratio rule is the optimum in the population but is not estimable
from this much data.** It selects an arbitrary subset of 16×16 = 256 histogram bins —
hundreds of effective free parameters against a gate's two thresholds — so it wins in-sample
by construction and generalises worse. Measured overfit gap (train minus held-out):

| family | DLBCL | myeloma |
|---|---|---|
| single antigen | −0.001 | ~0 |
| AND / NOT gate | −0.001 | 0.020 |
| likelihood ratio | **0.020** | **0.261** |

The fix is to report both numbers: the in-sample value verifies the theorem, the held-out
value says what actually works.

### Results

**The theorem is verified exactly, in both diseases, in-sample — where it is a mathematical
claim.** 0 violations of SINGLE ≤ AND/NOT ≤ LR in **0 of 666** DLBCL pairs and **0 of 1,953**
myeloma pairs.

**The practical gain from a second antigen is real but small.** Using the best simple gate
(not the unusable LR estimate) against a both-directions single-antigen baseline, all on
held-out cells:

| | DLBCL (16 patients, 666 pairs) | myeloma (11 samples, 1,953 pairs) |
|---|---|---|
| pairs where two antigens beat one | 560 (84%) | 1,248 (64%) |
| median gain | **+0.015** | **+0.003** |
| maximum gain | +0.163 | +0.136 |
| pairs gaining more than 0.05 | 90 | — |
| NOT gate beats AND gate | 524 (79%) | 1,076 (55%) |

The first version of this section reported a median gain of +0.152 and a best pair of
"CD20 NOT IgD, 0.215 → 0.645." **Both were artifacts of bug 1.** The corrected figures for
that pair are 0.468 → 0.582, a gain of +0.114 — still real, still among the best pairs, but a
fraction of what was claimed.

**Corrected top DLBCL pairs** (held-out gain over the best single antigen in the same pair):
CXCR5 NOT IgD **+0.163** (0.426 → 0.589), CD21+CXCR5 +0.142, CD22+CD21 +0.139,
CD22+CXCR5 +0.139, CD40 NOT IgD +0.116, CD20 NOT IgD +0.114.

**For the approved myeloma target pair, combining buys nothing**: CD138 + BCMA gives a gain of
−0.000.

### What survives, and what it means

1. **The mathematics is verified.** The nested hierarchy holds exactly in-sample in both
   diseases, across 2,619 antigen pairs.
2. **Exclusion beats requirement.** NOT gates outperform AND gates in 79% (DLBCL) and 55%
   (myeloma) of pairs. This survived both bug fixes and replicates across diseases, and it
   runs against the field's emphasis on AND logic (e.g. the Boolean screen of Dannenfelser
   et al., *Cell Systems* 2020).
3. **Combining two antigens buys much less than the theory allows.** The theoretical ceiling
   is far above what any estimable rule achieves at this sample size. That gap — not the
   ceiling itself — is the practically important number, and it is the kind of thing only a
   held-out analysis reveals.
4. **A methodological finding worth stating on its own:** the Neyman–Pearson optimum is the
   right theoretical object and the wrong practical one. A two-threshold gate generalises;
   the fitted likelihood ratio does not. With single-cell data at this scale, chasing the
   optimum is counterproductive.

### Stated limits

- Markers used in the gating are partly circular and their pairs should be discounted: CD3,
  CD19, CD20 in DLBCL; CD3E, CD14, CD33, CD56, CD19, CD138, CD269 in myeloma. The strongest
  DLBCL pairs (CXCR5, IgD, CD21, CD22, CD40) are mostly not gating markers; the top myeloma
  pair (CD3E + CD138) is entirely circular and should be ignored.
- DLBCL covers 16 of 41 patients (compute time); myeloma all 11 samples.
- Normal is the pooled compartments; gains are held-out medians across patients.

## 18. The outcome flip: approved drugs versus failed drugs

Every earlier test compared "approved drug" against "broadly expressed antigen" — but the
second label came from expression biology, which is close to what the model measures, so the
test was partly circular. This section removes that. **Both groups are antigens against which
a real drug was built and given to patients. The only difference is what happened to it.**
Nothing about expression enters the labels.

**Ground truth, verified from the literature this session:**

| antigen | outcome | drug |
|---|---|---|
| CD19 | approved | axi-cel, tisa-cel, liso-cel, blinatumomab |
| CD20 | approved | rituximab, mosunetuzumab, glofitamab, epcoritamab |
| CD22 | approved | inotuzumab ozogamicin |
| CD79b | approved | polatuzumab vedotin |
| CD52 | approved | alemtuzumab (CLL) |
| CD23 | **failed** | lumiliximab — LUCID phase 2/3, n=627, **stopped early for lack of efficacy** (Awan et al., *Br J Haematol* 167:466, doi 10.1111/bjh.13061) |
| CD80 | **failed** | galiximab — CALGB 50402, never approved |
| CD40 | **failed** | dacetuzumab / SGN-40, never approved |
| HLA-DR | **failed** | apolizumab / Hu1D10, never approved |
| CD25 | **failed** | camidanlumab tesirine, never approved |
| CD200 | **failed** | samalizumab, never approved |

CD52 was included deliberately even though it was expected to *hurt* the result — alemtuzumab
is approved yet targets a pan-lymphocyte antigen and causes profound lymphopenia. Dropping an
inconvenient case is exactly the error this test exists to avoid.

**Two predictors compared, by exact permutation test:**

| predictor | approved mean | failed mean | separation | exact one-sided p |
|---|---|---|---|---|
| **B-lineage restriction in healthy donors** | 0.926 | 0.753 | **0.840** | **0.048** |
| tumor-vs-normal ceiling (DLBCL) | 0.494 | 0.458 | 0.650 | 0.278 |

**The predictor that distinguishes drugs that worked from drugs that failed is which normal
lineage the target lives on, measured in healthy people — not the tumor-versus-normal
ceiling.** The ceiling is a valid bound on what is achievable, and it correctly excludes
targets, but on this labelled set it does not predict clinical success. Lineage restriction
does, at p = 0.048.

That is the dispensable-lineage thesis with a clinical-outcome p-value attached, and it
sharpens the project's contribution: the ceiling says what a target *can* do; the lineage says
whether it is *worth building*.

**Honest limits, which matter a lot here.** n = 10, and p = 0.048 is fragile — reclassifying
one drug flips it. The labels are mine, from the literature, and should be independently
checked before this is claimed. The two predictors are computed on slightly different subsets
(CD23 is absent from the healthy panel; CD52 and CD200 from the DLBCL panel). The single
misordered pair is HLA-DR (0.902, failed) above CD79b (0.890, approved); HLA-DR is also
expressed on monocytes and dendritic cells, so a breadth criterion would likely fix it — but
adding one *after* seeing the result would be fitting to the answer, and it should be
pre-registered first. This test is now written into `OSF_Preregistration_DRAFT.md`.

---

# Fourth addendum (September 7, 2026): circularity removed, full cohort, and the label set pushed to its ceiling

Sections 17 and 18 were run again after three further problems were found. All three are
recorded rather than silently corrected.

## 19. Gating circularity removed, and the full 41-patient cohort

**The problem.** In section 17 the comparator pooled T cells, myeloid/NK cells and normal B
cells. But the B gate is "(CD19 high OR CD20 high) AND CD3 low", so CD19 and CD20 separate B
cells from T and myeloid cells *by construction*. Any pair containing them was partly circular.

**The fix, which is also the better science.** The comparator is now **batch-matched RLN normal
B cells only**. Tumor B cells and normal B cells are gated by the identical rule, so the gate
cannot manufacture a difference between them and **no marker in the panel is circular**. It is
also the compartment that sections 3–4 showed is the binding one for every B-lineage target.
The analysis now covers **all 41 patients** (a resumable store replaced the 16-patient subset).

**In myeloma**, the same problem is fixed differently: the seven antigens that define the gates
(CD3E, CD14, CD33, CD56, CD19 for the normal compartments; CD138, CD269 for tumor) are
**excluded from the pair set entirely**, leaving 56 antigens and 1,540 pairs.

## 20. A wrong diagnosis of my own, recorded

While checking the results I found that the best single-antigen J reported by the pair
pipeline (0.618) did not match a direct per-antigen calculation (0.410). I hypothesised a
third bug — that cells falling outside the training histogram's range were being dropped,
biasing the marginals — and fixed it by opening the outer bins to ±infinity.

**The hypothesis was wrong.** The number did not move, because dropped cells were never the
cause. The real explanation is **median-of-maximum versus maximum-of-median**: the pipeline
lets the single-antigen baseline pick whichever of the pair's two antigens is better *for each
patient* and then takes the median across patients, which legitimately exceeds the best
fixed-antigen-for-everyone median. Verified directly: for CD27+FasR the median of the
per-patient maximum is 0.616, against per-antigen medians of 0.384 and 0.412.

Both numbers are correct and they answer different questions. The gates are optimised per
patient in exactly the same way, so the gate-versus-single comparison is internally consistent.
The ±infinity change was kept because it is more correct in principle, but it fixed nothing.

## 21. The two-antigen result, corrected

**Two single-antigen baselines are needed, because they answer different questions.**

- **SINGLE_hi** — kill only where the antigen is HIGH. This is what **one receptor** can
  physically do: a CAR kills cells that *have* its antigen; no receptor fires on the *absence*
  of a protein. This is the honest **engineering** baseline.
- **SINGLE_any** — kill where the antigen is high *or* low. Not buildable with one receptor
  (killing on absence needs an inhibitory receptor plus an activating one, i.e. two), but it is
  the right **statistical control** for "is the gain just an inverted marker?"

| DLBCL, 41 patients, 666 pairs, normal-B comparator | value |
|---|---|
| in-sample theorem violations (SINGLE ≤ AND/NOT ≤ LR) | **0 of 666** |
| best one-receptor rule (kill-high) | J = 0.618 |
| best two-receptor gate | J = **0.646** (gain +0.028) |
| engineering gain, median across pairs | **+0.133** (98% of pairs improve) |
| statistical gain, median across pairs | **+0.003** (70% improve, only 11 pairs > 0.05) |
| NOT gate beats AND gate | **574 of 666 (86%)** |

| myeloma, 11 samples, 1,540 non-circular pairs | value |
|---|---|
| in-sample theorem violations | **0 of 1,540** |
| engineering gain, median | +0.027 (89% improve) |
| statistical gain, median | +0.003 (65% improve) |
| NOT beats AND | **792 of 1,540 (51%)** |
| overfit gap: likelihood ratio vs simple gate | 0.267 vs 0.023 |

**What this means, stated carefully.** The large engineering gain (+0.133) and the tiny
statistical gain (+0.003) are both real and they are not in conflict. Most of the benefit of a
two-receptor construct is **not** that two antigens carry more information than one. It is that
an inhibitory second receptor **unlocks information already present in a single marker** — the
markers that are LOW on tumor — which no activating CAR can use alone. Measured directly, the
best single markers in the "kill-low" direction are IgD (J = 0.552), CD21 (0.517) and CD23
(0.443); in the buildable "kill-high" direction the best is FasR at 0.410.

The concrete design that falls out: **an activating CAR against a B-lineage antigen, paired
with an inhibitory receptor against IgD**, which spares naive B cells. Its ceiling is modest —
the best two-receptor J against normal B cells is 0.646 versus 0.618 for the best single
receptor.

**A finding that did NOT replicate.** "Exclusion beats requirement" is strong in DLBCL (86% of
pairs) but is a coin flip in myeloma (51%). It should be reported as a lymphoma finding, not a
general one. This is a disconfirmation of the generality I claimed in the third addendum.

**A methodological finding that did replicate.** The likelihood-ratio rule is the population
optimum and the wrong practical object: its overfit gap is 0.017 (DLBCL) and 0.267 (myeloma)
against ~0.02 for a two-threshold gate. Chasing the optimum is counterproductive at these
sample sizes.

## 22. The outcome test, generalised and pushed to its ceiling

**Generalisation.** The predictor is no longer "B-lineage restriction" but **restriction to the
lineage the cancer came from** — B cells for B-cell malignancies, myeloid cells for AML. Same
idea, stated so it covers more than one disease.

**Inclusion rule, fixed before computing:** a target is included only if (a) a real drug against
it was given to patients with a documented outcome, and (b) the malignancy has a clear normal
counterpart lineage present in blood. Rule (b) excludes myeloma targets (plasma cells are not
in blood) and Hodgkin targets such as CD30 and CD25. Those exclusions are judgements about the
disease, not about any number.

**A secondary predictor was pre-specified before being computed:** *margin* = restriction to the
disease lineage minus the highest restriction the antigen achieves for any other lineage. It was
specified in advance precisely because inventing it after seeing that it might fix the HLA-DR
miscall would have been fitting to the answer.

**Result (5 approved vs 5 failed, balanced):**

| predictor | approved mean | failed mean | separation | exact p |
|---|---|---|---|---|
| **restriction to the disease lineage** | 0.926 | 0.742 | **0.840** | **0.048** |
| margin over next-best lineage (pre-specified) | 0.298 | 0.045 | 0.800 | 0.075 |

**The pre-specified secondary predictor made things worse, not better.** Had it been added
post hoc it would have looked like a refinement; specified in advance, it is an honest negative.
It does not rescue the HLA-DR miscall either — HLA-DR's margin is 0.002 (correctly low), but
CD40's is 0.345 (high, yet the drug failed).

**Sensitivity analysis, flagged as unbalanced.** Adding CD47 (magrolimab, ENHANCE phase 3
programme discontinued) and CD70 (cusatuzumab), both failures with no matching approvals, gives
separation 0.886 at p = 0.015. This is reported separately and never as the headline, because
adding only failures can move the result in one direction. The imbalance is not cherry-picking:
the panel contains **no other approved blood-cancer surface target** meeting the inclusion rule
(CD23 and CD33 are absent from the panel; CD30 and CD25 are Hodgkin; CD38, BCMA and SLAMF7 are
myeloma).

**The ceiling on this test, stated plainly.** The number of *approved* surface targets in blood
cancer is small in the world, not just in this panel. **n cannot be raised much by more work.**
This test will remain underpowered however carefully it is done, and belongs in the write-up as
suggestive evidence with an exact p-value attached — not as a significant result.

---

## 23. The whole-body extension: the framework leaves blood cancer

Everything through section 22 rested on two blood-cancer datasets. The measure that actually
predicted clinical success — restriction to a normal lineage — needs only a normal reference,
not a matched tumour. So it can be computed for **any** cancer, provided a whole-body normal
reference exists.

**Reference A (primary).** Human Protein Atlas single-cell consensus: **154 human cell types**
spanning every major organ, ~20,000 genes, 3,087,080 rows.
**Reference B (independent replicate).** HPA bulk tissue consensus: 51 tissues.

**Two measures, both parameter-free:**

| measure | definition |
|---|---|
| `top_share` | share of a gene's total normal expression that sits in its single highest cell type |
| `n_types` | number of cell types expressing above 10% of that maximum |

**Panel.** 44 clinical targets as of September 10, 2026 (42 as of September 9; 40 as of
September 8; grown by two on each of the last two days — see the note below). Every label is
frozen in `src/aim13_scoping/target_labels.py`, a file that no analysis script may edit, and
every label is clinical safety history only — nothing about expression enters a label.

- **19 `safe`** — reached approval or documented tolerability with no dose-limiting on-target
  normal-tissue toxicity.
- **11 `toxic`** — documented severe, dose-limiting or fatal on-target toxicity.
- **12 `efficacy_fail`** — the programme failed, but for lack of efficacy. **Excluded from the
  toxicity test** and used only as a specificity control (section 23.4).
- **2 excluded before any number was computed**: PTPRC (no cytotoxic agent ever reached
  patients; retained only as a restriction negative control) and IL3RA (one approval *and* one
  toxicity failure on the same gene — ambiguous by our own rule).

**Panel growth, September 9, 2026.** The panel was checked against the full 72-gene ADC trial
registry (De Schepper et al. 2025, Supplementary Data 1/2 — the same registry already used as
an independent specificity control in section 27.3) for genes with a settled clinical verdict
not yet recorded here. Roughly 50 candidate genes were checked individually; the great majority
are still in active, undecided trials (which is exactly why most of them already appear in
`prospective_predictions.json` as forward predictions rather than as backward-fit labels).
Two had a clean, verifiable, settled verdict and were added:

- **F3 (tissue factor), toxic.** Tisotumab vedotin (Tivdak) — FDA accelerated approval 2021,
  full approval 2024, cervical cancer. The FDA label carries a **boxed warning for ocular
  toxicity**: up to 60% of patients affected (conjunctivitis, dry eye, corneal
  damage/ulceration, blepharitis), dose-limiting, requiring a dedicated eye-care REMS program.
  Same pattern as FOLR1/mirvetuximab — approved despite the toxicity, not because it is absent.
  **Worth stating plainly: the model's own computed home cell type for F3, from RNA data alone
  and with no toxicity information given to it, is "ocular epithelial cells"** — it
  independently locates the exact mechanism behind a documented boxed warning.
- **GPNMB, efficacy_fail.** Glembatumumab vedotin (CDX-011) — the METRIC phase 2b randomized
  trial (gpNMB-overexpressing triple-negative breast cancer) explicitly did not meet its
  primary endpoint (PFS). No on-target toxicity signal was reported; this is a clean efficacy
  failure, used only as a specificity-control addition (section 23.4), not in the toxicity test.

**Panel growth, September 10, 2026.** The same registry check was repeated, focused on the
remaining ~50 unlabelled genes (now including F3 and GPNMB as labelled), prioritising the
targets with the most registered trials on the theory that a mature programme is the most
likely to have reached a settled verdict. About 30 candidates were checked individually against
primary or near-primary sources (FDA filing news, sponsor press releases, peer-reviewed trial
reports); the same pattern held again — the large majority remain in active, undecided trials.
Two more had a clean, verifiable, settled verdict:

- **ERBB3 (HER3), efficacy_fail.** Patritumab deruxtecan's FDA Biologics License Application
  (EGFR-mutated NSCLC) was voluntarily withdrawn in 2026 after the confirmatory phase 3
  HERTHENA-Lung02 trial missed its overall-survival endpoint, despite an earlier positive
  progression-free-survival signal. Daiichi Sankyo and Merck's own public statements attribute
  the withdrawal to the survival miss, not to a safety finding. Two interstitial-lung-disease
  deaths occurred in the treatment arm but were not cited as the reason for withdrawal. Same
  pattern as GPNMB: a clean efficacy failure by the sponsors' own account.
- **AXL, efficacy_fail.** Genmab discontinued enapotamab vedotin in 2020 after its early-phase
  expansion cohorts; the company's own statement is that it "did not meet the stringent criteria
  for proof-of-concept," with "some evidence of clinical activity" that could not be improved by
  alternate dosing or biomarker selection. Explicitly an efficacy-driven decision; no toxicity
  concern was cited.

**A wrinkle on an existing label, documented rather than acted on.** SLC34A2 (NaPi2b) was
already labelled `efficacy_fail` in the panel (lifastuzumab vedotin). A second, more recent
NaPi2b agent — upifitamab rilsodotin — was checked this session and found to have a real,
serious complication: an FDA partial clinical hold after 5 grade-5 bleeding events (4
treatment-related) in roughly 560 patients. But the program's own stated ending, when it
happened, was efficacy (the UPLIFT trial missed its pre-specified response-rate threshold,
15.6% ORR), and no source found ever resolved *why* the bleeding occurred — every account
states the cause was "under evaluation" and moves on to the efficacy result without returning
to it. The label is kept as `efficacy_fail` rather than relabelled `toxic`, because that is
what the evidence actually supports (a confirmed efficacy failure; an unresolved, uncategorised
safety signal), but the ambiguity is recorded in `target_labels.py`'s evidence field so a future
reviewer does not read this label as a clean safety record.

Both ERBB3 and AXL were verified against sponsor statements, FDA filing news, or peer-reviewed
trial reports directly before being added, the same standard used for F3 and GPNMB. The panel
is not yet at the registry's full 72 genes — most of the remaining candidates simply do not have
a settled verdict yet. Growing it further means periodically re-checking those same genes for a
result, not inventing new ones. See `docs/Label_Review_Sheet.md` for the current label sheet
with citations.

### 23.1 Confinement alone separates safe from toxic

| predictor | separation | permutation p | n |
|---|---|---|---|
| `top_share` | **0.880** | **0.00014** | 19 safe vs 11 toxic |
| `n_types` | **0.876** | **0.00018** | 19 safe vs 11 toxic |

(Prior panel, 40 targets, Sept 8: 0.868, p = 0.00031/0.00033, 19 vs 10 — see the panel-growth
notes above. Unchanged by the September 10 additions, since both ERBB3 and AXL are
`efficacy_fail` and this test uses only the safe/toxic groups.)

This is a stronger result on more cases than the blood-only outcome test of section 22
(0.840, p = 0.048, n = 5 vs 5), and it is the first result in this project that covers solid
tumours at all. Representative values: CD20 0.850 / 1 cell type (approved) → BCMA 0.789 / 1
(approved) → PSMA 0.568 / 2 (approved) → CLDN18.2 0.359 / 8 (approved 2024) →
CEACAM5 0.283 / 7 (severe colitis) → F3 0.075 / 31 (boxed-warning ocular toxicity) →
EGFR 0.051 / 60 → HER2 0.041 / 72 (fatal CAR-T event) → CD47 0.030 / 120 (on-target anaemia).

### 23.2 The case that breaks pure confinement, and the rule that survives it

CEACAM5 is confined — 7 cell types, essentially colonocytes — and CEA CAR-T still caused severe
colitis. A confinement-only account cannot explain that, and it should not be patched over.

The project's own thesis already contains the missing factor, and it can be stated as a rule
with **no fitted parameter**:

> **A target is safe only if it is BOTH confined to few normal cell types AND those cell types
> are expendable or regenerating.**

Operationally, "expendable" is assigned from the **identity of the home cell type alone** — a
haematopoietic home means the lineage is continuously regenerated from marrow stem cells *and*
its loss has an established clinical replacement (immunoglobulin, transfusion, growth factor).
No drug outcome enters that assignment, so the test is not circular. The scoring rule is
lexicographic — expendable first, confinement second — implemented as `top_share + 1.0` when
the home is haematopoietic, since the `+1` cannot be beaten by any confinement value.

| rule | separation | permutation p |
|---|---|---|
| one factor (confinement only) | 0.880 | 0.00014 |
| **two factors (confined AND expendable)** | **0.914** | **0.00003** |

(Prior panel, Sept 8: 0.868/0.00031 one factor, 0.905/0.00004 two factors. Unchanged by the
September 10 additions, same reason as 23.1.)

The rule explains every case in the panel including the exceptions:
CD19/CD20 → B cells (confined + expendable) → approved.
BCMA → plasma cells (confined + expendable) → approved.
CLDN18.2 → gastric mucosa (confined + regenerating) → approved 2024.
CEACAM5 → colonocytes (confined, **not** expendable) → colitis.
HER2 / CD47 → everywhere (not confined) → catastrophic.

**A methodological bug found and fixed here, recorded.** The first implementation tested
haematopoietic membership by plain substring match, which classified *"conjunctival goblet
cells"* as haematopoietic because the string contains `"t cell"`. Matching is now
word-boundary anchored and `target_labels.py` carries an assertion-based guard table that
fails at import if the bug ever returns. The corrected two-factor separation is 0.905; the
buggy one read 0.916. **The bug inflated the result, and the number reported here is the
corrected, smaller one.**

### 23.3 Leave-one-out

Dropping each of the 30 outcome-bearing targets in turn (29 as of Sept 8; 30 since F3 was
added) and re-running the two-factor test: separation min **0.905**, median 0.909, max
0.947; p < 0.05 in **30 of 30** drops (p max 0.00008). No single target carries the result.
(Prior panel, Sept 8: min 0.895, max 0.942, p max 0.00022, 29 of 29 drops. Unaffected by the
September 10 additions, same reason as 23.1/23.2 — ERBB3 and AXL are `efficacy_fail`, outside
this safe/toxic test.)

### 23.4 Specificity: the measure is not just "was tried and did not work"

If confinement merely tracked "this programme failed", the efficacy failures would look
like the toxic group. Mostly, they do not — but this test has moved further toward significance
with the September 10 additions, and that movement is discussed honestly below rather than
smoothed over.

| group | mean `top_share` |
|---|---|
| safe (n = 19) | 0.370 |
| efficacy failures (n = 12) | 0.219 |
| on-target toxic (n = 11) | 0.121 |

safe vs efficacy-fail: separation **0.746**, p = **0.011** — clearly significant, having moved
from borderline (0.700, p = 0.043 on the 42-target panel; 0.684, p = 0.063 on the original
40-target panel). This moved because ERBB3 (top_share 0.083) and AXL (top_share 0.102), the two
targets added September 10, both sit well below the efficacy-fail group's prior mean.

**Why this happened, checked rather than waved away.** The efficacy-fail group is not
uniformly low — it splits into two real biological clusters, visible by simply sorting the
group's `top_share` values. Two members, FCER2/CD23 (0.688) and HLA-DRA (0.411), are
blood/immune-lineage antigens, confined and high-scoring, sitting alongside the safe group's
typical range. The remaining ten — CD40, ERBB3, MUC1, AXL, CD80, MUC16, SDC1, GPNMB, CD74,
SLC34A2 — are broadly-expressed **solid-tumour** surface antigens with low confinement, sitting
closer to the toxic group's range (EGFR, ERBB2, CD44 and EpCAM, the toxic group's own broadly-
expressed members, occupy the same low band). This split predates the September 10 additions —
the 10-target efficacy-fail set already had a median `top_share` of 0.183 before ERBB3 and AXL
were added — and both new targets happen to fall on the low side of that existing split rather
than creating a new one.

**The honest interpretation, stated as a limitation of the specificity test rather than of the
underlying thesis.** A broadly-expressed, low-confinement antigen may simply be hard to hit
with any wide therapeutic window, whatever the eventual clinical readout. Whether that shows up
in the record as on-target toxicity or as an efficacy failure may be two faces of one
underlying constraint — for instance, a dose capped for tolerability before a trial could ever
show efficacy would read, in the literature, as "failed for lack of efficacy," even though
confinement was the real limiting factor underneath. That is not a contradiction of the
project's central thesis (if anything it is a plausible corollary of it), but it does mean
this specific control — designed to show confinement is uninformative about efficacy failure in
general — no longer cleanly passes. **What still holds, and is the reason this is reported as a
limitation of one control rather than a threat to the headline result:** the efficacy-fail
group's mean (0.219) sits between the safe group's (0.370) and the toxic group's (0.121) rather
than merging into either one, and every headline safe-vs-toxic result in 23.1–23.3 is completely
unaffected by this section, because efficacy-fail targets never enter those tests. This section
exists to be read honestly, not to be read as undermining sections 23.1–23.3.

### 23.5 Sensitivity: germ and placental cell types removed

Several targets have their single highest cell type in germ or placental tissue
(CA9 → spermatids, DLL3 → oocytes, EGFR → syncytiotrophoblasts, CD52 → epididymal cells).
Those sites are immune-privileged and are not where the clinical toxicity fell. Removing all
20 germ/placental/reproductive cell types and recomputing from scratch:

| rule | with germ tissue | germ tissue removed |
|---|---|---|
| one factor | 0.880 | 0.856 (p = 0.00049) |
| two factors | 0.914 | **0.914 (p = 0.00005)** |

(Prior panel, Sept 8: 0.868 → 0.842/p=0.0011 one factor; 0.905 → 0.905/p=0.00006 two factors.
Unaffected by the September 10 additions, same reason as 23.1–23.3.)

Five targets change home cell type, and the changes are corrective rather than damaging —
CD52's home moves from epididymal principal cells to **B cells**, which is where alemtuzumab
actually acts. The two-factor result is unchanged.

### 23.6 The load-bearing check: solid tumours alone

The obvious objection is that the whole thing works only because blood targets happen to be
both confined and safe. So: **remove every blood-cancer target and re-run.**

- solid `safe` (n = 8): CLDN18.2, DLL3, PSMA, glypican-3, c-MET, mesothelin, Nectin-4, TROP2
- solid `toxic` (n = 10): CAIX, CD44v6, CEA, EGFR, EpCAM, HER2, tissue factor (F3), FAP, FRα,
  CD56

| predictor | separation | exact p |
|---|---|---|
| `top_share` | **0.812** | **0.0133** |
| `n_types` | 0.775 | 0.0253 |

(Prior panel, Sept 8, n = 8 safe vs 9 toxic: 0.792/p=0.0232 top_share, 0.764/p=0.0346 n_types.
Unaffected by the September 10 additions, same reason as 23.1–23.3.)

**It survives, and strengthened.** The finding is not a blood-cancer artifact — adding F3, a
solid-tissue toxic target whose own toxicity mechanism (ocular epithelium) the model correctly
locates from RNA data alone, moved this check in the stronger direction. (No haematopoietic
target remains in this subset, so the two-factor rule reduces to the one-factor rule by
construction and is not reported separately.)

### 23.7 An independent normal reference

Reference B (51 bulk tissues) is a genuinely separate reference, and it replicates —
`top_share`/two-factor separation 0.751, p = 0.0118 (was 0.726, p = 0.024 on the 40-target
panel). It is **weaker by construction**: 51 tissues cannot resolve a B cell from the lymph
node it sits in, which is precisely why the cell-type reference is primary. The haematopoietic
flag never fires on a tissue name, so reference B tests confinement only. Agreement between the
two references: Spearman ρ = 0.592 (p = 0.00002, n = 44) — correlated, but far from identical,
which is the honest statement that the two references are not interchangeable. (Prior panel,
Sept 8: ρ = 0.512, p = 0.0007, n = 40; Sept 9: ρ = 0.550, p = 0.00016, n = 42.)

### 23.8 Limitations of section 23, stated plainly

1. **This is RNA, not protein.** The field's own comparable tools (Dannenfelser, SCAN-ACT) use
   the same reference, but protein would be better. A protein-level cross-check is the named
   next step.
2. **`top_share` depends on how the reference partitions cell types.** A reference with finer
   annotation would push every value down. Only the *ordering* is used, and the ordering
   replicates in an independently partitioned reference (23.7).
3. **The safe/toxic labels are literature judgements made by one person.** They are frozen and
   citable, and pre-registration plus outside review is the fix (see the pre-registration doc).
4. **Modality is collapsed.** HER2 is fatal as a CAR-T and approved as a naked antibody; the
   label records the on-target toxicity, and the framework does not model modality.
5. **The specificity control in 23.4 is now a genuine, disclosed limitation, not a clean pass.**
   See 23.4 for the full honest account; it does not affect 23.1–23.3.

**Scripts:** `src/aim13_scoping/target_labels.py` (frozen labels),
`wholebody_v2.py` (all tests), `wholebody_v2.json` (all numbers),
`wholebody_v2_output.txt` (verbatim run log).
**Figure:** `fig9_wholebody_restriction.png` / `.svg`.

---

## 24. The pancreatic thread, computed end to end

Pancreatic cancer is not a late addition to this project. It is where the project started.
This section is the one place the whole thread is recomputed in a single run
(`src/aim13_scoping/pancreatic_thread.py`, writing `pancreatic_thread.json`).

### 24.1 The calibrations were already pancreatic

| calibration | data | result |
|---|---|---|
| **CLDN18.2** | 8 points across 3 pancreatic + 3 gastric lines (Barrett 2024, Fig 5D) | Spearman r = **0.964**, p = 0.00047, RMSE 0.237; fitted k_intra = 2.99e-05, n_half = 23,177 |
| **mesothelin** | 4 points (Chu 2024) | fitted k_intra = 1.31e-05, n_half = 53,049 (source paper's own r = 0.623, p = 0.035) |

Both are real molecules-per-cell dose-response curves calibrated with antibody-binding-capacity
beads, and both were fitted months before the normal-tissue atlas existed in this project.

### 24.2 CPTAC pancreatic ductal adenocarcinoma — an honest negative

n = 135 patients, 76 deaths, median follow-up 399 days. Splitting at the median mesothelin
protein level (24.314 log2): high-MSLN 67 patients / 39 deaths, low-MSLN 68 / 37.

**Log-rank χ² = 0.873, p = 0.350. The curves do not separate.**

This is reported as a negative result, and it is the right kind of negative to report: it is a
**prognostic** question (does the target's level track survival?), not a **response** question
(does the target's level predict whether a drug works?), and the framework never claimed the
former. Reporting it protects the project from the reasonable objection that only wins are shown.

### 24.3 Whole-body ranking of 14 real PDAC surface targets

Targets were chosen by "has been taken into a PDAC trial or is a named PDAC CAR-T/ADC target",
fixed before any confinement value was looked at.

| rank | target | confinement | cell types | clinical status |
|---|---|---|---|---|
| 1 | **claudin-18.2** | **0.359** | 8 | zolbetuximab **approved 2024**; CT041 CAR-T in PDAC |
| 2 | CEA (CEACAM5) | 0.283 | 7 | severe on-target colitis |
| 3 | PSCA | 0.270 | 8 | in trials |
| 4 | TROP2 | 0.258 | 7 | approved in other tumours |
| 5 | **mesothelin** | **0.239** | 11 | the long-standing tolerable PDAC target |
| 6 | CEACAM6 | 0.236 | 9 | preclinical |
| 7 | CDH17 | 0.230 | 9 | early |
| 8 | CA125 / MUC16 | 0.173 | 13 | failed |
| 9 | B7-H3 | 0.124 | 15 | in trials |
| 10 | integrin β6 | 0.098 | 23 | early |
| 11 | MUC1 | 0.093 | 25 | repeatedly failed |
| 12 | glypican-1 | 0.078 | 44 | diagnostic EV marker |
| 13 | EGFR | 0.051 | 60 | on-target skin/lung toxicity |
| 14 | HER2 | 0.041 | 72 | fatal on-target CAR-T event |

**A CORRECTION TO AN EARLIER, SMALLER RUN — recorded rather than quietly replaced.** An
earlier pass ranked only six PDAC targets and produced the claim *"the two top-ranked
pancreatic targets are the two this project independently calibrated."* On the fuller
14-target panel that claim is **false**: CLDN18.2 still ranks first, but mesothelin ranks
**fifth**, with PSCA, TROP2 and CEA above it. The smaller panel had made the concordance look
tighter than it is. The correct, weaker statement is the one used from here on:

> The single PDAC target that reached approval ranks **first of fourteen**; the three PDAC
> programmes that failed or caused on-target toxicity (MUC1, EGFR, HER2) rank **11th, 13th and
> 14th**; mesothelin, the field's tolerable workhorse, sits in the upper third.

**This is a concordance check, not a statistical test.** There are too few PDAC targets with a
settled verdict to attach a p-value, and doing so would be dishonest. What is checkable is that
the ranking was produced by a normal-tissue atlas containing **no pancreatic tumour data and no
clinical outcome data whatsoever**, and that CEACAM5's high placement is the framework's own
documented exception — confined, but to colon, which is not expendable.

**An incidental finding worth keeping.** Glypican-1, the pancreatic-cancer *exosome* marker that
made GPC1 famous, ranks 12th of 14 on confinement (44 normal cell types). A good diagnostic
marker is not automatically a good therapeutic target — the two jobs have different
requirements, and this framework says exactly why.

### 24.4 Why the thread is not retrofitted

The project's *original* headline — the CA19-9 versus EGFR ranking flip, the first demonstration
that a naive mean orders targets wrongly once the distribution is accounted for — was measured on
PANC-1 pancreatic tumour cells against HPDEC normal pancreatic cells. GPC1, the pancreatic
exosome marker, is the literal join between the EV foundation (Aims 1–9) and the ceiling work
(Aims 10–15). The calibrations are pancreatic. The survival data is pancreatic. The thread runs
from the first aim to the last.

**Figure:** `fig10_pancreatic_thread.png` / `.svg` (calibration, the negative survival result,
and the ranking).

---

## 25. A shuffled-label null for the two-antigen machinery

Section 21 reports two numbers for what a second antigen buys: a median **statistical** gain of
+0.003 over the best single antigen, and a median **engineering** gain of +0.133 (the part no
single activating receptor can physically reach, because a CAR cannot fire on the *absence* of
a protein). Both are modest-looking, and the right question is: *what would this pipeline
report if there were nothing at all to find?*

That has an exact answer. The pipeline is run unchanged, on the same cells, with one
difference — **the tumour/normal labels are shuffled.** Every biological difference is
destroyed and nothing else is. Whatever gain survives is manufactured by fitting alone, and it
is the floor beneath which a real gain cannot be distinguished from noise. This is a stronger
guard than the train/test split already in place, because it also catches bias hidden in the
gating, the subsampling and the binning.

**Run: 6 DLBCL patients, 200 randomly chosen antigen pairs per run, 2 shuffles per patient.**

| quantity | real data | shuffled labels |
|---|---|---|
| best single-antigen J | **+0.692** (0.203 … 0.893) | **+0.002** (−0.015 … +0.011) |
| statistical gain (2 antigens vs the best 1) | **+0.0000** (0.0000 … 0.0023) | +0.0003 (0.0000 … 0.0024) |
| engineering gain (2 antigens vs the best buildable 1) | **+0.106** (0.039 … 0.288) | +0.0008 (0.0000 … 0.0025) |

**The control validates itself first.** With labels shuffled the single-antigen J collapses from
0.692 to 0.002. If it had not, the control would be broken and nothing below it would mean
anything.

**Two conclusions, one of them negative, and the negative is the more useful one.**

1. **The statistical gain is at the noise floor.** Real median +0.0000 against a shuffled floor
   of +0.0003. Combining two antigens adds essentially **no new information** about which cells
   are tumour. Section 21 already said this with the small median and the 11-of-666 count; this
   states it in the strongest available form. Anyone claiming that logic-gated CAR-T works by
   *finding signal that one antigen misses* is, in this dataset, wrong.
2. **The engineering gain is 130× the floor.** Real median +0.106 against +0.0008. This is a
   real and large effect, and it is a different effect: the second receptor does not find new
   information, it **unlocks information already sitting in a single marker that no activating
   receptor can use alone** — the markers that are LOW on tumour. The concrete design that
   follows is an activating CAR against a B-lineage antigen paired with an *inhibitory*
   receptor against IgD, sparing naive B cells.

These two statements are not in tension, and the null control is what makes it possible to say
both with confidence rather than hedging between them.

**Script:** `src/aim13_scoping/two_antigen_null_control.py`; results in
`two_antigen_null_control.json`.

**A documentation correction made in the same file.** `two_antigen_ceiling.py`'s `edges()`
function carried a comment calling the ±∞ outer bins "the third bug" and claiming the fix moved
the best single-antigen J from 0.618 to 0.410. That diagnosis was wrong — the change is correct
and is kept, but it never moved the number. The real explanation is median-of-maximum versus
maximum-of-median, and the comment now says so. This mirrors the record already in §20; the
code and the prose now agree.

---

## 26. The protein-level cross-check: it half replicates, and the failure is the interesting half

Section 23's stated limitation was that the reference is RNA. Two Human Protein Atlas **protein**
references were added to test it, and the pass/fail line was registered **before** they were run
(pre-registration addendum, H5): *separation > 0.70 and p < 0.05, or the RNA result is reported
as reference-dependent.*

### 26.1 The registered result

| reference | groups | separation | p | H5 |
|---|---|---|---|---|
| **PROTEIN, immunohistochemistry** (tissue × cell type) | 76 | **0.785** | **0.0104** | **MET** |
| PROTEIN, IHC, uncertain antibodies dropped (sensitivity) | 76 | 0.775 | 0.0158 | MET |
| **PROTEIN, mass spectrometry**, 60 samples / 20 tissues | 20 | **0.650** | **0.1148** | **NOT MET** |
| PROTEIN, mass spectrometry, summarised to 20 tissues | 20 | 0.650 | 0.1148 | NOT MET |

So the finding **replicates at protein level on immunohistochemistry and does not replicate on
quantitative mass spectrometry.** Both halves are reported, and the second half is not buried.

Placed alongside the RNA references, the four results looked like a clean story about resolution:

| reference | groups | separation |
|---|---|---|
| RNA, single-cell consensus | 154 | 0.868 |
| PROTEIN, immunohistochemistry | 76 | 0.785 |
| RNA, bulk tissue consensus | 51 | 0.726 |
| PROTEIN, mass spectrometry | 20 | 0.650 |

Monotone in granularity. The tempting conclusion is *"the MS reference is simply too coarse."*
That is exactly the kind of after-the-fact excuse a pre-registration exists to prevent, so it
was **tested rather than asserted.**

### 26.2 The excuse was tested, and it failed

Two candidate explanations were checked. Both are recorded, including the one that looked
promising and did not survive.

**Candidate 1 — the MS table is sparse.** The Intensity column is empty for most gene × tissue
cells (39.5% coverage), and two targets were measured in exactly one tissue, which forces
`top_share` to 1.000 by arithmetic rather than biology. Restricting to targets measured in at
least K tissues, post hoc and labelled as such:

| K | n | separation |
|---|---|---|
| ≥ 1 | 14 vs 10 | 0.650 |
| ≥ 3 | 10 vs 9 | 0.600 |
| ≥ 5 | 6 vs 7 | 0.548 |
| ≥ 8 | 4 vs 7 | 0.393 |

Separation **falls** as the sparsely-measured targets are removed. Sparsity is not the
explanation.

**Candidate 2 — 20 groups is too coarse for the measure.** Testable directly: take the
154-cell-type RNA reference, randomly merge the cell types down to K groups, sum within groups,
and recompute. 60 random partitions per K.

| K | separation, mean [min, max] |
|---|---|
| 154 | 0.868 |
| 76 | 0.863 [0.800, 0.889] |
| 51 | 0.854 [0.795, 0.900] |
| **20** | **0.837 [0.753, 0.911]** |
| 10 | 0.831 [0.737, 0.911] |
| 6 | 0.813 [0.721, 0.953] |

**The measure barely degrades with coarsening.** At K = 20 the RNA reference still gives 0.837,
nowhere near the 0.650 mass spectrometry produced, and 0.650 sits outside the range of all 60
random partitions. **Granularity does not explain the failure either.** Both excuses are dead.

### 26.3 What the references actually agree about

If the four references were four measurements of one quantity, they would correlate strongly.
They do not.

| pair | Spearman ρ | p | n |
|---|---|---|---|
| MS samples vs MS tissues (same data, two summaries) | **0.963** | <0.0001 | 34 |
| IHC vs IHC strict | 1.000 | — | 33 |
| MS vs IHC | 0.615–0.702 | ≤0.0002 | 30–32 |
| RNA bulk vs IHC | 0.494 | 0.0026 | 35 |
| **RNA single-cell vs MS** | **0.346** | 0.045 | 34 |
| RNA single-cell vs IHC | 0.429 | 0.010 | 35 |
| RNA bulk vs MS | 0.178 | 0.31 | 34 |

The two protein references agree with **each other** (ρ ≈ 0.62–0.70) considerably better than
either agrees with RNA (ρ ≈ 0.35–0.43), and the MS internal-consistency check passes at 0.963,
so the disagreement is not noise in the MS table.

### 26.4 A WRONG TEST OF MY OWN, and the right one

The random-merging test in 26.2 was **the wrong test**, and it produced the wrong conclusion.
It is kept above rather than deleted, because the mistake is the instructive part.

**Why it was wrong.** Random merging preserves distinctness. A gene expressed in exactly one of
154 cell types still lands in exactly one random bucket, so its confinement stays near 1.0 no
matter how few buckets there are — which is why separation barely moved (0.837 at K = 20). Real
anatomical pooling does not behave that way at all. A B-cell gene is spread across every organ
containing lymphocytes — lymph node, spleen, tonsil, gut, lung, bone marrow — so pooling by
ORGAN dilutes it in a way random pooling never does. The test answered *"does having fewer
buckets hurt?"* (no). It never asked *"does mixing cell types into real organs hurt?"*, which is
the actual question.

**The right test.** The mass-spectrometry panel covers 20 named human tissues, and 18 of those
names appear verbatim in the HPA bulk RNA reference. So the identical measure can be computed on
the SAME organs, in RNA and in protein, with everything else held constant.

| reference | resolution | separation | p |
|---|---|---|---|
| RNA, single-cell consensus | **154 cell types** | **0.868** | 0.0003 |
| PROTEIN, antibody staining | **tissue × cell type (~96)** | **0.785** | 0.0104 |
| RNA, bulk tissue consensus | 51 whole organs | 0.726 | 0.0244 |
| **RNA, the same 18 organs the MS panel covers** | **18 whole organs** | **0.626** | **0.1425** |
| **PROTEIN, mass spectrometry** | **20 whole organs** | **0.650** | **0.1148** |

**RNA fails on those organs too, and by the same amount** — 0.626 against mass spectrometry's
0.650, a difference of 0.024. The two technologies agree once anatomical resolution is matched.

### 26.5 The conclusion, corrected

The pattern is a clean 2 × 2, and technology is not the variable:

| | cell-type resolution | whole-organ resolution |
|---|---|---|
| **RNA** | 0.868 ✓ | 0.626 ✗ |
| **PROTEIN** | 0.785 ✓ | 0.650 ✗ |

**The section 23 result is not RNA-specific. It is resolution-specific**, and it replicates at
protein level whenever the protein reference resolves cell types.

This is a *sharper* claim than section 23 alone stated, and it is one the framework itself
predicts. The whole thesis is that safety depends on **which cell type** carries a target and
whether that cell type is expendable. An organ homogenate deliberately destroys exactly that
information: "lymph node" contains B cells, T cells, endothelium and stroma at once, so a
B-cell-restricted target and a broadly expressed one look similar there. The measure failing at
organ resolution is therefore not a defect — **it is a prediction of the theory, and it is now a
confirmed one.**

The scope condition stands, but it is a statement about what to measure rather than a retreat:

> Restriction must be measured across normal **cell types**. Measured that way it separates
> tolerated targets from on-target-toxic ones in RNA (0.868) and in protein (0.785). Measured
> across whole organs it does not, in either technology (0.626 RNA, 0.650 protein) — because
> whole-organ pooling destroys the cell-type distinction the framework depends on.

**H5: MET at cell-type resolution, in both technologies. NOT MET at whole-organ resolution, in
both technologies.** The registered NOT-MET result on mass spectrometry stands as recorded; what
26.4 establishes is that it is not a protein-versus-RNA failure.

**H6 was NOT TESTABLE** on either protein reference: both are indexed by tissue or tissue × cell
type in a way that never marks a home as haematopoietic, so the flag fires on 0 of 24 labelled
targets and the two-factor rule collapses onto the one-factor rule by construction. Reported as
untestable, not as a failure.

**Also recorded: a NaN bug in the first version of `matched_resolution_check.py`.** The MS table
leaves Intensity blank where a protein was not detected, and averaging an empty list produced
NaN, which made the mass-spectrometry row read 0.000 instead of 0.650. Found and fixed before
any conclusion was drawn from it.

**Scripts:** `protein_cross_check.py`, `ms_coverage_diagnostic.py`, `resolution_sensitivity.py`
(the wrong test, kept), `matched_resolution_check.py` (the right one), `reference_agreement.py`.
Results in the matching `.json` files.

---

## 27. The ceiling in a solid tumour, against matched normal tissue from the same patients

Sections 1–22 computed the tumour-versus-normal ROC ceiling in blood cancers, because blood is
where a matched normal counterpart is easy to obtain. Section 22's most counterintuitive finding
was that this ceiling **does not** predict clinical success — approved targets sit near chance
against their own normal counterpart (CD19 0.576, CD20 0.554, CD79B 0.489, CD22 0.359) — and
that the measure which *does* predict success is restriction to a normal lineage.

That is the claim most likely to be disbelieved, so it was transplanted to a completely
different setting.

**Dataset.** De Schepper *et al.*, *Nature Communications* 2025 (s41467-025-67840-z): a
**rapid-autopsy** study of **30 patients with metastatic breast cancer**, 749 metastatic lesions
across 36 tissue sites, **plus normal tissue sampled from the same autopsies**, with **72
antibody-drug-conjugate target genes** measured by RNA-seq in every sample. Tumour and normal
come from the same body, the same protocol, the same batch — the confound that sections 19–21
had to control for by careful identical gating is removed here by the study design itself.

Different tumour type, different laboratory, different technology, different country, and a
matched-normal design this project could not have built.

**Scope note, stated before any number.** These are **bulk samples**, not single cells. The
separation computed here is between tumour *samples* and normal *samples*, so it is an upper
bound on the cellular ceiling and is **not** directly comparable to the CyTOF numbers in
sections 1–22. It is reported as its own quantity.

### 27.1 The ceiling does not predict safety here either

Per gene, separation = P(a random metastatic sample carries more than a random normal sample),
the same rank statistic used throughout. (Numbers updated September 29, 2026 for the 44-target panel: one target joined the toxic group as the panel grew, 17 → 18 labelled targets; the earlier values were 0.514/0.528/0.750 and 0.495 for the registry control.)

| | separation | exact p | n |
|---|---|---|---|
| pooled across all samples | 0.543 | 0.398 | 9 safe vs 9 toxic |
| computed **within each patient**, then averaged (16 patients) | 0.556 | 0.365 | 9 vs 9 |

Chance. And the scramble is instructive rather than merely null: TROP2 (approved) 0.781 and
Nectin-4 (approved) 0.724 sit high — but so do HER2 (fatal on-target event) 0.640 and CAIX
(dose-limiting liver toxicity) 0.632, while CLDN18.2 (approved 2024) sits at 0.262 and c-MET
(approved 2025) at 0.223 (within-patient values). Being easy to tell apart from normal tissue has nothing to do with
whether the drug was tolerated.

### 27.2 Head to head, on the same eighteen targets

The cleanest statement the project can make, because nothing differs except which quantity is
measured:

| measure | separation | exact p |
|---|---|---|
| tumour-vs-normal ceiling (this dataset, within-patient) | 0.556 | 0.365 |
| **confinement across normal cell types** | **0.778** | **0.025** |

Same targets. Same labels. Same statistic. One measure is uninformative and the other works.

**This supplies independent expression measurements for a comparator, but reuses the original outcome labels. It is not independent clinical replication. The paired AUC difference has an interval spanning zero (October 1 audit).** The project's
central and least intuitive claim survives the transplant.

### 27.3 A discriminant check using a table this project did not compile

The standing weakness of every outcome test here is that the safety labels are literature
judgements made by one person (§23.8). The same paper supplies **Supplementary Data 1**: a
registry of **623 registered clinical trials against 74 ADC targets**, with phase — compiled by
someone else entirely.

It cannot serve as a safety label: "reached phase 3" measures development progress, not
toxicity, and HER2 has 36 phase-3 trials *and* a fatal on-target event. But it makes an
excellent **specificity control**, and the prediction was stated before the run: *confinement
should NOT predict which targets reached phase 3*, because if it did, confinement would really
be tracking "targets people invest in".

| question, on overlapping targets | separation | exact p |
|---|---|---|
| does confinement predict reaching **phase 3**? (24 targets: 11 reached, 13 did not) | **0.448** | 0.675 |
| does confinement predict **on-target toxicity**? (the 18 labelled targets) | **0.778** | 0.025 |

**As predicted.** Confinement is uninformative about whether a programme advanced and
informative about whether it poisoned normal tissue — checked against an externally compiled
table. This does not remove the need for an outside reviewer on the safety labels, but it does
show the measure is not simply tracking target popularity.

### 27.4 What this dataset does not settle

It is RNA-seq, with immunohistochemistry on only 8 of the 72 targets, so it does **not** speak
to section 26's protein question. It is breast cancer only. The normal samples are post-mortem,
so RNA degradation is a real concern the original authors address but this reanalysis inherits.
And 18 labelled targets (9 safe, 9 toxic) is a small panel: the estimate is 0.778 with a 95%
bootstrap interval of 0.52-0.96, so p = 0.025 is real but not far from 0.05, and it should be
read as consistent with §23 rather than as independent confirmation at full strength.

**Scripts:** `adc_autopsy_ceiling.py`, `trial_registry_control.py`. Results in the matching
`.json` files. Data: `ev_capture_data/adc_autopsy/`.

---

## 28. The two-antigen result in a solid tumour

Sections 17, 21 and 25 established the two-antigen findings in blood cancers. This section
repeats the whole analysis in metastatic breast cancer against normal tissue from the same
autopsies — 72 ADC targets, **2,556 antigen pairs**, 749 metastatic and 96 normal samples.

**Reconstruction check first.** The per-sample matrix is rebuilt from the paper's long-format
table and verified against their own wide table before anything is computed: **5,992 of 5,992
values matched exactly**. The check is a hard assertion in the code, not a printed warning.

**Scope note.** Bulk samples, not single cells, so these numbers are not comparable to the CyTOF
results in sections 17–21. Only 96 normal samples, so held-out halves are small and the
shuffled-label null is what makes any of it readable.

### 28.1 The theorem replicates

| check | result |
|---|---|
| SINGLE ≤ max(AND, NOT) ≤ LIKELIHOOD RATIO, in-sample | **0 violations in 5,112 checks** |

The hierarchy is a mathematical consequence of Neyman–Pearson, so violations would mean a coding
error rather than a biological discovery. Zero, in a second disease on a third technology.

### 28.2 What a second antigen buys, against the shuffled-label null

| quantity | real | shuffled labels |
|---|---|---|
| best single buildable antigen, held out | **J = 0.398** | −0.014 |
| statistical gain over the best antigen *in that pair*, median | +0.0000 | +0.0062 |
| engineering gain over the best kill-high antigen *in that pair*, median | +0.0000 | +0.0051 |

**The statistical result replicates exactly.** A second antigen adds no new information about
which samples are tumour — median gain at, and in fact below, the noise floor. This is now the
same answer in two diseases, two technologies and two laboratories.

**The engineering result does NOT replicate at the median, and that is a real difference.** In
DLBCL the median engineering gain was +0.106; here it is 0.000. The gain has not vanished — it
has become **concentrated**: 34.8% of the 2,556 pairs still show a gain above 0.05, and the 90th
percentile is +0.277. In the lymphoma cohort nearly every pair benefited, because most of the
panel is B-lineage and many markers sit *lower* on tumour than on normal B cells. Across organs
that is no longer true of most pairs, so the benefit is real for a minority and absent for the
rest.

### 28.3 The number that actually matters, and a trap avoided

The per-pair "engineering gain" figures look enormous at the top of the ranking — ALCAM + PD-L1
reaches +0.693 — but that headline is **misleading, and the reason is worth stating** because it
is the same trap section 20 recorded. The one-receptor baseline *inside those pairs* is
**negative** (−0.16, −0.135): both members are individually poor or inverted, so the comparison
is against something worse than chance. A large ratio over a terrible baseline is not a large
result.

The honest cohort-level comparison is against the best buildable single receptor **across all 72
antigens**, not within the pair:

| design | held-out J |
|---|---|
| best single activating receptor, any of 72 antigens | 0.398 |
| **best two-receptor design** | **0.533** |
| net gain | **+0.135** |

Which is strikingly close to DLBCL's +0.133 — arrived at from completely different data.

Every one of the top-ranked designs is a **NOT gate**: an activating receptor against one
antigen, paired with an *inhibitory* receptor against a second that is higher on normal tissue.
That is the same concrete design conclusion section 21 reached in lymphoma (activating CAR plus
an inhibitory receptor against IgD), reached independently in a solid tumour.

### 28.4 Limitations of section 28

Bulk samples; 96 normal samples split in half for training and testing; breast cancer only;
post-mortem RNA. The per-pair rankings in `adc_two_antigen_top_pairs.json` should be read as
hypothesis-generating, not as design recommendations — several involve antigens (PD-L1, EPHA2)
whose within-pair baseline is negative, and a real design would start from the cohort-level
comparison in 28.3.

**Scripts:** `adc_two_antigen.py`. Results in `adc_two_antigen.json` and
`adc_two_antigen_top_pairs.json`.

---

## 29. The weighted multi-compartment ceiling: the two measures were the same theorem all along

> **Correction (September 29, 2026): the heading and the "they are the same statement" claim below go too far.** The theorem (ceiling = ROC hull against a weighted
> normal mixture; 0 violations in 16,400 checks) stands. But the bridge in §29.3 item 3 does not hold as written. The linear form H ∝ Σ_c w_c s_c (weight × share of the gene's
> expression) does not reward confinement: if all weights are equal it equals 1 for every gene, and with binary weights it depends only on how much expression sits in
> expendable cell types, which is measure V3 in §37.2 (separation 0.699, far below top share's 0.880). So `top_share` is not the "leading term" of that formula. Confinement follows from the
> *threshold* form of the same harm: a therapy kills a cell when its antigen level passes a threshold, so H(t) = Σ_c w_c π_c q_c(t) (π_c = the cell type's share of the body, q_c = the fraction
> of its cells above t); the more cell types carry the antigen, the more terms, so breadth matters. Top share is a peak-based stand-in for breadth, and confinement is an empirical finding
> *motivated by* the theorem, not derived from it. §38.5 tests the threshold form directly. Wording to use in the paper: the theorem gives the ceiling and shows why expendability matters;
> the confinement rule is an empirical result that the threshold form of the theorem makes plausible.

The project has been carrying two measures that both predict clinical safety and had no formal
relationship to each other:

- **the ROC ceiling** (sections 1–22): for any kill rule monotone in antigen, the achievable
  (tumour killed, normal killed) points lie in the convex hull of the tumour-versus-normal ROC.
  This is a theorem.
- **confinement plus expendability** (section 23): a target is safe only if concentrated in few
  normal cell types *and* those cell types are ones the body can spare. This was an empirical
  rule with a mechanistic story attached.

A judge is entitled to ask why the project needs both, and whether the second is a heuristic
bolted onto the first. **They are different: a response bound under explicit assumptions and an empirical ranking hypothesis. The bound motivates investigation but does not derive the score.**

### 29.1–29.3 Corrected scope (October 2, 2026)

These mathematical statements supersede the earlier formulation; the original is archived. With fixed nonnegative weights and a shared response, H/W = E_pW[f] for W>0. For nondecreasing antigen responses, the feasible region is the convex hull of antigen-threshold ROC points. Arbitrary responses have the larger all-subsets feasible region, whose upper boundary is ordered by likelihood ratio. AUC alone determines neither full frontier. W=0 needs separate treatment. The harm functional is a specified model, not a calibrated measure of clinical injury.

One-hot weighting recovers each compartment's reference, but the minimum of individually optimized ceilings need not be attainable under simultaneous compartment constraints. The expression-share identity does not derive the empirical two-factor score. The corrected proof, assumptions, counterexamples, and verification are in `strengthening_Oct2/mathematical_scope.md` and `frontier_results.json` in the outputs root.

### 29.4 Numerical verification

Proofs do not execute, so the algebra is checked against brute force on random problems, exactly
as section 8 does for the single-compartment case.

| check | result |
|---|---|
| largest amount by which any randomised kill rule exceeded the weighted ROC ceiling | **0.000e+00** |
| largest gap between likelihood-ratio thresholding and brute-force optimum over all subsets | **0.000e+00** |
| total checks (1–6 compartments, random weights, random rules) | **16,400** |
| single-compartment case reproduces the original theorem | ok |
| worst-compartment case recovered | ok |
| confinement/expendability identity | ok |

**A bug in the verifier, recorded.** The first version of the hull construction used an ad-hoc
"upper envelope" loop that was wrong and reported violations of up to 0.96. The bug was in the
verification code, not the theorem. It is noted here because a verification script that can
produce false violations is worse than no verification at all; the hull is now a standard
monotone-chain upper convex hull.

### 29.5 Why this matters for the write-up

This is the section that converts the project from *"a measure that predicts drug safety"* into
*"a theorem about what any antigen-targeted therapy can achieve, with an empirical rule that the
theorem makes plausible."* It also makes the expendability weight a **modelling choice a reader can
disagree with explicitly** — set your own w, get your own ceiling — rather than a judgement
buried in a label file.

**Script:** `weighted_ceiling_theorem.py`; results in `weighted_ceiling_theorem.json`.

---

## 30. Discovery: what the framework says the field should be targeting

Every result so far concerns targets other people chose. The question a reader asks once the
framework works is the obvious one: **so what should we be targeting that nobody is?**

This section answers it across all 19,475 expressed human genes, using four filters applied in a
fixed order, with no fitted parameter — the two-factor rule and its 0.25 threshold are the frozen
ones from §23.2 and the locked predictions.

| filter | genes remaining |
|---|---|
| start | 19,475 |
| **1.** observed cell-surface localisation (HPA immunofluorescence plasma membrane, or CD-marker list) | 1,889 |
| **1b.** above the CAR-T sensitivity floor in cancer cell lines (≥ 20 nTPM) | — |
| **2.** confined **and** expendable: `top_share + 1.0 if haematopoietic home` ≥ 0.25 | 821 |
| **3.** elevated in cancer **cell lines** | 482 → 242 after the floor |
| **4.** not already pursued (623-trial ADC registry, the frozen label file, known CAR-T antigens) | **214** |

### 30.1 Two failures of my own, both instructive, both kept

**Bulk tumour RNA is contaminated, and the first screen rediscovered its own filter through the
contamination.** Using HPA's bulk-tumour cancer specificity, the top of the list was haemoglobin
(HBA1/HBA2/HBB) "elevated in glioblastoma", immunoglobulins (IGHG1–4, IGHA1) "elevated in lung
adenocarcinoma", and microglia and neutrophil markers (P2RY12, CX3CR1, TREM1/2, FPR1) "elevated
in glioblastoma". None of those are tumour proteins. A bulk tumour sample contains blood and
infiltrating immune cells, so **any gene confined to a blood lineage looks cancer-elevated
precisely because it is confined to a blood lineage** — filter 3 was re-detecting filter 2
through a contaminant. Fixed by requiring elevation in cancer **cell lines**, which are pure
tumour cells.

**Perfect confinement of nothing is still nothing.** A gene expressed at 1.6 nTPM in exactly one
cell type scores `top_share` = 1.000. With no absolute floor the list filled with olfactory
receptors. Aim 10's calibration puts the CAR-T sensitivity floor near 10²–10³ molecules per
cell, so an expression floor is not an arbitrary cut-off — it is the model's own requirement.
Also, "predicted membrane protein" turned out to be too weak a surface criterion: it admitted
haemoglobin, which sits inside the red cell, and secreted immunoglobulin constant regions, which
the plasma cell exports rather than displays. Observed localisation replaced prediction.

### 30.2 The positive control

A discovery screen with no positive control is a wish list. There is no ground truth for "targets
that should exist but don't" — but there is a checkable proxy: **immunotherapy targets already in
active clinical development that are absent from the ADC registry the screen excluded.** If the
screen works, it should surface them without having been told anything about them. The list was
written from the clinical literature before ranks were looked at.

| target | rank of 214 | what it is |
|---|---|---|
| **LILRB4** | **6** | AML/myeloma CAR-T and antibody (IO-202), active trials |
| **GYPA** | **9** | glycophorin A, erythroleukaemia |
| **TNFRSF13C** | **16** | BAFF-R CAR-T, B-cell malignancy trials |
| **CD27** | **20** | varlilumab |
| **TNFRSF13B** | **23** | TACI, myeloma bispecific/CAR target |
| **CSF3R** | **24** | G-CSF receptor, AML |
| **CD79A** | **29** | partner of CD79B, which polatuzumab targets |
| **CLEC12A** | **38** | CLL-1, AML CAR-T and bispecific trials |
| CD37 | 85 | naratuximab |
| SIGLEC6 | 181 | AML CAR-T |
| CD5, CD7, KIT, FLT3, HAVCR2, ADGRE2 | absent | — |

**10 of 16 recovered; 8 in the top 50.** The screen keeps 214 of 19,475 genes (1.10%), so a
random screen would recover **0.18** of these by chance. Observing 10 is a **57-fold enrichment,
hypergeometric P(X ≥ 10) = 1.6 × 10⁻¹⁶.**

**And the misses are arguably correct calls, not failures.** CD5 and CD7 CAR-T cause T-cell
aplasia; FLT3 and KIT CAR-T cause myeloablation. Those are exactly the on-target toxicities the
framework exists to flag, and the screen declined to recommend them.

> **Update, September 28, 2026 (candidate hygiene).** The 214 above is the raw screen output. Inspecting
> it showed it contained genes that are not cell-surface targets (immunoglobulin kappa variable genes,
> multi-copy GAGE genes, secreted proteins). A filter with rules fixed from HPA's own protein-class
> fields (`candidate_hygiene.py`) was added and *measured on known targets*: it falsely excludes 6 of 107
> clinically pursued surface targets (at most 6%; the first version excluded 27 of 107 and was
> abandoned). Funnel: 214 → 208 (widened exclusion set) → **134 pass hygiene → 128 also not already
> pursued.** All 10 positive-control recoveries survive. Below about rank 10 the ordering is sensitive to
> the expendability weights (only 9 of the top 25 stay in the top 25 in ≥ 90% of perturbed draws), so
> ranks are not precise. The paragraph below describes the raw list, including the immunoglobulin
> genes the hygiene filter now removes. This screen is a hypothesis generator; it knows nothing about
> internalisation, shedding, epitope access or tumour dependence.

### 30.3 What comes out, and what it is worth

Of the 214 candidates, **134 sit in an expendable blood lineage** — the pattern every approved
blood-cancer target follows. Beyond the rediscovered ones, the list is topped by immunoglobulin
kappa variable genes in plasma cells (IGKV3D-15, IGKV3-20, IGKV3-11) at myeloma cell-line levels
of 1,300–5,500 nTPM. That is not a coincidence with this project's own earlier work: §15's
light-chain clonality analysis measured κ:λ restriction directly, and light-chain-restricted CAR-T
is a real clinical strategy. The screen reached the same place from the opposite direction.

**What this list is not.** Expression is necessary and nowhere near sufficient. The screen knows
nothing about internalisation, epitope accessibility, shedding into serum, or whether the tumour
depends on the protein. It is hypothesis-generating with a fully checkable derivation — which,
with a 57-fold enrichment for known targets behind it, is the honest description of a lead.

**Scripts:** `target_discovery.py`; results in `target_discovery.json` and
`target_discovery_positive_control.json`.

---

## 31. The protein-at-cell-type question, and why it stays open

Section 26 left one question: the result holds on RNA at cell-type resolution (0.868) and on
antibody staining at tissue × cell-type resolution (0.785), and fails at whole-organ resolution
in **both** technologies (0.626 RNA, 0.650 protein). The open case was **mass-spectrometry
protein at cell-type resolution**, which did not exist when §26 was written.

It exists now. HPA v25.1 integrated **Deep Visual Proteomics**: protein intensity per cell type,
27 cell types across 14 tissues, ~19,000 genes, and its per-gene summary ships inside
`proteinatlas.tsv.zip`.

**It cannot answer the question, and the reason is a property of the data.**

| test | separation | p |
|---|---|---|
| PROTEIN cell-type specificity (ordinal ladder) | 0.583 | 0.219 |
| RNA cell-type specificity (same ladder, matched control) | 0.553 | 0.254 |
| RNA cell-type `top_share` from the same table | 0.716 | 0.031 |
| PROTEIN whole-tissue ladder | 0.636 | 0.141 |
| RNA whole-tissue ladder | 0.616 | 0.131 |

**Two things make the protein row uninformative rather than negative.**

1. **The safe class is below the detection limit.** 8 of 29 labelled targets are *Not detected*
   in the cell-type proteome, and **7 of the 8 are safe targets** — BCMA, SLAMF7, CD30, GPRC5D,
   CD52, DLL3, GPC3. Those are precisely the confined blood-lineage antigens the framework is
   built on. DVP is 27 cell types from a **single donor**; plasma cells and activated lymphocytes
   are rare, so the antigens that define the safe class are simply not seen. Dropping them does
   not fix the problem — it removes the class being tested.
2. **The ordinal ladder is too coarse to detect the effect where the effect is known to exist.**
   On the RNA cell-type data, where the full continuous computation gives 0.868, HPA's four-level
   ladder gives only 0.553. The same table's continuous `top_share` gives 0.716 (p = 0.031). So a
   null on the ladder carries almost no information either way.

**A scoring error of mine, recorded.** The first version mapped *Not detected* to the bottom of
the ladder — the opposite of what it means — which put BCMA, SLAMF7, CD30, GPRC5D and CD52 at
"least specific of all" and produced a meaningless 0.413. *Not detected* now returns a missing
value and those targets are excluded rather than mis-scored.

**Conclusion, and it is a statement about the data rather than the framework.** No protein
reference at cell-type resolution currently covers the cell types this question needs. The
antibody-staining result (0.785, p = 0.010, ~96 tissue × cell-type groups) remains the best
protein-level evidence, and §26.5's scope condition stands unchanged. Downloading the full DVP
tables would not change this — it is the same single-donor 27-cell-type dataset, with the same
detection limit.

**Script:** `hpa_summary_reference.py`; results in `hpa_summary_reference.json`.

---

## 32. The protein question, answered as far as any existing data allows

Section 26 left one question open and section 31 could not close it. This section closes it as
far as the world's data permits, and the answer is a limitation rather than a win.

### 32.1 What exists, after a proper search

The question needs **protein** measured across **normal human cell types**, covering both blood
lineages and solid tissue, for the targets in the panel. What exists:

| reference | what it is | why it cannot answer the question |
|---|---|---|
| HPA Deep Visual Proteomics (May 2026) | mass spec, 27 cell types, 14 tissues | single donor; 7 of the 8 undetected targets are the safe blood-lineage antigens (§31) |
| Immune-cell proteomes (e.g. Rieckmann 2017, 28 immune cell types with copy numbers) | mass spec, blood only | no solid tissue, so the toxic class is absent |
| HPA immunohistochemistry | antibody, ~96 tissue × cell-type groups | **usable, and already run: 0.785, p = 0.010 (§26.1)** — ordinal, four levels |
| HPA mass spectrometry | protein, 20 whole organs | wrong resolution; fails, as RNA also does there (§26.4) |
| **This project's own CITE-seq panel** | **antibody, 228 surface proteins, 161,764 cells, 8 donors, 31 blood cell types** | **blood only — but it is a genuine protein cell-type reference and had not been used for this** |

**No single dataset measures protein across both blood and solid-tissue cell types with the
coverage this question needs.** That is a statement about the field in September 2026, not about
this analysis.

### 32.2 Using the reference the project already owned

The CITE-seq panel from §16 is protein, measured by **antibody** — the same physical modality a
CAR or therapeutic antibody uses to bind — at single-cell resolution across 31 blood cell types.
It cannot run the safe-vs-toxic test (the toxic class is almost all solid-tissue targets), but it
can answer the question underneath: **does confinement measured on protein track confinement
measured on RNA, over the same cell types?**

185 genes measured both ways, same lineages, same measure:

| comparison | Spearman ρ | p |
|---|---|---|
| **PROTEIN (antibody) vs RNA, same blood cell types** | **0.167** | 0.024 |
| number of cell types, same comparison | 0.171 | 0.020 |
| *for scale:* RNA cell-type vs RNA bulk-tissue | 0.512 | §23.7 |
| *for scale:* RNA cell-type vs protein antibody staining | 0.429 | §26.3 |
| *for scale:* RNA cell-type vs protein mass spec | 0.346 | §26.3 |

**ρ = 0.167 is the weakest agreement of any pair this project has measured**, and it is measured
under the most favourable conditions available — same cell types, same measure, matched design.

Separating the two claims a confinement value makes:

- **The home** (which cell type carries the antigen): agrees for **9 of 13** therapy targets —
  CD19, CD20, CD22, CD79B, CD52 → B cells; CD38, SLAMF7 → plasmablasts; NCAM1 → NK; CD44 →
  monocytes. Overall agreement across all 185 genes is only **37%**, but much of that is a
  **roster** problem rather than a disagreement: the PBMC panel contains no tissue macrophages,
  neutrophils or Kupffer cells, so a neutrophil protein has nowhere correct to land.
- **The share** (how concentrated it is there): does not transfer, ρ = 0.167.

### 32.3 The honest conclusion, and what it costs the project

**The confinement measure is not technology-independent.** Antibody-measured and
transcript-measured confinement give substantially different numbers for the same protein in the
same cell types. The project cannot claim that its RNA-derived confinement values *are* protein
confinement values.

**What survives, stated precisely:**

1. The headline result (0.868, p = 0.0003) is a statement about **RNA-derived** confinement, and
   should be written that way throughout.
2. The result **does** replicate on an independent protein reference — antibody staining across
   ~96 tissue × cell-type groups, 0.785, p = 0.010. That replication is not weakened by anything
   here; a measure can predict an outcome without being numerically identical to its protein
   counterpart.
3. What cannot be claimed is the stronger version: *"confinement is a property of the protein,
   and RNA is just a convenient proxy."* On present evidence RNA-confinement and
   protein-confinement are correlated but genuinely different quantities.

**How much this costs.** Less than it might appear, and the reason is the outcome test. The
framework is validated by whether it predicts **clinical toxicity**, not by whether two
measurement technologies agree. It predicts toxicity from RNA (0.868) and from antibody staining
(0.785), which is two independent technologies reaching the same clinical conclusion. The
disagreement in 32.2 says the two are not measuring an identical physical quantity — which is
interesting, and belongs in the limitations, but does not undo a result that was never
conditioned on them being identical.

**What would settle it** is a protein atlas covering rare blood lineages and solid tissue at
cell-type resolution in the same donors. That does not exist yet. Deep Visual Proteomics is the
first attempt and is one donor deep; a multi-donor successor would answer this directly, and is
the obvious thing to revisit if one appears.

**Scripts:** `protein_celltype_citeseq.py`, `adt_gene_map.py`; results in
`protein_celltype_citeseq.json` and `protein_celltype_home_agreement.json`.


# Fifth addendum (September 9, 2026): an independent healthy-marrow reference for the myeloma side

Section 16 closed a gap for CD19/DLBCL: an external, independent, healthy-donor dataset (161,764
PBMCs, 8 donors, not this project's own patients) confirming the four approved lymphoma targets
are lineage-restricted in genuinely healthy people. No equivalent existed for the myeloma/BCMA
side -- every myeloma comparison so far has been **within-patient** (tumour plasma cells versus
normal compartments from the *same* bone marrow biopsy, sections 15 and 19). That is real data,
but it cannot rule out a batch- or patient-specific artifact the way an independently collected
healthy cohort can. This section closes that asymmetry, as far as the available public data allows.

## 33. An independent healthy bone marrow cohort (Boiarsky et al. 2022)

**Data.** Boiarsky et al., *Nature Communications* 2022 (doi:10.1038/s41467-022-33944-z, PMID
36396631), GEO accession **GSE193531** -- CD138-sorted (plasma-cell-enriched) scRNA-seq from the
Ghobrial lab at Dana-Farber, spanning normal bone marrow (NBM), MGUS, smouldering myeloma (SMM),
and myeloma (MM) donors. **29,387 total cells, 9 NBM (genuinely healthy, non-myeloma) donors
contributing 9,329 of them.** Freely downloadable from GEO, no controlled-access application
required -- the strongest candidate found across a three-round search of four public datasets (the
others: Zhang et al. GSE245108, too few plasma cells -- 135 of 72,179 total cells; Duan et al., no
accession located with a public protein layer; Da Vià et al. 2025, larger and protein-annotated
but EGA controlled-access, accession EGAC50000000232, not an instant download).

**Honest scope limit, stated before anything else.** This is scRNA-seq (transcript counts), not
CITE-seq or flow cytometry -- RNA-level expression, not surface protein or absolute molecule
counts. It is one more data point in the same category as section 15/19's RNA-level within-patient
comparisons, not a protein-level replication. It is independent in the sense that matters here --
different donors, different lab, different sequencing run, never in contact with this project's
own cohort -- but it is not a different measurement modality.

**A second, unplanned scope limit, found on inspection.** This matrix's cell composition is not
whole bone marrow. Checking canonical lineage markers across all 29,387 cells: plasma-cell markers
CD38 and SDC1/CD138 are detected (>0 counts) in 67.3% and 60.5% of cells respectively, while
T-cell (CD3D 0.45%, CD3E 0.33%), myeloid (CD14 1.04%), B-cell (MS4A1 3.96%, CD19 5.43%) and NK
(NCAM1 11.63%) markers are all in the low single digits or near-zero. **This confirms the file is
already CD138-sorted/plasma-cell-enriched**, consistent with the study's own design, not whole
marrow requiring a T/myeloid/NK/B gating step the way this project's own myeloma CITE-seq data
does (`myeloma_ceiling.py`'s gating). That rules out repeating section 16's exact
restriction-AUC test here -- there is no meaningful population of other lineages left in the file
to compute restriction against.

**What this dataset can still test, and what it found.** Restricting to the 9,329 NBM
(healthy-donor) cells: **BCMA (TNFRSF17) is detected in 82.87%** of them (median 2 counts, 90th
percentile 7, max 37) -- more frequently detected than the plasma-cell markers used to confirm the
population's identity (CD38 73.78%, SDC1/CD138 57.44%). The 15.0% of NBM cells negative for both
CD38 and SDC1 are not a distinct, non-plasma-cell population: their median sequencing depth is
roughly half that of the CD38+/SDC1+ majority (5,627 vs. 10,288 median UMI counts; 907 vs. 1,302
median genes detected), the signature of dropout at lower depth, not a different cell type -- and
BCMA is still detected in 72.6% of them, along with near-total absence of T/myeloid/NK/B markers
(all under 7%, most under 1%).

**Read plainly: in an independently collected healthy human bone marrow cohort never seen by this
project before, BCMA is a robust, near-ubiquitous marker of the plasma-cell compartment it was
sorted for -- not a project-internal artifact of this project's own 11-sample cohort.** That is a
real, if modest, form of external validation: the antigen this project has used throughout the
myeloma analysis behaves the same way in someone else's data.

**For contrast**, the same NBM cells were checked for the two other myeloma antigens with
extracted markers: **SLAMF7 (50.28% positive, median 1, p90 2)** -- moderate, consistent with
SLAMF7's known broader expression (also on NK cells and some activated lymphocytes, which this
sorted dataset cannot test for) -- and **GPRC5D (5.67% positive, median 0, p90 0)** -- sharply lower,
consistent with GPRC5D's known status as a more plasma-cell-specific but lower-abundance target in
the literature. The pipeline distinguishes near-ubiquitous (BCMA) from genuinely sparse (GPRC5D)
expression within the same sorted population, which is some reassurance the 82.87% BCMA figure is
a real signal and not a ceiling effect of the extraction method.

**What this does *not* show, stated as plainly as the positive finding.** This is not a
restriction test -- it cannot show BCMA is *absent* from non-plasma-cell lineages, because the
sorted design leaves too few of those cells in the file to test against (the CD38-/SDC1- subset is
almost certainly under-sequenced plasma cells, not an independent other-lineage sample). It closes
the "is this an artifact of this project's own cohort" question, not the "is BCMA safe to target"
question -- section 11a's calibration-gap finding and the within-patient tumour-vs-normal
comparisons (sections 15, 19) remain the load-bearing analyses for target selectivity; this section
adds one more independent confirmation that the underlying antigen biology is real, not an
artifact.

**Scripts:** `extract_boiarsky_markers.py` (one-time extraction of 11 target-gene rows and cell
metadata from the raw GEO matrix into `extracted_markers.json`), `boiarsky_healthy_reference.py`
(the analysis above). Data: `ev_capture_data/healthy_bone_marrow_boiarsky/` (GSE193531 UMI count
matrix and cell-level metadata, downloaded directly from GEO).

---

# Sixth addendum (September 10, 2026): the panel grown to 44, and an honest look at what that did to one control

Two more targets (ERBB3, AXL) were added to the frozen label file after a systematic re-check of
the remaining unlabelled genes in the De Schepper et al. ADC trial registry (see section 23's
panel-growth note for the full account of both). The core safe-vs-toxic result (sections
23.1–23.3, 23.5–23.7) is completely unchanged by this, because both new targets are
`efficacy_fail`, outside every one of those tests.

## 34. What did change: the specificity control, and why it is now a disclosed limitation

Section 23.4 is the one place in the document these two additions land. Adding ERBB3
(`top_share` 0.083) and AXL (0.102) — both low-confinement, broadly-expressed solid-tumour
antigens — pushed the safe-vs-efficacy-fail specificity control from borderline (0.700,
p = 0.043 on the 42-target panel) to clearly significant (0.746, p = 0.011). Investigated
directly rather than left as a number: the efficacy-fail group has always split into two
biological clusters — a small immune-lineage subgroup (FCER2/CD23, HLA-DRA) that scores like the
safe group, and a larger broadly-expressed-solid-antigen subgroup that scores like the toxic
group — and that split predates this session's additions (the 10-target efficacy-fail group's
median `top_share` was already 0.183 before ERBB3 and AXL were added). Both new targets happen
to land in the low subgroup, which is why the control moved.

**The honest reading, stated once here and cross-referenced from 23.4 and 23.8:** a
broadly-expressed, low-confinement antigen may be difficult to build a wide therapeutic window
around at all, whichever way a given trial happens to fail — as on-target toxicity, or as an
efficacy failure that masks a dose that was capped for tolerability reasons before efficacy
could be shown. That is a plausible extension of the project's own thesis, not a contradiction
of it, but it does mean this one specificity control no longer cleanly demonstrates "confinement
is uninformative about efficacy failure in general." It does not touch the headline safe-vs-toxic
result, which never uses the efficacy-fail group at all.

## 35. A wrinkle recorded on an existing label, not acted on

While researching this round, a second, more recent drug against an already-labelled target
(SLC34A2/NaPi2b) was found to carry a real, serious, and never-fully-resolved safety signal (an
FDA partial clinical hold following fatal bleeding events). The label was not changed — the
drug's own program ended for efficacy, and the bleeding's mechanism was never confirmed on-target
— but the ambiguity is now recorded in `target_labels.py`'s evidence field for that gene, and
summarised in section 23's panel-growth note, so a future reviewer sees the full picture rather
than a clean-looking label that doesn't reflect everything found.

**Scripts:** no new scripts; `wholebody_v2.py` re-run unchanged against the updated
`target_labels.py`. `test_wholebody_restriction.py`'s panel-size tripwire test and
`verify_documented_numbers.py`'s panel-composition check were both updated to 19/11/12 and to
also check the section 23.4 separation value, so future drift is caught the same way panel
composition already was.

---

# Seventh addendum (September 29, 2026): continuous expendability weights, tested twice, and why the yes/no flag stays

Section 29 proves the ceiling theorem for a continuous expendability weight on every normal cell type, but the
empirical rule (§23.2) uses a yes/no flag. This addendum records the two attempts to replace the flag with something
continuous. Both schemes were fixed before they were run (`continuous_weights_v2_precommitment.md` for the second).

## 36. Attempt 1: weight only the home cell type

Per-cell-type weights (0 = expendable, 1 = catastrophic), assigned from regeneration, clinical replaceability and
survivability with a rationale recorded per entry (`expendability_weights.py`), applied to the target's home cell type:
separation **0.809** against the binary flag's 0.914 (a regeneration-only variant: 0.842). The result did not beat the
same weights randomly reassigned across cell types (p = 0.15). Cause: for a gene spread over dozens of cell types the
computed home is close to arbitrary (EGFR's is placental tissue, scored fully expendable), so a toxic target collects a
large bonus. The yes/no flag is coarse, which makes it robust to that.

## 36.1 Attempt 2: use the whole expression distribution, as the theorem does

The theorem's harm is the sum over cell types of weight × share of the gene's expression in that cell type, so it does
not need a home at all. Four schemes were fixed in advance (E1 hand weights only; E2 fraction of expression in blood
lineages; E3 confinement plus expression-weighted hand weights; E4 confinement plus the blood-lineage fraction):

| scheme | separation | 95% bootstrap CI | paired difference vs binary | p vs shuffled weights |
|---|---|---|---|---|
| binary flag (§23.2) | 0.914 | 0.789–0.995 | — | — |
| confinement only (§23.1) | 0.880 | 0.732–0.986 | — | — |
| E1 | 0.818 | 0.651–0.947 | −0.096 (−0.230 to +0.010) | 0.004 |
| E2 | 0.699 | 0.498–0.871 | −0.215 (−0.440 to −0.010) | 0.007 |
| E3 | 0.885 | 0.742–0.986 | −0.029 (−0.096 to +0.019) | 0.30 |
| E4 | 0.880 | 0.742–0.976 | −0.033 (−0.139 to +0.053) | 0.29 |

**Decision, by the rule fixed in advance:** a continuous scheme would replace the binary flag only if it scored at least
0.914 and beat its shuffled-weight null (p < 0.05). None did, so **the binary flag stays primary.**
Predictions: P1 (E3, E1 beat 0.809) confirmed, E1 only marginally; P2 (E4 within 0.03 of binary) narrowly missed (−0.033);
P3 (E1 weaker than E3 by ≥ 0.02) confirmed; P4 (no hand-weight scheme beats binary by > 0.03) confirmed.

**What the pattern says.** E3 and E4 land at the confinement-only level (0.880) and do not beat their shuffled-weight
nulls: a soft expendability term adds nothing beyond confinement. The informative fact is categorical: 8 of 19 safe targets
have a blood-lineage home and 0 of 11 toxic ones do. A share-based measure loses it, because broadly expressed toxic
targets still have some expression in blood cells (CD44 0.31 of its expression is in blood lineages, CD47 0.16), which a
soft scheme scores as partly expendable, while confined solid safe targets have almost none (0.003–0.04). "Has its home in
blood" and "has some expression in blood" are different statements, and only the first separates the groups.

**Scope.** These four schemes plus the first attempt were the first seven looks at a 30-target panel; none was selected
for scoring well. The hand weights were assigned by me from stated biology. Weights from an outside dataset were tried
next (§36.2). The theorem holds for any weights, so the generality remains a statement about scope.

**Scripts:** `continuous_weight_sensitivity.py`, `continuous_weights_v2.py`; results in `continuous_weight_sensitivity.json`
and `continuous_weights_v2.json`; pinned by `tests/test_continuous_weights_v2.py`.

## 36.2 Attempt 3: weights from an outside dataset (measured cell lifespans)

Rounds 1 and 2 used weights I assigned. To remove that judgement, the third round used a published measurement made by
other people for another purpose: the lifespan of each human cell type (Sender & Milo, *Nat Med* 2021; data file kept in
`ev_capture_data/cell_turnover_sender_milo/`). Idea: a fast-turnover cell type is replaced by the body, so losing it is
expendable. Weight = log10(lifespan in days) / log10(300,000), clipped to 0–1 (neurons = 1, colon lining, 3.4 days, = 0.10).
The dataset has about two dozen cell categories, so I mapped the 154 reference cell types onto them by lineage and tissue
family only, and fixed the mapping before scoring (`external_turnover_weights.py`, precommitment and input hashes in
`external_turnover_precommitment.md`): 77 of 154 types map; the rest get the median weight of the mapped ones. The mapping
was not blinded in the strict sense (I knew the labels and earlier results when I wrote it), which is why it was fixed in
writing first. Six schemes, fixed in advance, same adoption rule as before:

| scheme | score | separation | 95% bootstrap CI | paired difference vs binary | p vs shuffled weights |
|---|---|---|---|---|---|
| binary flag (§23.2) | confinement + blood-home flag | 0.914 | 0.789–0.995 | — | — |
| confinement only (§23.1) | top share | 0.880 | 0.732–0.986 | — | — |
| X1 | 1 − Σ share × turnover weight | 0.541 | 0.301–0.770 | −0.373 (−0.579 to −0.177) | 0.57 |
| X2 | 1 − turnover weight of home cell type | 0.498 | 0.261–0.727 | −0.416 (−0.665 to −0.177) | 0.51 |
| X3 | confinement + X1 | 0.833 | 0.656–0.971 | −0.081 (−0.187 to −0.010) | 0.77 |
| X4 | confinement + X2 | 0.809 | 0.617–0.962 | −0.105 (−0.239 to +0.000) | 0.32 |
| X5 | confinement + [home lifespan ≤ 30 d] | 0.766 | 0.517–0.962 | −0.148 (−0.340 to −0.005) | 0.68 |
| X6 | confinement + [home lifespan ≤ 365 d] | 0.718 | 0.512–0.885 | −0.196 (−0.397 to −0.043) | 0.49 |

**Decision, by the rule fixed in advance: none of the six comes close, so the binary flag stays primary.** Turnover alone
(X1, X2) does not separate the groups at all (p = 0.37 and 0.52). Adding it to confinement makes the ranking *worse* than
confinement alone (0.72–0.83 against 0.880), and none beats randomly shuffled weights.
Predictions: P1 (none reaches 0.914) confirmed; P2 (X3/X4 not above 0.880 by more than 0.02) confirmed; P3 (X5 the weakest
binary flag) **wrong**, X6 was weaker (0.718 against 0.766); P4 (none passes the shuffled null) confirmed.

**Why it fails.** Turnover measures how fast a cell type is *replaced*, not whether losing it during the replacement gap is
survivable. The safe blood-lineage targets sit on B cells and plasma cells (63-day lifespan, weight 0.33), which is
*less* expendable by this measure than the colon lining (0.10) or the oesophageal surface (0.27) where the toxic targets
CEACAM5, EPCAM and ERBB2 have their home. So the measure ranks three toxic epithelial targets as more expendable than
the safe B-cell targets, which is the opposite of the clinical facts (B-cell loss is managed with antibody replacement;
loss of gut lining is not). Coverage is also limited: 11 of the 30 targets have a home cell type (for example EGFR, DLL3,
GPC3, MET) that has no lifespan in the dataset and default to the neutral weight.

**Scope.** Twelve alternative schemes had been tried at this point (two in round 1, four in round 2, six in round 3) against the
binary flag on a 30-target panel (further rounds are in §37). None was adopted. The binary flag was chosen before any of them, on grounds fixed in §23.2,
and the methods should record that these alternatives were tried. The finding that separates the groups is categorical (blood home 8/19 safe,
0/11 toxic); what would test whether it generalises is new labelled targets, not more weighting schemes.

**Scripts:** `external_turnover_weights.py`, `external_turnover_weights_test.py`; results in `external_turnover_weights_test.json`;
pinned by `tests/test_external_turnover_weights.py`.

# Eighth addendum (September 29, 2026): cleaning up the blood flag, and a search for a more effective rule

## 37. Round 4A: define the blood flag from the whole distribution, not one arbitrary "home" cell type

The §23.2 flag adds +1 when the single top-share cell type is on the frozen blood list. For a gene spread over dozens of cell
types the top one is close to arbitrary, so a broad gene can collect the bonus by luck. Two whole-distribution definitions
were fixed in advance (`flag_cleanup_precommitment.md`, `lineage_map.py`, `flag_cleanup.py`), with a cleanup decision rule:
adopt if separation is at least 0.900 and the paired interval against the home-cell rule includes 0.

| definition of "expendable, blood" | targets flagged (safe / toxic) | separation | paired difference vs home rule | p vs shuffled flags |
|---|---|---|---|---|
| home cell type is on the blood list (§23.2) | 8 / 0 | 0.914 | — | — |
| M1: blood lineages hold at least half of the gene's expression | 6 / 0 | 0.904 | −0.010 (−0.043 to +0.000) | 0.001 |
| M2: blood is the largest of ten lineage families (`lineage_map.py`) | 8 / 0 | **0.914** | +0.000 (identical flags on all 30 targets) | 0.008 |

M1 loses CD19 and CD22 (blood holds 47% and 40% of their expression, below the half-way mark). **M2 reproduces the home-cell rule
exactly on the panel while removing its dependence on one arbitrary cell type, so every headline number is unchanged.**
Predictions: P1 (M1 between 0.85 and 0.914) confirmed; P2 (M2 at least M1) confirmed; P3 confirmed. Both definitions meet the cleanup rule; M2 is
preferred because it changes nothing on the panel (a choice made after seeing that both passed, so it is recorded here as such).
A bug in the shuffled null for M2 was found and fixed before the reported result (originally blood-named types stayed in the
blood family after their flag was shuffled away); the observed separations were not affected.

## 37.1 Round 4B: does any other confinement or abundance measure do better?

Reasoning first (`rule_improvement_precommitment.md`): fractions of expression hide absolute abundance (DLL3 is nearly absent
from adult tissue, yet its share looks diffuse); top share depends on how finely cell types are split; one maximum summarises a
154-value profile. Eighteen schemes (nine measures, each alone and with the M2 flag) were fixed in advance and charged for
selection with a family-wise permutation test (shuffle labels, take the maximum over all schemes).

| measure | alone | with flag | | measure | alone | with flag |
|---|---|---|---|---|---|---|
| C0 top share (current) | 0.880 | 0.914 | | C5 −(number of types above 10% of peak) | 0.876 | 0.909 |
| C1 top lineage share | 0.545 | 0.742 | | C6 −log10(1 + total expression) | 0.761 | 0.799 |
| C2 1 − entropy | 0.794 | 0.847 | | C7 −log10(1 + highest expression) | 0.622 | 0.722 |
| C3 tau (Yanai 2005) | 0.880 | 0.914 | | C8 rank-mean of C0 and C6 | 0.856 | 0.890 |
| C4 top-3 share | 0.818 | 0.866 | | | | |

**No scheme beats the current rule; none was promoted.** Two findings. (i) Tau, the standard tissue-specificity index, is a
monotone function of top share for a fixed number of cell types, so it ranks identically (0.880 and 0.914); it was therefore not an independent look. Top
share *is* the field's specificity index in disguise, which is a useful thing to say in the paper. (ii) Measures that use more
of the profile (entropy, top-3, lineage-level) are worse, and adding absolute abundance lowers the ranking.
Predictions: Q1 (abundance alone at least 0.70) held for total expression (0.761) but not for highest expression (0.622); Q2 (entropy, tau, top-3 within 0.05 of top share) was wrong for entropy (−0.086) and top-3 (−0.062); Q3 confirmed; Q4 (some scheme beats 0.914) not met.

## 37.2 Round 4C: theory-consistent "vital-mass" measures and absolute-threshold breadth

The theorem counts harm only in vital cell types (weight 0 for expendable ones), and a therapy responds to antigen above a
threshold. Seven measures, each alone and with the flag, were fixed in advance (`rule_improvement_c_precommitment.md`).

| measure | alone | with flag |
|---|---|---|
| V1 −log10(1 + expression summed over vital types) | 0.775 | 0.794 |
| V2 −log10(1 + highest vital-type expression) | 0.703 | 0.718 |
| V3 share of expression in expendable lineages | 0.699 | 0.699 |
| V4 rank-mean of top share and V1 | 0.861 | 0.885 |
| D1 −(types with nCPM ≥ 1) | 0.775 | 0.823 |
| D2 −(types with nCPM ≥ 10) | 0.854 | 0.880 |
| D3 −(types with nCPM ≥ 50) | 0.823 | 0.861 |

None was promoted. Predictions: R1 confirmed; R2 confirmed; R3 (D1–D3 within 0.05 of 0.876) missed for D1 and D3; R4 confirmed.

## 37.3 What the search shows, and what it does not

**Forty-four alternatives have now been tried against the flag on this panel (42 distinct: tau is duplicated); the current rule (top share plus the blood flag) beat all
of them.** That is reassuring about the rule, and it is also the point at which further searching on these 30 targets stops being informative: with a null
maximum of about 0.74 across the family, differences of a few hundredths cannot be told apart.
Where the remaining errors are: the current rule misorders 16 of 209 safe/toxic pairs, 8 of them because of CEACAM5 (confined to
colon lining, toxic) and 4 because of CAIX; the safe targets it ranks too low are mostly broadly expressed epithelial ones (MET, Nectin-4, mesothelin, glypican-3).
Two things in the labels, not in the RNA, plausibly explain this: (1) *modality*: several toxic verdicts come from cell therapies and
T-cell engagers (HER2, CEA, CAIX and FAP CAR-T), which respond to very low antigen levels, whereas several safe verdicts come from antibody-drug
conjugates; (2) *tissue criticality*: loss of gut lining, bile ducts or cornea is dangerous in a way the blood flag cannot express. Neither is testable without more labelled
targets that record the modality. Any measure that reaches 0.944 on these 30 targets would need confirmation on new targets before it could be used.

**Scripts:** `lineage_map.py`, `flag_cleanup.py`, `rule_improvement.py`, `rule_improvement_c.py`; results in `flag_cleanup.json`, `rule_improvement.json`,
`rule_improvement_c.json`; pinned by `tests/test_rule_search.py`.

# Ninth addendum (September 29–30, 2026): new labelled targets, and whether "more realistic" rules could be losing to noise

## 38. New labelled targets (validation set)

Two independent literature searches, given no expression data, looked for labelled targets outside the 44-gene panel, applying the same frozen rules. They found few that
meet the rules cleanly, because most clinical on-target attributions are hedged by the authors (payload, lymphodepletion, cytokine release, or off-target cross-reactivity).
Labels, tiers, confidence and sources were frozen before any expression score was computed (`target_labels_validation.py`, `validation_precommitment.md`, hashes recorded).
**Primary set (pre-specified): safe MAGEA4, SSTR2, CCR4, CTAG1B, IL2RA; toxic GUCY2C, PSCA.** Hedged toxic (reported separately): TFRC, SLC44A4, MLANA, MAGEA12, EPHA2.
Excluded: PMEL (an approval and an on-target toxicity record on the same gene), SLC6A2, CD276, PRAME.

| result | value |
|---|---|
| frozen rule on the primary set (4 safe v 2 toxic scorable) | **0.875**, exact p = 0.133 (smallest attainable p with 8 pairs is 0.067) |
| top share alone, same set | 0.750 |
| both toxic targets score below the median safe target | yes |
| hedged toxic (4 scorable) v primary safe | 0.812 |
| post-hoc: CTAG1B and MAGEA12 have no measured expression in any reference cell type, so a share score is undefined; scoring an absent gene as safest gives | 0.900, p = 0.095 (5 v 2) |
| pooled with the 30 development targets (not independent) | 0.893 (23 safe v 13 toxic) |

Predictions V1 (at least 0.80), V2 and V3 (hedged toxic v safe at least 0.60) were all confirmed. **This set is far too small to be evidence on its own** (7 to 12 targets); it shows the rule was not
obviously overfitted, and it shows how scarce clean labels are. Two things follow for the paper: state the label rule and its hedging problem plainly, and treat the OSF-registered
prospective predictions as the real out-of-sample test.

## 38.1 Label-consistency check and flip sensitivity

The boxed-warning cross-check turned up an inconsistency an outside reader could raise: NECTIN4 (labelled safe) carries a boxed warning for fatal SJS/TEN; the ocular boxed warning appears on
F3 and FOLR1 (labelled toxic) and also on a BCMA antibody-drug conjugate (TNFRSF17, labelled safe from CAR-T and bispecific evidence). Labels are agent-specific verdicts, not gene-level facts.
Labels in the repository were not changed. Flipping each disputed label and re-scoring the frozen rule:

| labels flipped | safe / toxic | two-factor separation | one-factor |
|---|---|---|---|
| none (as frozen) | 19 / 11 | 0.914 | 0.880 |
| NECTIN4 to toxic | 18 / 12 | 0.926 | 0.889 |
| F3 and FOLR1 to safe | 21 / 9 | 0.884 | 0.852 |
| both of the above | 20 / 10 | 0.890 | 0.855 |
| TNFRSF17 to toxic | 18 / 12 | 0.838 | 0.810 |

The result survives every plausible relabelling (all p < 0.001); the largest drop is 0.076.

## 38.2 Can a truer rule lose to a coarser one because of noise? (simulation)

Yes, and the panel cannot tell. Simulated truth: two rules whose scores are correlated (0.9, as observed between the frozen rule and the whole-distribution alternatives), the truer one with a
population separation higher by 0.02 or 0.04 than a coarse one at 0.90. Share of simulated panels on which the truer rule scores LOWER, and share on which a paired test (DeLong, p < 0.05) detects it as better:

| panel (safe / toxic) | truer by 0.02: scores lower / detected | truer by 0.04: scores lower / detected |
|---|---|---|
| 19 / 11 (this project) | 22% / 3% | 8% / 6% |
| 95 / 55 | 5% / 32% | 0% / 89% |
| 190 / 110 | 1% / 60% | 0% / 100% |

At the size of this panel a genuinely better rule (by 0.02 to 0.04) loses one time in five to one in twelve and is almost never *detected*. So a lower score on 30 targets does not show that a rule is
worse, and it equally cannot show that a more realistic rule is better. Detecting a 0.04 improvement needs roughly 150 to 300 labelled targets, which do not exist. Script: `rule_comparison_noise_simulation.py`.

## 38.3 Judging realism with a different, richer test: does the model put the harm where it fell?

Because the safe/toxic bit carries little information per target, realism was also judged on tissue location (`independent_tests_precommitment.md`, `damage_tissue_map.py`, frozen first). For each of the
13 toxic targets (11 panel, plus GUCY2C and PSCA), the model assigns every cell type a harm contribution; the score is the best rank of any documented-damaged cell type (1 = perfect).

| model | harm contribution | median best rank | targets with damaged tissue in top 3 | vs raw shares (better / worse) | sign-test p |
|---|---|---|---|---|---|
| N1 | share of expression | 7 | 4 of 13 | — | — |
| N2 | share, expendable lineages set to zero | 7 | 5 of 13 | 6 / 1 | 0.125 |
| N3 | share × hand-assigned weight | 4 | 5 of 13 | 8 / 2 | 0.109 |
| N4 | share × measured-lifespan weight | 5 | 4 of 13 | 5 / 7 | 0.774 |

**The hand-weighted model locates the damaged tissue better than raw shares (median rank 4 against 7, better on 8 targets and worse on 2), which is the direction your hypothesis predicts, but with 13 targets it is
suggestive, not significant (p = 0.11).** Measured lifespan does not help. CD47 is the informative failure: erythrocytes are "expendable" by lineage yet anaemia was dose-limiting, so zeroing them pushes the damaged tissue
from rank 34 to rank 137. Prediction check: N2 improved CA9 and EGFR and worsened CD47 (confirmed); N4 was the worst (confirmed by the sign test); N3 within one target of N1 (missed, it improved eight).
Candidate rules built on the vital-share idea, R1 = 1 − (largest vital-cell-type share) and R2 = R1 + blood flag, failed on outcomes (separation 0.435 and 0.507), because concentration in one cell type is not harm; **neither was adopted.**

## 38.4 A rule derived from the threshold form of the theorem

Because of the correction at §29, harm was also written in its threshold form and tested directly (`threshold_harm_precommitment.md`): TH1(t) = −(number of vital cell types with nCPM ≥ t); TH2(t) = the same with hand weights. Thresholds t = 1, 10, 50 nCPM.

| scheme | separation | paired difference vs frozen rule | validation (4 v 2) |
|---|---|---|---|
| frozen rule (top share + blood flag) | 0.914 | — | 0.875 |
| TH1(1) | 0.804 | −0.110 (−0.254 to −0.014) | 1.000 |
| TH2(1) | 0.785 | −0.129 (−0.278 to −0.019) | 1.000 |
| TH1(10) | 0.856 | −0.057 (−0.160 to +0.014) | 1.000 |
| TH2(10) | 0.837 | −0.077 (−0.191 to +0.005) | 1.000 |
| TH1(50) | 0.859 | −0.055 (−0.167 to +0.029) | 1.000 |
| TH2(50) | 0.833 | −0.081 (−0.201 to +0.010) | 1.000 |

None was promoted (family-wise null maximum about 0.77). The theory-faithful threshold rules score 0.05 to 0.08 below the frozen rule; for t = 10 and 50 the paired interval includes 0, so they are statistically tied with it, and
they score a perfect 1.000 on the 8 validation pairs, which is not informative at that size. Predictions S1 (all between 0.78 and 0.91), S2 and S3 were confirmed.

## 38.5 What this means

(1) The frozen two-factor rule remains the pre-specified primary rule: it was fixed before any alternative, it is simple, and nothing beats it; 52 alternatives had been tried at this point (50 distinct); section 39 brings the count to 74.
(2) The data cannot separate it from more realistic rules; the tissue-location test hints that hand-weighted versions represent the biology better, and the threshold-form rules are statistically tied with it. Say so plainly in the paper rather than claiming
the binary rule is the truth.
(3) Two fixes to the story: the theorem-to-rule bridge is motivational (see the correction at §29), and labels are agent-specific verdicts with hedged attributions.
(4) The single most valuable next data would be many more labelled targets, each recording the modality; failing that, the registered prospective predictions.

**Scripts:** `target_labels_validation.py`, `validation_test.py`, `validation_addendum.py`, `damage_tissue_map.py`, `independent_tests.py`, `threshold_harm.py`; results in `validation_test.json`, `validation_addendum.json`,
`independent_tests.json`, `threshold_harm.json`; pinned by `tests/test_validation_and_independent.py`.


---

## 39. Tenth addendum (September 29, 2026, second session): are the weights the bottleneck, and what is the most realistic rule?

Prompted by two questions: *"maybe the expendability weights are messing things up"* and *"find the rule you truly believe is the most realistic, even if it is much more complex."* Everything below was
precommitted before it was run (`weights_diagnostic_precommitment.md`, `faers_precommitment.md` + addendum, `relative_exposure_precommitment.md`, each with input hashes). One line of honesty first:
**the alternatives tried against the frozen rule now number 74** (52 before this round, 13 measured-severity variants in §39.4, 9 relative-exposure variants in §39.5), which the methods section must state.

### 39.1 Where the frozen rule's information actually sits

The blood flag alone (8 safe targets flagged, 0 toxic) separates at **0.711**. Among the 22 targets that are not blood-flagged (11 safe, 11 toxic), top share alone separates at **0.851**. Every misordered pair of the frozen rule
(0.914 overall) is inside that 22-target group, so that is where any better weighting would have to show up.

### 39.2 Could any weights help? Upper bound and honest cross-validation (round 6A)

Family offsets (10 free numbers, one per lineage family of `lineage_map.py`) were fitted to the labels with hindsight, which is an optimistic upper bound, and again with leave-one-target-out cross-validation.

| model | in-sample best | optimiser on shuffled labels (mean / 95th pct) | nested leave-one-target-out |
|---|---|---|---|
| frozen rule (blood = 1, others 0) | 0.914 | — | — |
| F: offset for the largest lineage family (10 numbers) | **0.923** | 0.669 / 0.804 | 0.756 |
| G: offset weighted by lineage composition (10 numbers) | 0.928 | 0.706 / 0.837 | 0.799 |
| system-level weights (28 numbers) inside the relative-exposure sum of §39.5 | 0.978 | 0.694 / 0.835 | 0.775 |

**With hindsight the best possible family-level weights add 0.009** (0.914 to 0.923), and the same 16 pairs stay misordered (CEACAM5 7, CA9 5, FOLR1 2, FAP 1, EPCAM 1 on the toxic side; MET 5, NECTIN4 3, MSLN 2, GPC3 2, TNFRSF8 2,
GPRC5D 1, TACSTD2 1 on the safe side). Finer weights fit these 30 targets better (0.978) but do not carry over (0.775 cross-validated, below the frozen rule): they memorise noise. So the answer to "are the weights the bottleneck"
for the safe/toxic question is **no**: the misorderings are between epithelial targets on both sides, and no weighting of lineages can pull them apart.

### 39.3 A measured expendability weight (round 6B)

Hand weights are judgement and lifespan is only a proxy, so a measured weight was built from FDA adverse-event reports (openFDA / FAERS, data last updated 2026-07-30, fetched 2026-09-29; `faers_counts.json`, `faers_severity.py`).
For each of 28 organ systems, the reactions that signal loss of that system's function (2 to 5 MedDRA terms each, fixed in the precommitment) were pooled; **severity = share of those reports that record death**, and weight = severity / largest severity.
Cell types map to systems by identity (147 of 154 mapped; 7, such as adipocytes and fibroblasts, get the median). Because I saw that the death weights disagreed with the hand weights before scoring anything, a second weight was added the same way
(share recording death, life-threatening or disabling), declared in `faers_precommitment_addendum.md` before any score existed.

| system | reports | deaths | weight (death) | weight (severe) | hand weight (mean of member cell types) |
|---|---|---|---|---|---|
| lung | 171,161 | 38.2% | 1.000 | 1.000 | 0.722 |
| heart muscle | 115,163 | 28.9% | 0.756 | 0.759 | 0.850 |
| liver | 94,244 | 28.1% | 0.737 | 0.736 | 0.600 |
| skin | 62,637 | 13.6% | 0.357 | 0.581 | 0.525 |
| brain | 48,537 | 16.1% | 0.421 | 0.566 | 0.783 |
| skeletal muscle | 56,518 | 10.3% | 0.270 | 0.550 | 0.600 |
| platelets | 114,175 | 17.7% | 0.463 | 0.539 | 0.450 |
| kidney tubule | 349,345 | 17.5% | 0.457 | 0.536 | 0.700 |
| marrow stem cells | 110,248 | 16.7% | 0.439 | 0.534 | 0.600 |
| vascular | 23,046 | 14.0% | 0.366 | 0.491 | 0.725 |
| neutrophils | 208,519 | 13.1% | 0.344 | 0.451 | 0.600 |
| colon and gut | 57,420 | 14.7% | 0.384 | 0.429 | 0.586 |
| pancreatic islets | 28,679 | 5.9% | 0.153 | 0.390 | 0.900 |
| red cells | 191,242 | 12.5% | 0.326 | 0.390 | 0.550 |
| adrenal | 13,581 | 7.7% | 0.201 | 0.370 | 0.700 |
| pituitary | 3,614 | 8.3% | 0.218 | 0.366 | 0.629 |
| lymphoid (B, T, other immune) | 36,932 | 10.5% | 0.274 | 0.339 | 0.434 |
| upper GI | 111,719 | 8.9% | 0.233 | 0.322 | 0.606 |
| biliary | 26,646 | 10.7% | 0.279 | 0.312 | 0.700 |
| thyroid | 33,956 | 7.0% | 0.184 | 0.307 | no cell type in the 154 |
| pancreas, exocrine | 73,826 | 7.5% | 0.196 | 0.272 | 0.650 |
| peripheral nerve | 109,046 | 4.8% | 0.126 | 0.272 | 0.700 |
| reproductive | 3,200 | 1.7% | 0.045 | 0.265 | 0.473 |
| kidney glomerulus | 26,389 | 7.3% | 0.192 | 0.248 | 0.900 |
| retina and optic nerve | 169,345 | 2.3% | 0.060 | 0.237 | 0.825 |
| salivary | 81,203 | 4.9% | 0.130 | 0.203 | 0.640 |
| melanocytes (vitiligo) | 2,275 | 2.4% | 0.063 | 0.160 | 0.150 |
| eye surface | 54,795 | 1.0% | 0.027 | 0.119 | 0.700 |

**The measured weights are uncorrelated with the hand weights** (Spearman −0.008 for death, 0.041 for severe; the two measured weights agree with each other at 0.900). The sensible ordering (lung, heart, liver most lethal; vitiligo, dry eye least) is reassuring for
the measurement. The disagreement is a difference of concept: the hand weights encode *irreversible loss* (islets, glomeruli, retina, eye surface all at least 0.7), the report data encode *acute lethality*. Dose-limiting toxicity in the labels includes
severe non-fatal damage (corneal keratopathy, neuropathy), which a death share under-counts by construction; this is why the second weight was added. Caveats: spontaneous reports, fatal cases over-reported, sick patients, and a human choice of which terms belong to which system.

### 39.4 Measured weights in the earlier rule forms (13 variants)

Seven rule forms (A1 to A7) were fixed in advance. A1 = top share + (1 − weight of the home cell type); A2 = top share + (1 − expression-weighted mean weight); A3 = −(worst organ) with absolute expression, thresholds 1, 3, 10, 30, 100 nCPM integrated;
A4 = A3 with a sum; A5 = A3 with the old per-cell-type hand weights (control); A6 = A3 on protein staining (organ level); A7 = RNA and protein together. Weight-shuffled null: the 28 weights permuted 1,000 times.

| rule | weights | separation (95% interval) | paired difference vs frozen | shuffled-weight p |
|---|---|---|---|---|
| A1 | death / severe | 0.727 (0.531 to 0.895) / 0.775 (0.593 to 0.923) | −0.187 / −0.139 | 0.658 / 0.498 |
| A2 | death / severe | 0.823 (0.651 to 0.962) / 0.842 (0.675 to 0.967) | −0.091 / −0.072 | 0.674 / 0.589 |
| A3 | death / severe | 0.708 (0.526 to 0.868) / 0.699 (0.512 to 0.864) | −0.206 / −0.215 | 0.777 / 0.817 |
| A4 | death / severe | 0.789 (0.593 to 0.943) / 0.799 (0.612 to 0.947) | −0.124 / −0.115 | 0.839 / 0.811 |
| A5 (hand control) | hand | 0.799 (0.622 to 0.938) | −0.115 | 0.275 |
| A6 (protein, 23 targets) | death / severe | 0.708 / 0.688 | −0.200 / −0.221 | 0.337 / 0.348 |
| A7 (RNA + protein, 23 targets) | death / severe | 0.721 / 0.704 | −0.188 / −0.204 | 0.633 / 0.658 |

None was promoted (frozen rule 0.914; promotion needed 0.944 and a family-wise p below 0.05). **The absolute-level, worst-organ rules that resemble industry practice (A3, A6, A7) score 0.69 to 0.72 and are clearly below the frozen rule (paired intervals exclude 0; A7, the rule I had declared most realistic in this round, is 0.721 with a paired interval of −0.385 to −0.022).
Only A2 is statistically tied with the frozen rule (paired intervals −0.254 to +0.048 and −0.187 to +0.014). Measured weights are indistinguishable from shuffled weights in every form.** The reason is in the labels: safe and toxic targets alike are expressed in vital organs (TROP2, MET and PSMA as much as HER2 and EGFR); what separates them is how the normal level compares with the tumor level
and how the drug is built, and an absolute cutoff cannot see that.

### 39.5 The theory-faithful form: harm relative to the tumor level (round 6C)

The theorem's harm is measured against the tumor comparator. Per-target tumor levels are not available, so the tumor level was set to the gene's own normal peak, with a floor of 10 nCPM: x_ref = max(peak, 10). A cell type counts as exposed at
r = 0.03, 0.1 and 0.3 of x_ref (the range from cell-therapy sensitivity to antibody-drug-conjugate sensitivity) and exposure is the mean over the three. Harm is the MAX (worst organ) or SUM over cell types of weight × exposure; score = −harm. This form contains
both parts of the frozen rule (the home cell type is always exposed, so its weight enters; more cell types above r × peak add breadth). Nine variants, fixed in advance:

| weights | form | separation (95% interval) | paired difference vs frozen | shuffled-weight p | family-wise p | non-blood targets (frozen 0.851) |
|---|---|---|---|---|---|---|
| binary (blood = 0, else 1) | sum | **0.895** (0.766 to 0.981) | −0.019 (−0.086 to +0.053) | 0.012 | 0.209 | 0.835 |
| hand | sum | 0.876 (0.732 to 0.976) | −0.038 (−0.120 to +0.033) | 0.244 | 0.659 | 0.826 |
| measured severe | sum | 0.861 (0.703 to 0.976) | −0.053 (−0.153 to +0.029) | 0.543 | 0.950 | 0.810 |
| measured death | sum | 0.852 (0.684 to 0.971) | −0.062 (−0.196 to +0.048) | 0.644 | 0.994 | 0.810 |
| hand | max | 0.727 | −0.187 | 0.508 | 1.000 | 0.628 |
| measured death | max | 0.694 | −0.220 | 0.717 | 1.000 | 0.603 |
| measured severe | max | 0.675 | −0.239 | 0.729 | 1.000 | 0.579 |
| binary | max | 0.632 | −0.282 | 0.004 | 1.000 | 0.500 |
| measured death, RNA + protein (23 targets) | max, system level | 0.742 | −0.167 | 0.575 | 1.000 | 0.609 |

**Predictions P1 to P3 were all confirmed.** P1: the binary-sum form is within 0.05 of the frozen rule (difference −0.019, interval includes 0). P2: no graded weight set beats binary in the same formula (0.895 against 0.876, 0.861, 0.852). P3: nothing meets the promotion rule
(family-wise null 95th percentile 0.942). Among the sum forms, only the binary weights differ from shuffled weights (p = 0.012); the graded sum forms have shuffled-weight p of 0.24 to 0.64, i.e. indistinguishable from arbitrary weights. The floor T0 changes nothing at 3, 10 or 30 nCPM because no panel gene has a peak below 18.6 nCPM (CD79B), so this panel cannot test the floor (it exists for tumor-specific genes with almost no normal expression, which are absent here). Compared with §38.4, where the same sum with an absolute threshold scored 0.80 to 0.86, measuring exposure relative to the gene's own peak lifts the binary form to 0.895. The worst-organ (max) form is poor because it ignores breadth, and outcomes track breadth; the SUM form, which is the theorem's own, is 0.15 to 0.26 better.

### 39.6 Do graded weights at least put the damage in the right place? (naming test, same metric as §38.3)

Best rank (of 154 cell types) of the documented-damaged tissue when each cell type is ranked by weight × relative exposure, 13 toxic targets (exploratory; precommitted metric, not a new rule):

| weights | median best rank | damaged tissue in top 3 | vs no weights (better / worse) | sign-test p |
|---|---|---|---|---|
| none (exposure only) | 9.0 | 2 of 13 | — | — |
| binary | 8.0 | 2 of 13 | 2 / 1 | 1.000 |
| hand | 4.5 | 5 of 13 | 9 / 4 | 0.267 |
| measured (death) | 4.0 | 4 of 13 | 7 / 4 | 0.549 |
| measured (severe) | 5.0 | 3 of 13 | 7 / 4 | 0.549 |

Graded weights, hand or measured, move the damaged tissue from about rank 9 to about rank 4 to 5, and binary weights do not; none of it is significant with 13 targets. It is the same direction as §38.3, and it is the one place where graded weights look more realistic.
The failures are informative: eye surface has the lowest measured weight of all 28 systems (0.027 for death, 0.119 for severe), so F3 (ocular boxed warning) drops from rank 8 with no weights to rank 49.5 (death) and 41.5 (severe).

### 39.7 What this means, in plain terms

(1) **The weights are not what limits the safe/toxic score.** Best-possible family weights add 0.009; fitted system weights do not generalise (0.775); independently measured weights are no better than shuffled ones; graded weights lose to binary in the same formula.
(2) **Realism is not the same as predictive power here.** The rule I had pre-declared as most realistic (worst-organ, relative exposure, measured death weights; §39.5) scores 0.694, clearly below the frozen rule, with a paired interval (−0.397 to −0.053) that excludes 0. The rule I had declared most realistic in the previous round (A7, absolute level, worst organ, protein and RNA) also scored lower (0.721). On this panel
outcomes follow breadth and lineage more than the lethality of the single most exposed organ, so my a-priori belief in a worst-organ aggregation was wrong here, and the theorem's own sum is the better form.
(3) **The frozen rule now has a theory-faithful restatement.** The threshold form of the theorem with binary weights and a tumor-proxy threshold gives 0.895 without using top share at all, as predicted before it ran. That closes part of the gap opened by the correction at §29: confinement follows from summing over
cell types above a threshold measured against the tumor level. Two cautions: it is one of 74 alternatives, and its family-wise p is 0.209, so present it as a concordance check, not a discovery.
(3b) Where graded weights help is *where* the damage lands (§39.6), not *whether* a target is dose-limiting. A two-part description is therefore honest: whether a target is dose-limiting is predicted by confinement and lineage; where it lands is described better by graded weights (not significant on 13 targets).
(4) **Why the safe/toxic labels resist a normal-tissue rule.** They are agent-specific verdicts about dose-limiting toxicity, which depends on modality and on the tumor-to-normal ratio, and includes non-fatal damage. What the rule lacks is a per-target tumor level and modality, not better organ weights.
(5) The frozen rule stays the primary, pre-specified rule. State the count of alternatives (74).

**Scripts:** `weights_diagnostic.py`, `faers_severity.py`, `faers_rules.py`, `relative_exposure.py`, `relative_naming.py`; results in `weights_diagnostic.json`, `faers_counts.json`, `faers_rules.json`, `relative_exposure.json`, `relative_naming.json`;
pinned by `tests/test_faers_and_relative.py`.

---

## 40. Eleventh addendum: how fragile is the result to wrong labels, and a proved link between top share and exposure breadth

Two checks that add no new rule (the count of alternatives stays 74). Both were run after the Round 7 finding that some labels are weakly supported. Precommitted in `label_fragility_precommitment.md` (hash-locked before the run); scripts `label_fragility.py`, `breadth_bound_lemma.py`; pinned by `tests/test_fragility_and_bound.py`.

### 40.1 Label fragility: how many wrong labels would erase the result?

The rule and its scores are unchanged (top share plus the lineage-level blood flag, separation 0.914 on 19 safe v 11 toxic). Only labels are swapped (safe to toxic or toxic to safe), in every possible combination. The p-value is exact (the 30 scores have no ties, so the exact Mann-Whitney distribution is used; for the frozen labels it gives 2.8e-5, against the stored permutation estimate of 3.5e-5). This replaces the normal approximation named in the precommitment, which was not needed.

| labels swapped (k) | sets tried | lowest separation | median separation | share of sets still p < 0.05 | worst p in any set |
|---|---|---|---|---|---|
| 1 | 30 | 0.833 | 0.880 | 100% | 0.001 |
| 2 | 435 | 0.765 | 0.855 | 100% | 0.007 |
| 3 | 4,060 | 0.704 | 0.824 | 100% | 0.032 |
| 4 | 27,405 | 0.643 | 0.796 | 98.7% | 0.100 |
| 5 | 142,506 | 0.583 | 0.769 | 92.1% | 0.232 |
| 6 | 593,775 | 0.526 | 0.742 | 80.5% | 0.416 |

**Fragility index = 4:** no set of 3 swapped labels removes significance; significance is first lost at 4 swaps (349 of the 27,405 four-swap sets). **All 349 of those sets move at least one of the eight approved blood-cancer targets (CD19, CD22, CD33, CD38, CD79B, CD20, SLAMF7, BCMA) from safe to toxic**, which no reviewer has proposed; the single worst set swaps CD47, ERBB2, CD20 and BCMA. If 4 of the 30 labels were wrong in a way unrelated to the scores, 98.7% of possible sets would still leave p < 0.05.

**The specific weak labels (declared before the run):** EGFR, EPCAM, NCAM1, FAP, CD47, ERBB2 and CEACAM5 swapped to safe, NECTIN4 swapped to toxic, all 256 subsets. Separation ranges from 0.672 to 0.957 and 97.3% of subsets keep p < 0.05. Any 4 of the 8 swapped together leave p at most 0.043; with 5 or more the worst subset reaches p = 0.079 to 0.126. Swapping **all eight** gives separation 0.744 (p = 0.048), that is, borderline significant even if every one of the eight were wrong (the groups become 25 safe v 5 toxic, so the test is weaker as well). Swapping CEACAM5 alone raises separation to 0.950.

Precommitted predictions: P1 (1-2 swaps keep separation at or above 0.80 and p < 0.05) was half met: p stayed below 0.05, but the lowest separation with 2 swaps was 0.765. P2 (index at least 4) and P3 (all-eight swap drops below 0.75) were met.

What this does and does not show: it measures how many mistakes the panel can absorb; it does not say which labels are wrong, and it assumes the taxonomy (dose-limiting on-target toxicity) is the right axis.

### 40.2 A proved link between top share and exposure breadth

Let a gene have expression x_c in N cell types, peak = max x_c, top share t = peak / total. Call a cell type exposed at level r if x_c is at least r times the peak (the top type always is), and let E(r) be their number.

**Lemma 1.** E(r) is at most 1 + (1/t - 1)/r. *Proof:* the total is at least peak + (E - 1) r peak, so 1/t is at least 1 + (E - 1) r. **Lemma 2.** t is at least 1 / (E + (N - E) r) (each unexposed type is below r times the peak). **Tightness:** for a given t and r the bound on E is reached, up to rounding down, by a profile with E - 1 types exactly at r times the peak and the rest spread below it (checked by construction, 2,000 cases). So if only t is known, the largest number of exposed cell types that can be excluded is exactly this bound: top share is the sufficient statistic for the worst-case breadth.

**Consequence for the theorem's threshold form:** with weights between 0 and 1, the summed harm at level r is at most E(r), which falls as top share rises. Ordering genes by top share is therefore ordering them by a guaranteed ceiling on threshold-form harm. This turns the correction at §29 into a partial repair: confinement is not the leading term of the share-weighted formula, but it is the tightest guarantee on the breadth term of the threshold form. It is one-directional: a high top share guarantees limited exposure; a low top share does not guarantee wide exposure.

Numerical check: 0 violations of either lemma in 132 real gene-level checks (44 genes, r = 0.03, 0.1, 0.3) and 0 in 20,000 random expression vectors. On the real genes the bound is loose, as expected of a worst case: the actual number of exposed types is a median of 19%, 28% and 30% of the bound at r = 0.03, 0.1, 0.3, and the bound is below the 154-type ceiling (so informative) for 19, 26 and 30 of the 30 panel genes. But top share tracks actual breadth closely: Spearman correlation of top share with the number of exposed types is -0.95, -0.95, -0.92 at r = 0.03, 0.1, 0.3. At r = 0.1 the median bound is 29 exposed types for safe targets and 125 for toxic targets.

**Scripts:** `label_fragility.py`, `breadth_bound_lemma.py`; results `label_fragility.json`, `breadth_bound_lemma.json`; pinned by `tests/test_fragility_and_bound.py`.

## 41. Locked-predictions check (Sept 29, 2026): a scoping pass

**Set-up.** `prospective_predictions.json` (locked 2026-09-08) made 32 calls (14 TOLERATED, 18 ON-TARGET TOXICITY) plus 20 no-calls on targets in registered trials. The labelling rule, the A/B/C evidence classes and the no-rescoring rule were written and hashed before any search (`locked_predictions_check_precommitment.md`). Four literature agents searched PubMed abstracts and, for approved drugs, the openFDA label; five key abstracts were re-checked by hand.

**Class A (genuinely forward): none.** No called target has a verdict first public after 2026-09-08. This was expected three weeks after the lock. The only post-lock item, a B7-H3 ADC phase 3 (CD276, 2026-09-12), concerns a no-call and reports the tolerability of a payload-limited ADC, so it is not a verdict.

**Class B (evidence existed before the lock but the rule never used it): 5 of 32 called targets scored.**

| target | locked call | evidence | result |
|---|---|---|---|
| CD70 | TOLERATED | CTX130 CAR-T, 16 patients, no dose-limiting toxicity (2024); ALLO-316 DLTs tied to lymphodepletion (2026) | hit |
| LRRC15 | TOLERATED | ABBV-085, safe at 3.6 mg/kg (2021); target is mostly tumour stroma, so a weak test | hit |
| ITGB6 | TOXIC | sigvotatug vedotin, acceptable safety, no organ attributed (2026) | miss |
| ROR1 | TOXIC | zilovertamab vedotin, manageable safety, but 1 responder and ROR1 in 3 of 17 tumours: a poor test (2025) | miss |
| LY6E | TOXIC | DLYE5953A, no dose-limiting toxicity from 0.2 to 2.4 mg/kg (2020) | miss |

**Score: 2 hits, 3 misses of 5.** Hits and misses are reported at equal weight. Five is below the precommitted floor of 8 for drawing any conclusion, so this is not evidence for or against the rule.

**Why so few could be scored.** 13 of the 32 calls are HEDGED and 14 more have no verdict. Almost all of these targets have only been tested as vedotin or deruxtecan ADCs, where toxicity is dominated by the payload, and authors seldom say whether a toxicity is on-target. This is a measurement limit, not a rule failure: a payload-limited ADC can stop before the dose that would reveal on-target damage. The three misses share the same weakness, since all are tolerated ADCs with unattributed toxicity. I did not remove them for that reason after seeing them, because the precommit defined no evidence tiers. If tiers were defined afterwards, only LY6E (no DLT at all) would be a clean "safe" evidence point.

**No-call targets (informational, not scored).** STEAP1, MUC16 and CLDN6 look tolerated; FGFR2 (bemarituzumab, corneal toxicity, inferred) looks toxic; ST8SIA1 counts as toxic only through the GD2 antigen (dinutuximab boxed warning), a gene-to-antigen mapping I would not score.

**Overlap.** GUCY2C, TFRC, SLC44A4 and CD276 (also in the validation set) and AXL, ERBB3 and F3 (also in the main panel) were all HEDGED or unresolved, so nothing here is double-counted.

**What this means for the paper.** Say plainly that the September lock exists, that no forward test has resolved, and that a pre-lock check on the 5 scoreable targets went 2 of 5. Present the lock as a standing design for future resolutions, not as a result. No label enters the frozen panel from this check.

**Limits.** Abstract-level evidence only; full texts were often blocked. Registries were not reachable. LYPD3, CEACAM6 and DPEP3 had thin searches. Files: `locked_predictions_check.json`, `locked_predictions_check_precommitment.md(.hashes)`.
