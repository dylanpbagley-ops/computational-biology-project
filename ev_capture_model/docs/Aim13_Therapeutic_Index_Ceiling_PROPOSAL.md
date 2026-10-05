# Proposal: Aim 13-15 — the therapeutic-index ceiling

**Date: September 6, 2026. Status: proposal, with every load-bearing claim already checked against real code and real data in this session. Nothing below is speculative unless it is explicitly labelled as a next step.**

This document proposes the next phase of the project. It has three parts:

1. **A problem** — two real, unaddressed confounds in Aim 12, the project's current strongest result, found this session by re-deriving Aim 12's numbers from the local raw files rather than trusting the handoff.
2. **A simplification** — an exact algebraic identity that collapses the model's entire parameter grid to a single, physically meaningful number, verified bit-for-bit against the existing code.
3. **A new claim** — a genuinely new, falsifiable, previously-unstated result that the identity makes possible, already run end-to-end on the 41 real patients and validated against real clinical ground truth.

The three are connected: the fix for the problem is the same operation that produces the new claim.

---

## Part 1 — two real confounds in Aim 12

Aim 12 computes, for each of 41 real DLBCL patients, a per-cell CD19 distribution and ranks patients by it. Every number in `REAL_PATIENT_CD19_SUMMARY` was independently reproduced this session from the local FCS files (P5904 mean 62.701, P6216 29.518, P6389 101.409, P7676 109.106, O6357 348.50 — all exact). The parser and the pipeline are correct. The problem is what the numbers are computed over.

### Confound A: the distributions are whole-biopsy, not tumor

`load_marker_values_from_fcs_file` returns every event in the file. `$TOT` for P5904 is 47,568 and Aim 12 records `n_cells = 47568`, so no gating was applied anywhere. A DLBCL biopsy is not a tumor-cell suspension — it is tumor B cells plus tumor-infiltrating T cells plus myeloid cells plus debris.

Measured across all 41 files (CD19/CD20-positive, CD3-negative for B; CD3-positive for T):

| | min | max |
|---|---|---|
| B-cell fraction | 0.050 (R6589) | 0.868 (Q5971) |
| T-cell fraction | 0.11 (Q4644) | 0.93 (R6589) |

So a patient's "CD19 distribution" is a mixture whose composition varies 17-fold across the cohort. The consequences are concrete:

- Spearman correlation between the ungated per-patient mean CD19 and the B-cell-gated mean CD19 is **0.825** — correlated, but **15 of 41 patients move more than 5 rank positions**, and one moves **18**.
- The ungated mean CD19 correlates with B-cell fraction at rho = **0.405**: part of what the current ranking measures is tumor purity, not CD19 density.
- The documented worked example is affected. Aim 12's top-5-by-naive-mean is `O6357 > Q7998 > Q5971 > O8155 > Q7915`; this session reproduced that exactly from the ungated data. After gating to B cells it becomes `Q7998 > O6357 > R5097 > Q7915 > P7901`. In particular, the specific pair the write-up highlights — P7676 (mean 109.1) versus Q7915 (mean 111.8), where the point is that P7676 has the *lower* mean — reverses under gating: Q7915 is 227.0 and P7676 is 139.0.

The 59.3% reordering finding is still a real finding. But "the naive mean and the distribution-aware score disagree" currently rests on a naive mean that is partly a measure of immune infiltrate. A judge who asks "were these gated?" needs a better answer than the current one.

### Confound B: cross-patient signal scale varies ~10x by batch

The 41 samples were run in 5 barcoded batches. Median 193Ir-DNA signal — a standard internal reference for overall CyTOF signal intensity — by batch:

| batch | n | median DNA |
|---|---|---|
| CytB2 | 9 | 34.0 |
| CytB3 | 9 | 23.9 |
| CytB4 | 9 | 234.6 |
| CytB5 | 12 | 288.4 |
| CytB6 | 2 | 301.5 |

That is a clean ~10x split between the first two batches and the last three, and a 22x range overall (19.1 to 421.4). Batch identity alone explains **19.2%** of the variance in log ungated mean CD19 (10.4% after gating). Any *cross-patient* comparison of absolute CD19 on this data is partly a comparison of run conditions.

Neither confound is fatal and neither was hidden — nobody had checked. Both are fixable, and the fix leads directly to something better.

---

## Part 2 — the identity: escape fraction is a CDF evaluated at one point

The model's kill probability is

```
theta(n) = z*q / (1 + z*q),    q = exp(n * sigma * K) - 1
```

and `escape_fraction` counts cells with `theta(n) < tau`. Because `x -> x/(1+x)` is strictly increasing on `x >= 0` and `exp` is strictly increasing, that inequality solves exactly:

```
theta(n) < tau   <=>   n < n*,   where   n* = ln(1 + tau/(z*(1 - tau))) / (sigma * K)
```

**Therefore `escape_fraction(population) = F(n*)`, the empirical CDF of the population's antigen counts evaluated at the single point n\*.** Verified numerically against the real P5904 data across 60 combinations of `k_intra` (1e-5 to 1e-2), `tau` (0.1 to 0.9) and `z` (0.5 to 2.0): the computed escape fraction and the computed `fraction(n < n*)` agree to better than 1e-12 in every single case, with no exceptions.

Three consequences.

**1. The 27-point grid is one-dimensional, not two.** `k_intra` and `kill_threshold` are not independent axes; each pair just picks a value of n\*. Sampling 27 points of a 1-D curve is unnecessary — the whole curve is computable exactly.

**2. n\* is not an abstract parameter, it is the field's own quantity.** At `tau = 0.5` the formula gives n\* = ln(2)/(sigma*K), which for the three real Aim 10 calibrations returns exactly 81, 843 and 8,690 molecules/cell — the measured activation thresholds themselves. n\* *is* the antigen-density threshold for CAR/TCE triggering, a quantity that is directly measured in the literature. The 27-point grid spans n\* from about 12 to about 28,900 molecules/cell.

This is a substantial upgrade in how the existing result can be stated. "59.3% of a 27-point grid reorders" invites "why that grid?" The identity replaces it with: *ranking depends on the agent's antigen-density threshold; here is the exact set of thresholds at which the ranking changes; here is where real, measured clinical agents sit relative to them.*

**3. It names the real criterion.** Ranking two populations by escape fraction is ranking them by `F_A(n*)` versus `F_B(n*)`. Mean-based ranking agrees with the model at *every* n\* if and only if one CDF lies entirely below the other — first-order stochastic dominance. When the CDFs cross, mean-based ranking is provably wrong for some n\*, and the crossing points say exactly which. The project's empirical "% of grid reorders" is a coarse sample of a sharp, classical criterion.

---

## Part 3 — the new claim: a computable ceiling on any target's therapeutic index

### The idea

Escape depends on the tumor distribution through `F_T(n*)`. Every real agent also acts on normal cells that carry the same antigen — that is on-target/off-tumor toxicity, the field's central problem outside B-cell malignancies. Toxicity is `1 - F_N(n*)` for the normal-cell distribution.

As n\* sweeps from 0 to infinity, the pair `(1 - F_N(n*), 1 - F_T(n*))` traces exactly the **ROC curve** of antigen density as a classifier of tumor versus normal cells.

This gives a theorem-shaped statement:

> For any agent whose per-cell kill probability is a monotone increasing function of target antigen count, the achievable (tumor killed, normal spared) operating points lie on or below the ROC curve of the tumor and normal per-cell antigen-density distributions. Consequently the maximum achievable therapeutic index for a given target is a fixed, computable property of that antigen's tumor-versus-normal distribution overlap, and cannot be improved by any amount of affinity tuning, CAR-density engineering, or costimulatory redesign that leaves the kill function monotone.

The engineering knobs the field spends most of its effort on — scFv affinity, CAR surface density, hinge/spacer, costimulatory domain — move n\*. They move you *along* the curve. They do not move the curve. Where the curve itself is bad, the target is bad, and no construct fixes it.

That is a genuinely new, decision-relevant output: **an upper bound on how good any therapy against a given target could ever be, computable in advance from single-cell data, before anyone builds a construct.**

It also produces a clean actionable pairing: the ROC gives the ceiling, and the measured n\* of a real agent gives the operating point actually reached. The gap between them is the part that *is* fixable by engineering.

### Why this fixes Part 1's confounds rather than inheriting them

- **Batch scaling cancels.** ROC and AUC are invariant under any strictly monotone transform applied to a channel. Tumor and normal cells being compared come from the same file, same run, same staining, so a batch scale factor applies identically to both and cancels exactly.
- **The molecule-calibration gap largely evaporates.** Aim 12's most-cited limitation is that the CyTOF signal is not confirmed molecules-per-cell. For a *ceiling* stated as an AUC, that does not matter — AUC needs only the ordering of cells, not their units. Absolute calibration is still needed to place the *operating point* of a specific real drug, which should be stated as such.
- **Gating is required by construction,** so Confound A gets fixed on the way in rather than left open.

### It has already been run, on the real data, and it works

The 41 local FCS files carry a 40-antibody B-cell panel, not just CD19. Present in it are **CD19, CD20, CD22, CD79B, CD38, CD40, CD21, CD72, CD23, CD24, CD25, CD27, CD43, CD44, CD48, CD80, CD83, CD84, CD86, CD49D, CD10, CD45RA, HLA-DR, HLA-ABC, HLA-E, CXCR4, CXCR5, CCR4, PD-L1, FasR, TRAIL-R1, IL-4R, BTLA, IgM, IgD, IgKappa, IgLambda** — 37 usable antigens measured on the same cells. It also carries **CD3 and CD45**, so a non-B leukocyte population is available inside every sample as a matched bystander comparator.

Pilot run, all 41 patients, QC-gated (intact by DNA, live by 195Pt, bead-excluded, CD45-positive), tumor = CD19/CD20-positive CD3-negative B cells, bystander = CD3-positive T cells, up to 20,000 cells per group per patient, AUC by exact rank statistic:

Ranking the 37 antigens by median tumor-versus-bystander AUC gives, top-8: **CD19 (0.976), CD20 (0.972), HLA-DR (0.931), CD72 (0.928), CD40 (0.924), CD22 (0.918), CXCR5 (0.912), CD21 (0.875)**.

Ranking the same 37 antigens by median tumor mean expression — a stand-in for what target-selection practice actually does — gives top-8: **HLA-DR, CD44, HLA-ABC, CD19, CD20, CD38, CD22, CD72**.

Spearman correlation between the two rankings: **0.372**. They substantially disagree.

**Test against real clinical ground truth.** Ground-truth labels were fixed before looking at the ranks. Validated = an antigen with an approved immunotherapy in B-cell lymphoma: **CD19** (axi-cel, tisa-cel, liso-cel, blinatumomab), **CD20** (rituximab, mosunetuzumab, glofitamab, epcoritamab), **CD22** (inotuzumab ozogamicin, CD22 CAR-T), **CD79B** (polatuzumab vedotin). Broad = antigens known to be widely expressed on normal leukocytes and therefore expected to carry on-target/off-tumor risk: CD44, HLA-ABC, CD48, CD27, CD43, CD38, CD84, HLA-E, CD25, FasR, CD45RA.

| | overlap (AUC) ranking | mean-expression ranking |
|---|---|---|
| mean rank of the 4 validated targets (of 37) | **5.0** | 8.5 |
| mean rank of the broad/risky antigens (of 37) | **29.9** | 17.4 |
| all 4 validated ranked above all broad ones | **yes, 1.000** | no, 0.705 |

Perfect separation of 4 validated targets above 11 broad ones has probability 1/C(15,4) = **7.3e-4** under random ranking.

The individual cases are the interesting part:

- **CD38** — mean rank 6, AUC rank 28, median AUC **0.538**. CD38 is nearly useless for distinguishing tumor B cells from T cells in this data. CD38 CAR-T is well known to suffer T-cell fratricide because T cells express CD38. The model recovers the mechanism of a real, published clinical failure from data alone.
- **CD44** — mean rank 2, AUC rank 31, AUC 0.475 (no discrimination). Pan-leukocyte adhesion molecule; CD44v6-directed therapy has a real toxicity record.
- **HLA-ABC** — mean rank 3, AUC 0.544. MHC class I, on every nucleated cell. A mean-expression rule puts it third; an overlap rule correctly discards it.
- **CD27, CD48, CD43** — AUCs of 0.171, 0.312, 0.218: *inverted*, higher on T cells than on tumor. These are the natural candidates for a NOT-gate, which connects cleanly to Aim 5's existing AND-gate machinery.
- **CD79B** — mean rank 18 but AUC rank 11. A validated ADC target that mean-expression ranking buries below CD44, HLA-ABC and CD38.
- **CD72 and CD40**, ranked 4th and 5th by overlap, are both real active development targets rather than noise — CD72 CAR-T in particular is in development for B-cell malignancies.

**Robustness.** Repeating the whole pipeline with lineage-gate cutoffs at the 40th, 50th, 60th and 70th percentile: validated-above-broad separation stays **1.000** at every setting, and CD19/CD20 stay ranked 1 and 2 at every setting. The result is not an artifact of one gating choice.

### Why this is the thing the project has been missing

The project's own novelty diagnosis concluded that the winning pattern is *validate a model against independent published data, then use the validated model to produce one specific new claim*, and that this project had never done the second half — every result was a sensitivity or reordering finding checked against data whose answer already existed. This closes that gap on both counts at once:

- It produces a **new, previously-unstated claim** (a computable ceiling on achievable therapeutic index per target, and a ranking of 37 real antigens by it) rather than another consistency check.
- It is **validated against genuinely independent, external ground truth** — which antigens actually produced approved drugs — and the model's ranking recovers that ground truth while the field's default heuristic does not.

It also, incidentally, supplies the healthy-versus-disease comparison at single-cell resolution that the diagnosis document named as Weak point 2, and the head-to-head multi-marker ranking comparison named as Weak point 4. Both were listed as needing new data that had never been found. Both are answerable from files already on disk.

---

## Honest limitations, stated now rather than discovered in January

1. **The comparator is currently bystander T cells, not normal B cells.** This is the biggest gap. For B-lineage-restricted antigens, T cells are an easy comparator, and CD19/CD20/CD22 score near-perfectly against them — but the real dose-limiting toxicity of CD19 CAR-T is B-cell aplasia, i.e. normal B cells. Against a normal-B comparator, CD19 should score *poorly*, and it should, because every patient on CD19 CAR-T loses their normal B cells. That would be a strong retrodiction of a universally known clinical fact, and it is currently untested. Two routes, both real: (a) the same Mendeley release contains **31 reactive lymph node (non-malignant) control samples** that have never been downloaded — this session's network is blocked from Mendeley, so getting them needs a browser session; (b) residual polyclonal normal B cells inside each tumor sample, identified by light-chain restriction using the panel's IgKappa and IgLambda channels — checked and clearly present (sample O6357 shows a kappa:lambda ratio of 416:1 among B cells, textbook monoclonality, with a measurable polyclonal remainder), though surface-Ig-negative DLBCL will limit how many patients this works in. Route (b) is scientifically better where it works, because the control is same-patient and same-tube. Both should be attempted.
2. **The monotonicity assumption is load-bearing, and this project already owns a counterexample.** The theorem requires kill probability to increase monotonically with antigen count. The project's own CLDN18.2 pancreatic-line finding is a real measured non-monotonicity. That should be stated as the honest boundary of the theorem, not buried; owning a counterexample from your own data is stronger than pretending none exists.
3. **The gating is mine, not an expert's.** Percentile-based lineage gates are defensible and demonstrably robust across a 40-70th-percentile sweep, but they are cruder than clustering or expert manual gating. A proper version should compare at least two gating strategies.
4. **The ground-truth labels are a judgement call.** "Validated" and "broad" were fixed before looking at ranks, but the label set is small (4 and 11) and one case is arguable — CD45RA is high on naive B cells, so calling it purely "broad" is not clean. The claim needs a pre-registered label list, ideally with the labels assigned by criteria written down in advance and, if possible, checked by someone else.
5. **AUC is a ceiling, not a prediction of clinical success.** A target can have a great ROC and fail for reasons this framework says nothing about — trafficking, exhaustion, the microenvironment, antigen shedding, manufacturing. The claim is a necessary condition, not a sufficient one, and must be phrased that way every single time.
6. **The novelty check is preliminary, not complete.** The nearest prior art found so far is Dannenfelser et al., *Cell Systems* 11:215-228 (2020), which screens Boolean AND/NOT antigen combinations against bulk tumor and normal-tissue transcriptomes — same goal (tumor-versus-normal discrimination for CAR targets), different data (bulk RNA, not per-cell surface protein), different method (binary logic gates, not a continuous dose-response threshold), and no ceiling statement. Also relevant: single-cell on-target/off-tumor expression atlases, and the antigen-density-threshold literature (e.g. functional-avidity work establishing that constructs have measurable density thresholds). None found so far states the ROC-as-achievability-bound result or the stochastic-dominance ranking criterion. **A systematic novelty search must be completed before this is claimed as new** — this is the single highest-priority thing to do before building anything.

---

## Proposed structure

**Aim 13 — the identity and the criterion.** Prove `escape = F(n*)`, prove the stochastic-dominance criterion, replace the 27-point grid with the exact continuum, and re-express the existing Aim 11/12 results in terms of where real measured agent thresholds fall. Small, mostly analytic, and it strengthens work that already exists.

**Aim 14 — gated re-analysis of Aim 12.** Redo Aim 12 on properly gated tumor cells with an explicit batch-effect analysis. Report honestly how much of the original finding survives. This is the integrity move, and a judge who spots the confound independently will be far more impressed that it was caught and fixed than that it was never mentioned.

**Aim 15 — the therapeutic-index ceiling.** The 37-antigen ranking, both comparators (bystander leukocytes and normal B cells), the ROC ceiling per target, the operating point of real agents from Aim 10's calibrations, the ground-truth validation, and the NOT-gate candidates that fall out of the inverted antigens.

Rough sequencing over the ~6 months available: novelty search first (it can kill or reshape the whole thing, so it goes first); then Aim 13, which is cheap; then Aim 14; then Aim 15, with the normal-B comparator as its main data dependency; then the paper and poster with roughly two months of slack. The pilot in this document took under an hour of compute on data already on disk, which is a reasonable indication that the remaining work is a few sessions rather than a new project.

---

## What is needed from Dylan

- A decision on whether to pursue this at all, and at what scope.
- A browser session at some point to fetch the 31 reactive-lymph-node control files from the Mendeley release (this session's network cannot reach Mendeley from either the sandbox or the local shell).
- A read of the honest-limitations list above, particularly items 1, 2 and 6, before any of this gets written into a claim.
