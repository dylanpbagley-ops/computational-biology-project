# Project Handoff v5: Computational Model of Multivalent EV Capture

**Supersedes:** `EV-Capture-Project-Handoff-v4.md` (v4, August 2026). Where the two disagree, this document wins. Key changes are listed in "What changed since v4" below.

**Purpose:** hand a new chat everything it needs to continue this project without re-deriving anything or repeating dead ends. Read "What NOT to revisit" carefully — this project reached its current form after three prior conversations of failed shapes, a fourth (v2→v3) that finished the figure-digitization step v2 had left blocked, a fifth (v3→v4) that resolved v3's top-priority open question (n vs. areal receptor density), and a sixth (v4→v5) that tackled v4's top-priority open question: Aim 2 leg 3, external validation against real published signatures.

---

## Who I am and how I work

High school student, working solo, ISEF-level science fair project in **computational biology**. Regional fairs are typically February–March. Working timeline is roughly six months from late August 2026.

How to work with me:

- **Plain language, minimal jargon.** Define terms the first time. Unfamiliar vocabulary actively hinders my understanding; it doesn't impress me.
- **Explain fully.** Concise is good, but not at the cost of leaving something out. If it needs saying, say it.
- **I'm a beginner coder.** Python basics only. Walk me through both the physics and the debugging. Don't assume I know what a class is or how to structure a project.
- **Stress-test before executing.** Evaluate from multiple angles and flag problems before implementing. I've been burned by premature commitment three times.
- **Be honest about weaknesses.** I'd rather hear "this could fail because X" now than find out in January. Don't reassure me reflexively.
- **I catch errors and expect them acknowledged.** I've caught conflated numbers, unnecessary anchoring on a cell line, and overstated claims. Acknowledge corrections directly rather than glossing over them.

---

## The project in one sentence

Affinity capture of extracellular vesicles is **nonlinear** in per-vesicle receptor count, so achievable tumor-vs-normal selectivity is governed by the **overlap of two receptor distributions**, not by the fold-change of their means — which means the standard way the field picks capture targets can rank them incorrectly.

**Extracellular vesicles (EVs)** are particles roughly 30–150 nm across that cells release into biofluids, carrying proteins from their parent cell on their surface. "Affinity capture" means coating a surface with molecules that grab those proteins, to pull tumor EVs out of a blood sample — the basis of a liquid-biopsy cancer diagnostic.

---

## What changed since v4

v4's blocking item — "Immediate next steps" item 1, Aim 2 leg 3 (reproduce known external experimental signatures) — has been substantially advanced, along with an honest complication that surfaced immediately and was resolved by asking rather than assuming.

This session worked from a handoff document (v4) that referenced two folders, `ev_capture_model/` (the code) and `ev_capture_data/` (the digitized Lennon figures). Both were attached from the start this time — the v4 handoff's own closing note (to attach both folders up front) was followed, and it worked: no round-trip was lost re-requesting a missing folder, unlike the v3→v4 session.

Five substantive updates:

1. **A real physics mismatch was found and flagged before any code was written, not glossed over.** The DNA nanostar "optimal valency = 3" validation target (Linne et al. 2024, *Soft Matter*, DOI 10.1039/d4sm00704b) turns out to depend on mechanics our EV model does not have: a fixed, small number of discrete binding arms (valency `k`, a real design parameter) meeting *variable* receptor density, where the optimum arises from a *surface-saturation timing* trade-off — more arms means faster binding from solution, which saturates the system before the interesting superselective transition becomes visible. Our EV model has the geometry flipped (a vesicle with variable receptor count meeting *fixed* ligand density) and has zero saturation/depletion mechanics at all — every vesicle captures independently in our model, with no cap on total capture. Relabeling our existing model's variables would not have tested this mechanism at all. This was surfaced directly and a decision was asked for (build the real mechanics vs. skip vs. do a lighter check) rather than silently picking one path.
2. **The real mechanics were built, as a new, separate, honestly-scoped module.** `src/nanostar_validation.py` implements Linne et al.'s own published equations (their Eq. 2–4: a monovalent Langmuir isotherm for k=1, and an avidity extension for k>1 with two equilibrium constants, K_A for the first bond and K_intra for subsequent bonds) directly — not a repurposing of `binding_model.py`. Getting a working reproduction took real iteration: the first parameter guess produced no optimum at all (selectivity was highest at k=1, a sign the model wasn't in the genuinely dilute regime the theory requires), and fixing that surfaced the same "must be dilute enough to see the effect" condition that turned out to matter for our own `selectivity_alpha()` too (see point 4). The paper's own *exact* fitted numbers (their Figure 2e table) were not machine-readable from either the PMC full-text API or the RSC-hosted PDF — both were checked directly. Rather than inventing numbers and presenting them as the paper's own, representative parameters were chosen to match the paper's *stated qualitative relationship* (K_A roughly constant across valency, K_intra decreasing with valency) and the result is reported honestly throughout as "mechanism-faithful, not digit-faithful."
3. **The reproduction works, with two honest limitations surfaced and locked into the test suite rather than hidden.** At representative parameters (K_A=0.01, K_intra(k) ∝ k⁻³·⁵), valency k=3 achieves the highest peak selectivity among k=1,3,6,10 (α=1.49 vs. 1.35 for k=6 and 0.69 for k=10), qualitatively matching the paper's reported experimental optimum. Two limitations found and reported directly rather than smoothed over: (a) this specific optimum is **not robust to K_A** — a ten-fold weaker K_A shifts the winner to k=6 instead, consistent with (not contradicting) the paper's own point that valency, binding strength, and concentration jointly determine the optimum; (b) at these representative parameters, k=10 doesn't exceed the monovalent Langmuir ceiling at all (α=0.69 < 1) — checked systematically, there is no tested K_intra(k) decay rate that makes k=3 win *and* keeps k=10 above the ceiling simultaneously. Both limitations are now regression-tested, not just mentioned in a comment.
4. **`binding_model.selectivity_alpha()` was exercised for the first time since it was written in Aim 1, and this surfaced a genuinely new, previously-undocumented finding about the existing model.** At z=1.0 — the activity/concentration value used in *every* existing demo script in this codebase (`compare_analytical_vs_mc.py`, `demo_population_effect.py`, `patch_regime_sensitivity.py`) — α never exceeds 1 anywhere across a wide receptor-count sweep, meaning none of the project's existing demo figures are technically showing superselective behavior by the strict Martinez-Veracoechea/Frenkel definition (α > 1), even though the underlying exponential nonlinearity is still doing real work in those figures. This is not a bug — it traces to the same "must be in a sufficiently dilute regime" condition discovered independently while building the nanostar reproduction (point 2). At a much smaller z (e.g. 1e-4), the same sigma and k_intra recover the textbook superselective shape exactly: α starts near 1 at low n, rises to a real peak (~6.4 in one tested case) at intermediate n, then decays toward 0 as capture probability saturates. This is now documented directly in the function's own docstring and is a real consideration for Aim 3/4 calibration: whichever way z ends up calibrated will determine whether the superselectivity mechanism is actually operating in the final model, or present as unused machinery.
5. **Seven new unit tests added (15 total, all passing)**, covering both the `selectivity_alpha()` dilute-regime finding and the nanostar reproduction (including its two honest limitations, which are now checked facts rather than one-off script output that could silently stop being true). All 8 original tests, both original demo scripts, and the patch-regime-sensitivity script were re-run this session and reproduce exactly (1.69 SE, 2.2×10⁻¹⁶, 28.1%, and the Lennon flip/no-flip numbers all unchanged), confirming the new code didn't silently change anything it wasn't supposed to.

---

## Current state: done vs. not done

### Done

- **Project shape settled** after three failed versions (see "What NOT to revisit").
- **Novelty verified.** Direct search for "superselectivity" + "extracellular vesicle" returns nothing. The application is open.
- **Aim 1 complete**, including the n-vs-areal-density resolution (see v4 for full detail — unchanged this session). Code lives in `ev_capture_model/` (`src/binding_model.py`, `src/monte_carlo.py`, `src/population_model.py`, `tests/test_aim1.py`, plus demo scripts and figures).
- **Verified by independent execution**, not just by reading the docs:
  - Analytical vs. Monte Carlo agree within **1.69 standard errors** over 200,000 trials.
  - Population demo produces a **28.1% relative yield difference** between Poisson and negative binomial at *identical means* — the central claim, demonstrated numerically.
  - All 15 unit tests pass (8 original + 7 new this session).
- **Aim 2, legs 1 and 2 complete.** Analytical vs. simulation ✓. Limiting cases ✓ — Poisson theorem holds to 2.2×10⁻¹⁶ at the bond-count partition function level; valency=1 recovers Langmuir.
- **NEW — Aim 2, leg 3 substantially advanced (was fully unstarted in v4).** The DNA nanostar optimal-valency-3 signature (Linne et al. 2024) has been reproduced at a mechanism-faithful (not digit-faithful) level, via a new, honestly-scoped, separate model (`src/nanostar_validation.py`) built specifically because the EV model's own mechanics could not represent the relevant physics (finite valency + surface saturation). `binding_model.selectivity_alpha()` has also finally been exercised (never tested before this session), surfacing a real finding about the z-dependence of superselectivity in the existing model. See "What changed since v4" above for full detail. **Still not done: the Mizenko 10x/5x rows of the Validation Targets table** — these remain correctly assigned to Aim 3's K_INTRA/z/k_geom calibration, not to this validation leg, since Mizenko's paper reports fold-underperformance numbers, not a superselectivity α.
- **Poisson theorem refined.** It is *exactly* true at the level of the bond-count partition function `q_b`, and only approximately true (1–5%) one step later after the Langmuir conversion to a probability. This refinement is mine and is more precise than the original claim. The tests deliberately check the specific claim that is exact.
- **Monte Carlo is non-circular.** It simulates the primitive Poisson bond-forming process and checks resulting frequencies against the closed form, rather than drawing against the formula's own output.
- **Figures produced:** `population_shape_effect.png`, `analytical_vs_mc.png`, `patch_regime_sensitivity.png`, and **NEW** `nanostar_valency_validation.png` (bound-fraction curves, selectivity curves, and peak-selectivity-vs-valency bar chart, analogous to Linne et al.'s own Figure 2b/2d).
- **Mentor/BSL-2 no longer blocking.** Fully computational on public data requires no Qualified Scientist and no lab access. An ISEF Adult Sponsor is still required; confirm current forms with the fair's Scientific Review Committee.
- **Lennon Figures 2 & 3 fully digitized.** Per-vesicle point clouds for PANC-1/HPDEC × EGFR/CA19-9 (Figure 2, cell lines) and healthy/PDAC patient plasma × EGFR/CA19-9 (Figure 3c). Saved as CSVs with source-image crops and full digitization scripts. See "Data situation" below. Unchanged this session.
- **Sign-inversion decomposition, run on real digitized data.** Preserved from v3/v4, still current — see "The decomposition" below.
- **`n` vs. areal receptor density resolved** (v3→v4 work, unchanged this session). `binding_model.py` has a `patch_exponent` parameter with two named regimes, a physical justification for the default, and a working analysis run against the real digitized data.
- **Finding: the ranking flip depends on population size-similarity, not just marker choice.** Unchanged this session — re-verified to reproduce exactly.
- **Test suite: 15 tests, all passing** (8 from Aim 1/2 legs 1-2, 7 new this session for Aim 2 leg 3).

### Not done

- **Aim 2, leg 3 — Mizenko's 10x/5x fold-underperformance numbers specifically.** Still not attempted as a validation target in the sense of testing `selectivity_alpha()` or any other function against it — correctly, since these numbers are fold-underperformance ratios, not a superselectivity α, and belong to Aim 3's calibration step (K_INTRA, z, k_geom), not to this validation leg. Worth stating explicitly so a future session doesn't waste time trying to force-fit Mizenko's numbers into an alpha-style test.
- **`K_INTRA`, `z`, and `k_geom` are placeholders**, correctly labeled. Current figures show behavior driven by invented constants. These get calibrated in Aim 3.
- **NEW — whichever `z` value gets used in Aim 3 calibration will determine whether the calibrated model is actually superselective (α>1 achievable) or not** — this is now a concrete, documented consideration for that calibration step, not an abstract worry (see "What changed since v4," point 4, and the updated Immediate Next Steps below).
- **NEW — the nanostar reproduction is mechanism-faithful, not digit-faithful**, and its k=3 optimum is not robust to every possible K_A (documented and tested, not hidden). If the paper's exact fitted Figure 2e values become available (e.g. from the ESI/4TU.ResearchData repository referenced in the paper's Data Availability statement), they should replace the representative parameters used here for a tighter, digit-level validation. Worth a specific check before treating this validation leg as fully closed.
- **Aims 3 and 4 entirely** — no target ranking, no multi-marker joint distribution, no independence test, no sensitivity analysis. (Partial exception: the raw material for Aim 3's PANC-1/HPDEC and healthy/PDAC cases exists as digitized point clouds — see Done above — but no model has been run against them yet.)
- **Corona/accessibility not in the model.**
- **Incubation-time axis not implemented.**
- **Cancer type / cell line not finally decided** — and the basis for the earlier lean toward ovarian changed with the Lennon finding.
- **ISEF forms, Adult Sponsor, AI-use log** not started.
- **The soft-adhesion contact-patch scaling (`p=2` default) is a physically-motivated approximation, not a direct EV measurement.** Unchanged from v4 — no published study was found measuring actual EV-surface contact-patch area as a function of vesicle diameter.
- **`k_geom` (the contact-patch geometric constant) is uncalibrated**, same status as `K_INTRA` and `z`. Unchanged from v4.
- **Figure 3c gate polygon was hand-approximated, not pixel-digitized.** Unchanged from v4 — low priority, optional tightening.

---

## The four aims

### Aim 1 — Build the model *(complete, including the n-vs-areal-density resolution)*
Multivalent binding model for a vesicle carrying `n` receptors contacting a surface with ligand density `σ`. Analytical closed form plus independent Monte Carlo. Population-level version integrating over a receptor distribution P(n) rather than a fixed `n` — this population step is the core novelty and is what the field's fold-change approach skips. Unchanged this session — see v4 for full detail.

### Aim 2 — Validate three ways *(all three legs now have substantive work, leg 3 not fully closed)*
1. Analytical vs. Monte Carlo. ✓
2. Limiting cases: Poisson → mean-field exact; valency 1 → Langmuir. ✓
3. **Reproduce known experimental signatures.** ~ — DNA nanostar valency-3 optimum reproduced at a mechanism-faithful level (new this session, see "What changed since v4"); Mizenko's 10x/5x numbers remain assigned to Aim 3 calibration, not this leg. `selectivity_alpha()` exercised for the first time, surfacing the z-dependence finding.

### Aim 3 — Re-rank real targets, single and multi-marker *(not started, but raw data now exists for the primary case, and a new calibration consideration has emerged)*
Run the model on published receptor distributions; compute achievable yield and selectivity; compare to fold-change ranking. **The deliverable is the disagreement.**

For the PANC-1/HPDEC and healthy/PDAC-plasma cases specifically, real digitized per-vesicle joint (diameter, copy-number) distributions are now available (`ev_capture_data/`) rather than needing reconstruction from summary statistics alone. Every other system in the literature table still needs summary-statistic reconstruction.

**NEW consideration for this Aim, surfaced this session:** whatever `z` value calibration lands on will determine whether the calibrated model is actually in the superselective regime (α>1 achievable) or not (see "What changed since v4," point 4). This should be checked explicitly once K_INTRA/z/k_geom are calibrated against Mizenko's numbers — it would be a notable finding either way (if the calibrated model turns out non-superselective, that's itself worth reporting and reconciling with the project's central "nonlinearity matters" claim; the two are not the same claim, and the difference should be stated clearly in the write-up).

**Multi-marker capture is a core sub-aim, not a footnote.** Unchanged from v4 — see below for full detail preserved.
- Insurance if single-marker ceilings are low (likely — expect 2–4×, not 20×).
- The field already recommends combinations without quantifying them (Mizenko et al. say "use a combination," never which or at what cost).
- **The independence assumption is itself a finding.** The 9× math assumes the proteins appear independently across vesicles. If they cluster, the multiplication fails. Lennon et al. observed non-uniform EGFR distribution on EV membranes — regions with no detectable EGFR — which *is* clustering. That makes the independence question empirical, not hypothetical.

### Aim 4 — Sensitivity analysis: tell the field what to measure *(not started, but two axes now have preliminary answers)*
Which unknown parameter dominates uncertainty in the selectivity prediction? If distribution *shape* dominates the *mean* by a large factor, the recommendation is concrete: stop reporting mean copy number, report the distribution. Changing what a field measures is a bigger contribution than any single device. The accessibility/corona factor is a natural additional axis here.

**The raw-count-vs-areal-density axis** (from v4, unchanged): the ranking flips for the patient-plasma comparison (large size difference between populations) and does not flip for the cell-line comparison (matched sizes) — a strong, reportable Aim 4 finding.

**NEW — a second preliminary sensitivity axis from this session: whether the model is even in the superselective regime depends sensitively on the activity parameter `z`.** This is conceptually a sensitivity finding of the same flavor as the n-vs-density one: a modeling/calibration choice that isn't obviously important turns out to determine a qualitative property of the result (here: whether α ever exceeds 1 at all). Worth treating as a first-class Aim 4 axis alongside distribution shape and n-vs-density, once z is calibrated in Aim 3.

---

## Three physical mechanisms the model must include

Unchanged this session. Preserved in full from v4:

A naive calculation says capture shouldn't work: a vesicle is a sphere touching a flat surface, so only ~3% of its surface is in the contact patch, giving well under one receptor in contact. But capture demonstrably works. Three mechanisms resolve this, and **all three are real and they interact**:

**(a) Lateral diffusion / receptor recruitment.** EV membranes are fluid, so receptors drift into the contact zone after the vesicle settles. Jahnke et al. 2024 (*J Extracell Vesicles*, DOI 10.1002/jev2.12436) calculate this takes about a second or less for a 100 nm vesicle. Real capture runs 15 minutes to overnight — 10³ to 10⁵ times longer than needed. The static 3% snapshot is the wrong quantity. Independent confirmation from virology (Brownian Adhesive Dynamics, PMC1305223): binding is much greater when receptors can diffuse into the contact zone than when held fixed. **Models that hold receptors still systematically underpredict capture.**

**(b) Vesicle deformation.** Vesicles are soft and flatten on contact. Jahnke et al.: a vesicle engaging 5–50 bonds at ~4 k_BT each, against a bending rigidity of 20–30 k_BT, means a few bonds already deform it. **This is a feedback loop** — bonds form → vesicle flattens → patch grows → more receptors fit → more bonds. **This mechanism is the explicit basis of the resolved n-vs-areal-density question** (see Aim 1): the same low bending rigidity that makes this feedback loop possible is also the physical reason the model's default treats contact-patch area as scaling with diameter² (soft-adhesion/"liquid" contact-mechanics limit for small, floppy bodies).

**(c) Ligand grip strength.** Antibodies bind roughly 1000× more tightly than RGD peptide (nanomolar vs. micromolar Kd, measured head-to-head, *Anal. Chem.* DOI 10.1021/acs.analchem.7b00554). Working capture devices all use antibodies. Reference values: linear RGD binds αvβ3 at ~89 nM; heptapeptide GRGDSPK at ~12.2 nM (*Biophys Rev* DOI 10.1007/s12551-023-01084-3).

**(d) Accessibility / protein corona.** EVs in biofluids acquire an adsorbed plasma-protein layer that occludes surface epitopes. Implement as a single multiplier `f_accessible` on `n`, not as a mechanistic corona model. If `f_accessible` differs between tumor and normal EVs — plausible, since corona composition depends on surface charge and lipid composition — it contaminates the selectivity ratio, which was supposed to be robust because errors cancel. Carry it as an uncertainty and test sensitivity. Source: Försönits et al. 2025, *J Extracell Vesicles*, DOI 10.1002/jev2.70124.

**A testable prediction falls out of (a):** incubation time is a selectivity knob nobody is turning. Short incubation → no time to reshuffle → only richly decorated vesicles stick → selective, low yield. Long incubation → sparse vesicles gather enough bonds → high yield, worse selectivity. Mapping this trade-off appears to be novel.

---

## My prior result (preserve this — it's the logical spine)

I proved that **for Poisson-distributed receptors, mean-field theory is exact** (exactly, at the `q_b` partition-function level; approximately, within 1–5%, after the Langmuir conversion). Restated: *if* receptor counts were Poisson, the mean would be sufficient and fold-change ranking would be correct.

Real EV distributions are measurably **not** Poisson — zero-inflated and overdispersed (over half of vesicles carry none of the standard markers).

Argument structure:
1. Under Poisson, fold-change ranking is provably correct — **my result, the null model**
2. Real distributions are non-Poisson — **literature**
3. Therefore fold-change ranking has a quantifiable error — **my contribution**
4. Here is the corrected method and how much the ranking changes — **my result**

Step 1 doubles as a unit test.

---

## Validation targets (specific numbers the model must reproduce)

| Target | Value | Source | Status |
|---|---|---|---|
| Analytical vs. Monte Carlo | agreement within noise | internal | ✓ 1.69 SE |
| Poisson limiting case | mean-field exactly recovered | my derivation | ✓ 2.2e-16 |
| Valency = 1 | Langmuir binding recovered | standard | ✓ |
| DNA nanostar optimal valency | 3 | Linne et al., *Soft Matter* DOI 10.1039/d4sm00704b | **~ mechanism-faithful reproduction (NEW this session) — see "What changed since v4"; not a digit-level match to the paper's own fitted constants, which were not recoverable as text this session** |
| Single-tetraspanin capture underperformance, EpCAM | ~10× fewer vs. non-specific | Mizenko et al. 2021 | ✗ (belongs to Aim 3 calibration, not this validation leg — see Not Done above) |
| Same, Her2 | ~5× fewer | Mizenko et al. 2021 | ✗ (same) |

The last two are Aim 3's calibration targets: `K_INTRA`, `z`, and `k_geom` get calibrated against the Mizenko pair — three free parameters (as of v4), two measured targets. If physically sensible values can't fit both simultaneously, that's informative: a mechanism is missing.

---

## Key literature

**The anchor paper — Mizenko et al. 2021**, *J Nanobiotechnology* 19:250, DOI 10.1186/s12951-021-00987-1. Open access, ~171 citations. Unchanged from v4 — see full detail preserved below.

They measured experimentally that capturing EVs by a single surface protein biases which tumor markers you find — roughly tenfold fewer EpCAM-positive and fivefold fewer Her2-positive vesicles via single-tetraspanin capture versus non-specific capture. Their interpretation is my thesis, stated by them: a large fraction of non-tetraspanin-expressing vesicles may carry the bulk of these markers. Their recommendation is my multi-marker idea: study multiplexing efficiency before choosing a single capture protein.

**Cite prominently as motivation.** What they did NOT do, and what remains mine:
- No theory — they report 10×, can't explain *why* 10× rather than 3× or 50×, can't predict for untested proteins
- No prediction — every new protein/cell line/surface needs a new experiment
- No quantified trade-off — "use a combination" without saying which or at what yield cost
- Only three proteins on one cell line

Framing: *the bias is known and measured — here is the physical theory that explains its magnitude, predicts it for untested proteins, and quantifies the yield cost of the fix everyone recommends but nobody has sized.* **Explaining a known phenomenon is safer than discovering a new one** — it hands you real numbers to validate against.

**NEW — the nanostar validation paper: Linne, Heemskerk, Zwanikken, Kraft & Laan 2024**, *Soft Matter* 20:8515-8523, DOI 10.1039/d4sm00704b, PMC11484159. Open access.

DNA-origami nanostars with a fixed number of branched arms (valency `k`), each with a single-stranded "sticky end" that binds a complementary receptor on a supported lipid bilayer, imaged by TIRFM. Measured bound fraction Y as a function of receptor density σ_R, for k=1,3,4,10 (main valency sweep) and for three different sticky-end binding strengths at fixed k=6. Key findings, all directly relevant to this project's theoretical foundation:
- Selectivity α = d ln(Y)/d ln(σ_R) peaks at an *intermediate* valency, k=3 — not the largest valency tested (k=10). This is the headline external signature this project's Aim 2 leg 3 aimed to reproduce.
- The mechanism: fitting their model separately to each valency shows the first-bond equilibrium constant K_A stays roughly the same across k, while the subsequent-bond constant K_intra *decreases* as k increases — larger nanostars are less likely to successfully engage additional arms, plausibly due to steric/geometric constraints between arms. This trade-off (more arms = more combinatorial boost, but each individual arm engages less well) is what produces the interior optimum.
- Binding strength also has an optimum (their k=6 sweep): too strong and the transition happens at too-low receptor density to resolve; too weak and multivalent effects wash out.
- Their base "independent arms" model does not, by itself with one shared K_A/K_intra pair, predict the k=3 optimum — reproducing it required their extended model with a ligand-ligand cooperativity/competition free-energy term (Δ_G_int), fit via kinetic Monte Carlo. This project's reproduction (`src/nanostar_validation.py`) does NOT implement that extended kinetic model — it instead uses their base independent-arms equations (their Eq. 2-4) with representative, per-k K_intra values chosen to match their *reported qualitative trend* (K_A flat, K_intra decreasing with k), which is sufficient to reproduce the *qualitative* k=3 optimum but not a claim of matching their full mechanistic explanation. This is a real scope limit, stated honestly — see "What changed since v4" and `src/nanostar_validation.py`'s module docstring for the full account.
- Their exact fitted numeric values (Figure 2e's table, giving K_A and K_intra for k=1,3,6,10) were not recoverable as machine-readable text from either the PMC full-text API or the RSC-hosted PDF this session — both were tried directly. If these become available later (e.g. via the paper's ESI or its 4TU.ResearchData deposit, referenced in the paper's Data Availability statement), they should replace the representative parameters currently used for a tighter, digit-level validation.

**The PDAC data source: Lennon et al. 2019**, *J Extracell Vesicles* 8:1685634, DOI 10.1080/20013078.2019.1685634, PMC6844376. Open access. Unchanged from v4 — see full detail preserved below.

Quantitative single-molecule localization microscopy (qSMLM) counting individual molecules on individual vesicles. What it contains:
- **PANC-1 (PDAC) vs. HPDEC (normal pancreatic ductal epithelial)** — a matched tumor/normal pair, structurally analogous to MDA-MB-231/MCF-10A
- Per-vesicle scatter plots: each point is one EV, y-axis molecules detected, x-axis apparent diameter. **This is distributional data in digitizable form — and it is now digitized.** See "Data situation" below.
- PANC-1 EVs averaged 28 detected EGFR per EV at 51 nm apparent diameter (via cetuximab-AF647); average diameter 72 nm by membrane lectin stain, 75 nm by TEM
- HPDEC yielded fewer EVs with lower EGFR content; no significant diameter difference
- PANC-1 yielded significantly more EVs with higher CA19-9 content than HPDEC
- Patient plasma (5 PDAC, 6 healthy): 5-fold and 15-fold increases in EGFR- and CA19-9-enriched EV numbers overall; **35-fold and 130-fold in the negative-gated PDAC-enriched region**
- **The PDAC-enriched region contained larger EVs — see "The decomposition" below for the corrected, digitized-data version of what happens to target-molecule count in that region.**
- Observed non-uniform EGFR distribution on membranes — regions with no detectable EGFR. That is clustering, relevant to the independence question.
- EVs isolated by UF-SEC (Izon qEV), EV-depleted FBS used, MISEV-consistent marker verification (CD63, TSG101)
- Proteomics deposited: PRIDE PXD012660. RNA: exRNA Atlas EXR-KJENS1PANCANC1-AN.

**Important caveat the authors state:** they report **detected** densities, not absolute ones, because labeling and detection efficiency are imperfect. Their counts are systematic underestimates by an unknown factor. Partly cancels in a tumor/normal ratio — another argument for ratios over absolutes. This belongs in the limitations section.

**The decomposition, run on the actual digitized per-vesicle point cloud from Figure 3c (unchanged since v3/v4):**

| Axis | What the digitized data shows | Direction |
|---|---|---|
| Fraction/count of marker-positive vesicles | Whole plasma population: 5× (EGFR) / 15× (CA19-9) more positive EVs in PDAC per the paper's text. Gated PDAC-enriched region specifically: 35× (EGFR) / 130× (CA19-9) per the paper's text. | **Strongly toward tumor** |
| Raw copies per vesicle, gated population vs. all healthy | Digitized: gated PDAC-enriched mean 49.8 copies/EV (EGFR) and 58.5 copies/EV (CA19-9), vs. all-healthy mean 42.7 and 30.0 respectively. **This is 1.2–2× MORE, not fewer.** | **Toward tumor** |
| Surface density (copies per nm² of membrane), gated population vs. all healthy | Digitized: gated PDAC-enriched surface density is ~0.079× (EGFR) and ~0.119× (CA19-9) of the healthy mean — an 8- to 13-fold **decrease**. Driven by gated-population mean diameter (176–192 nm) being roughly 3.5× the healthy mean diameter (48–55 nm). | **Away from tumor** |

**What this means:** the paper's own phrase "fewer target molecules" is a **density** statement, not a **raw-count** statement — the digitized data shows the gated PDAC-enriched population actually carries as many or slightly more absolute target molecules per vesicle than a typical healthy EV, but spread across so much more membrane (larger diameter) that the areal density is far lower. This directly motivated the n-vs-areal-density model resolution in v4 — unchanged this session.

Cross-checked using only exactly-resolved (pixel-exact, non-density-estimated) digitized points, independent of the density-reconstruction methodology used for the saturated ink regions — gives the same qualitative and closely the same quantitative answer (surface-density ratio 0.070–0.080× vs. 0.079–0.119× for the full reconstruction). Not an artifact of the digitization method.

**What full pixel-level access would still add:** the gate polygon itself was hand-approximated rather than pixel-digitized — a more precise digitization of the dashed gate line would tighten the exact fold-numbers, though the qualitative/order-of-magnitude conclusion is not expected to change given the cross-check above.

**Theory foundation (cite, don't reinvent):**
- Martinez-Veracoechea & Frenkel, *PNAS* 2011;108(27):10963 — founding superselectivity paper. Cites αvβ3 targeting as motivation, never connects to EV capture.
- Dubacheva, Curk, Frenkel & Richter, *JACS* 2019;141:2577 — defines α = d ln θ / d ln n_R; α > 1 is the superselectivity signature.
- Curk, Dobnikar & Frenkel — the physics the current code was built from.
- Linne, Heemskerk, Zwanikken, Kraft & Laan, *Soft Matter* 2024, DOI 10.1039/d4sm00704b — nanostar optimal valency; **now reproduced at a mechanism-faithful level, see above and "What changed since v4."**
- *Accounts Chem Res* DOI 10.1021/acs.accounts.2c00672 — superselectivity as a design principle for biomarker surfaces; generalizes beyond nanostars; EVs conspicuously absent.

**Other single-vesicle data sources:** (unchanged from v4)
- Liu et al. 2026, *Advanced Science*, DOI 10.1002/advs.202515129 — only single-EV integrin copy-number study. 12-plex DNA-PAINT at single-copy resolution. Integrins heterogeneously expressed at low copy number; most exosomes carry mis-paired or unpaired integrins. Measured *combinations*, so directly relevant to multi-marker. **Paywalled and bot-blocked** — needs library access.
- Tetraspanin copy numbers by super-resolution on seminal exosomes: CD63 ~12.6, CD81 ~1.6, CD9 ~16.6 average copies. Reviewed in *Chem Biomed Imaging* DOI 10.1021/cbmi.3c00095 — **use this peer-reviewed review rather than the US Patent 12,313,519 cited in v1.**
- Same review reports the Lennon EGFR/CA19-9 numbers (~24 and ~31 in healthy plasma EVs, dropping to ~16 in PDAC) — useful cross-check.
- Deconvolution work: **53.5% of membrane-stained vesicles showed no significant expression of CD63, CD81, or CD9.** The zero-inflation fact.
- *Cytotherapy* 2025, DOI S1465-3249(25)00128-8 — calibrated nano-flow cytometry, copies per particle across three stem-cell sources: 50–70% CD81+, 40–60% CD63+, 15–35% CD73+; ~30–40 CD81 and ~41–53 CD63 copies per positive particle. Detection limit 4–10 molecules/particle.
- PMC9369185 — ExoView across HeLa, HepG2, H1650, SH-SY5Y.
- Mizenko et al. 2021 — SK-OV-3 (ovarian) with co-occurrence data.

**Contested — do not cite as established, and do not rank as a target:** GPC1 as a pancreatic EV biomarker. Unchanged from v4 — see below for full detail preserved.

Melo et al. 2015 *Nature* reported near-perfect discrimination; roughly ten subsequent studies produced contradictory results. Specifically checked: one study found GPC1 enriched in PDAC circulating exosomes relative to source tumors, but could not distinguish PDAC from benign pancreatic disease using exosomal GPC1, and found no significant GPC1 difference between normal and PDAC tissue even in matched pairs (Oncotarget, PMC5922373). Another, using direct mass-spec quantification, found exosomal GPC1 "not diagnostic for PDAC" with levels relatively similar across normal, PDAC, and chronic pancreatitis (PMC5386003).

**No single-vesicle copy-number distribution exists for GPC1** — searched specifically, found nothing. Every GPC1 measurement in the literature is bulk (ELISA, mass spec, flow cytometry panels), never per-vesicle molecule counts of the kind Lennon et al. did for EGFR/CA19-9. This means GPC1 **cannot be entered into Aim 3's target ranking** — there is no distribution to rank it with, only a bulk-level controversy.

A "benign disease causes a moderate density bump, cancer causes an extreme one" mechanism was proposed externally to explain the controversy, with invented copy-number categories (healthy 0, benign 1–4, cancer 15–50+) presented as if measured. **This is not supported by the data.** Do not use this framing or these numbers anywhere in the write-up. GPC1's only legitimate role here is as a **discussion-section case study**: a high-profile bulk biomarker controversy that single-vesicle distributional data could in principle help resolve, but hasn't, because that data doesn't exist.

---

## Data situation (honest)

Unchanged this session. Preserved in full from v4:

**Ruled out:** the Mizenko Zenodo repository (DOI 10.5281/zenodo.5189322) is public and CC0 but contains 96.5 GB of raw `.pgm` microscope images, not extracted per-vesicle numbers. Figure 7.zip (7.2 kB) was downloaded and contains nothing usable. **Do not send me back to this dataset.**

**Lennon Figures 2 & 3 digitized and saved.** Location: `ev_capture_data/` in the workspace folder (persists across sessions; the source PDF pages do not need to be re-uploaded). Contents:
- Per-vesicle CSVs for Figure 2 panels a, d, e, f (PANC-1/HPDEC × EGFR/CA19-9) — exact vesicle counts matching the paper (1651, 197, 736, 104), with a rescaled version whose mean/CV exactly match the paper's reported values.
- Per-vesicle CSVs for Figure 3c (patient plasma, healthy vs. PDAC, pooled across all subjects, both markers) — shape-accurate reconstruction, no summary-statistic rescaling available since the paper doesn't publish per-point stats for this pooled panel.
- Source image crops for every digitized panel, plus full 400-DPI page renders of both PDF pages, so any panel can be re-examined or re-digitized differently without re-uploading anything.
- Full Python digitization pipeline (axis calibration, blob detection, density-surface reconstruction for the unavoidably-overlapping dense corner of each plot, rescaling, and the decomposition analysis) — reproducible and auditable.
- `README.md` in that folder documents the method, its limits, validation checks, and the full decomposition result in detail. **Read that file before using the CSVs for anything quantitative.**

**Primary method for every OTHER system in the literature table (no per-vesicle figure available):** digitize published per-vesicle scatter plots and histograms where they exist (WebPlotDigitizer, free, browser-based, for any future figure access), and reconstruct distribution families from summary statistics (percent positive + mean copy number + co-occurrence fractions) where only summaries are published.

**Why reconstruction is acceptable:** the claim was never "EV receptor counts follow distribution X." The claim is "capture depends on distribution shape, here is how strongly, and here is the range of answers across all plausible shapes." Testing sensitivity to the shape assumption *is* Aim 4.

**Prepared answer to the circularity objection.** A judge may ask: *"you assumed the shape, isn't that circular?"* No — assuming a family and then reporting sensitivity to that assumption is more honest than fitting one shape and presenting it as truth. This applies less forcefully now to the PANC-1/HPDEC and healthy/PDAC systems specifically, since those now have digitized real distributions rather than reconstructed ones — but it still applies to every other system in the comparison table.

**Physical sanity bound on any fitted distribution:** a 100 nm EV has roughly 31,000 nm² of outer surface; a large receptor ectodomain footprint is tens of nm². If a fitted distribution puts meaningful mass above ~500 copies, it is physically wrong.

**Open lead, not yet pursued:** the Carney lab (UC Davis) published 96 GB under CC0, signalling openness. Rachel Mizenko (first author, ORCID 0000-0002-9991-8799) and Randy Carney (senior, ORCID 0000-0001-8193-1664) could be emailed for tabular per-vesicle data. Also a plausible mentor lead. Deferred by choice, not ruled out. **NEW leads from this session:** Linne et al.'s DNA nanostar data (including, potentially, the exact Figure 2e fit values needed to sharpen the nanostar validation from mechanism-faithful to digit-faithful) is referenced as deposited at 4TU.ResearchData in their Data Availability statement — worth checking directly if a tighter validation is wanted before treating Aim 2 leg 3 as fully closed.

---

## Cancer type — decided by data, not by prior commitment

Unchanged this session. Preserved in full from v4:

**Do not anchor on any one cell line.** The framework is cell-line-agnostic: distributions in, yield and selectivity out. An earlier conversation drifted into anchoring on MDA-MB-231 out of inertia; I caught it. Don't repeat that.

| System | Data quality | Notes |
|---|---|---|
| **PANC-1 vs. HPDEC (pancreatic)** | **Best matched pair — now with digitized per-vesicle data** | Per-vesicle counts, EGFR + CA19-9, tumor AND normal, plus patient plasma. Open access. Figures 2 and 3 fully digitized as point clouds. |
| SK-OV-3 (ovarian) | Strong | Co-occurrence + three tumor markers; no matched normal line |
| HeLa, HepG2, H1650, SH-SY5Y | Good | Standard tetraspanins only |
| Stem cells (3 sources) | Good | Actual calibrated copy numbers |
| MDA-MB-231 (breast) | Thin | Scattered percentages |

**Current plan:** the Lennon data makes pancreatic a data-supported primary case rather than only a forward-looking one. Validate across everything with adequate data and report the comparison.

**Personal note:** my grandfather died of pancreatic cancer. It belongs in the motivation and conclusion, stated plainly once. The project should not structurally depend on it — but the data no longer forces a choice between the story and the evidence.

---

## Scope: what this project does and doesn't claim

Unchanged this session. Preserved in full from v4:

The project shifted from engineering (*build a thing that works*) to science (*determine what's possible and how to choose well*).

- **The model predicts both yield and selectivity** — same calculation. But **selectivity is the more trustworthy output** because it's a ratio: uncertainties in ligand density, contact area, incubation conditions, and detection efficiency affect tumor and normal predictions the same way and partially cancel.
- **The real deliverable is the trade-off curve**, not a single optimum. Few ligands → low yield, good selectivity. Many ligands → high yield, poor selectivity. The curve also shows which combinations are *physically impossible* regardless of engineering.
- **Relevance is asymmetric, state it honestly:** high specificity matters most for *diagnostics* (tumor EVs are a tiny minority in blood; low yield is fixable by running more sample, poor specificity isn't). For *therapeutic* filtering, yield matters much more.
- **Hard ceiling worth remembering:** if only 20% of vesicles carry a target, yield can never exceed 20%. Bulk measurements cannot see this.
- **A transferable idea worth developing:** a *computational* limit of detection. Given measured distributions, what is the minimum tumor-EV fraction that any affinity surface could detect?

**The decomposition that matters most.** Published copy numbers average over all vesicles including zeros, blending two different quantities:

> **copies per carrier = overall average ÷ fraction carrying**

Worked example (positivity fractions estimated — illustrative, not publishable): CD9 17.0 ÷ ~0.33 ≈ 52. CD81 1.6 ÷ ~0.10 ≈ 16. A tenfold apparent gap collapses to threefold.

**This decomposition has been run on the real digitized Lennon matched pair** — see "The decomposition" above for the full result. Unchanged this session.

---

## What NOT to revisit

Settled. Reopening these wastes time.

- ❌ **Hydrogel formulation optimization / DOE / GPR / Bayesian optimization.** Circular — no measured objective function exists. Killed two prior attempts.
- ❌ **Calibrating against published capture efficiencies.** Surveyed; metrics incompatible, ligand density essentially never reported. Dead end.
- ❌ **The Mizenko Zenodo dataset.** 96 GB of images. Checked.
- ❌ **Pure literature analysis with no model.** Explicitly rejected — I need to build something.
- ❌ **Anchoring on a cell line for non-data reasons.**
- ❌ **Modelling receptor counts as Poisson or Gaussian.** Poisson is the exact distribution under which my central finding vanishes by construction. Gaussian is wrong for count data. **Negative binomial** is correct and is what the existing code uses. This was recommended by an external chatbot; it is wrong.
- ❌ **Defining success as "the model captures cancer EVs and ignores healthy ones."** That is a predetermined outcome, not a hypothesis test.
- ❌ **Using GPC1 as a ranked target, or citing any "density category" numbers for it.** No single-vesicle GPC1 data exists.
- ❌ **Re-deriving from scratch whether contact-patch area scales with vesicle diameter.** Settled in v3→v4 — see `binding_model.py`'s module docstring. DO revisit specifically if an EV-specific contact-patch measurement turns up.
- ❌ **Re-litigating whether the Figure 3c sign-inversion finding is "real."** Checked twice independently (v2, v3), further resolved in the model itself in v4. Don't restart this check from scratch.
- ❌ **NEW — Re-litigating whether the nanostar valency-3 optimum can be reproduced by simply relabeling `binding_model.py`'s variables.** This was tried in reasoning (not in code, since the mismatch was caught before writing anything) and confirmed not to work: our EV model has no valency parameter in the sense Linne et al. use it, and no surface-saturation mechanics, which is the actual source of their optimum. `src/nanostar_validation.py` is a deliberately separate model built from their own equations — don't attempt to fold this back into `binding_model.py` or treat it as redundant with it.
- ❌ **NEW — Assuming `binding_model.selectivity_alpha()` at z=1.0 (the value used in every existing demo) shows superselectivity.** It doesn't — checked and documented this session (see "What changed since v4," point 4). Don't cite any existing demo figure as showing α>1 without checking the z value used; use a small z (e.g. 1e-4) specifically if a genuinely superselective example is needed for a figure.
- ❌ **NEW — Trying to force Mizenko's 10x/5x numbers into a `selectivity_alpha()`-style test.** These are fold-underperformance ratios (single-tetraspanin capture vs. non-specific capture), not a superselectivity α measurement. They belong to Aim 3's K_INTRA/z/k_geom calibration step. Confirmed this session while scoping Aim 2 leg 3 — don't re-investigate whether Mizenko reports an α; they don't.
- ⚠️ **"Capture efficiency" across papers is not comparable.** Denominators vary wildly, and ultracentrifugation itself recovers only 5–25%. Always check the denominator.
- ⚠️ **Superselectivity is not a threshold or an on/off switch.** It is a steep but continuous power-law response (α > 1 over an intermediate range, then saturation). Do not describe it as "unbreakable velcro" or a density threshold.
- ⚠️ **Classical superselectivity theory is not directly transferable.** The original theory describes a multivalent *probe* discriminating between *surfaces* of differing receptor density. This project inverts the geometry. **That gap is the contribution** — do not paper over it. **This gap is now concretely demonstrated, not just asserted**, by the fact that the nanostar validation required a genuinely separate model rather than a relabeling of our own (see above).

---

## Immediate next steps, in order

Ordering principle: do the thing that could most change the project's direction first. v4's step 1 (Aim 2 leg 3) has been substantially advanced this session (nanostar reproduction done at a mechanism-faithful level; `selectivity_alpha()` finally exercised). Numbering below restarts to reflect the new state.

**1. Calibrate `K_INTRA`, `z`, and `k_geom` against Mizenko's 10× and 5×.** This is now the highest-priority not-started item. Three free parameters, two measured targets — if physically sensible values can't fit all three simultaneously, that's informative and worth flagging rather than forcing a fit. Consider whether `patch_exponent` itself should be fixed at 2.0 (physically-motivated default) rather than treated as a free parameter — probably yes, since it represents a geometric regime choice, not a fitted bond-strength number. **NEW consideration this session: once `z` is calibrated, explicitly check whether the resulting model is in the superselective regime (α>1 achievable via `selectivity_alpha()`) or not** — this is now a known, documented sensitivity (see "What changed since v4," point 4) and should be checked and reported either way, not assumed.

**2. (Optional, worth a quick check before treating Aim 2 leg 3 as fully closed) Try to recover Linne et al.'s exact Figure 2e fitted values.** Check the paper's ESI (electronic supplementary information) and/or their 4TU.ResearchData deposit (referenced in the paper's Data Availability statement) for the actual K_A/K_intra numbers at k=1,3,6,10. If found, use them directly in `nanostar_validation.py` in place of the current representative parameters, for a genuinely digit-level (not just mechanism-level) validation. Not blocking — the current mechanism-faithful reproduction is a legitimate, honestly-scoped result on its own — but would strengthen this validation leg further.

**3. Add `f_accessible` to the model.** One multiplier on `n`. Don't model corona formation mechanistically; carry the uncertainty and test sensitivity in Aim 4.

**4. Aim 3 — single-marker ranking**, starting with the PANC-1/HPDEC and healthy/PDAC-plasma cases using the digitized real distributions in `ev_capture_data/` and the working `joint_population_yield()` function, then extending across every other system with adequate (reconstructed) data, compared against fold-change ranking. **Run every ranking under both `patch_exponent` regimes** (or at minimum, check whether the two populations being compared differ substantially in mean diameter) rather than silently picking one.

**5. Aim 3b — multi-marker joint distributions and the independence test.** Lennon's observed EGFR clustering makes independence an empirical question.

**6. Aim 4 — sensitivity analysis**, including distribution-family choice, `f_accessible`, the raw-count-vs-areal-density sensitivity (re-run with calibrated parameters and across the full literature table), and **NEW: the z/superselectivity-regime sensitivity** discovered this session (does the calibrated model actually land in the α>1 regime, and how sensitive is that to reasonable z uncertainty?).

**7. (Optional, low priority) Tighten the Figure 3c gate-polygon digitization.** Unchanged from v4 — currently hand-approximated, qualitative conclusion not expected to change.

**8. (Optional, worth a quick check before Aim 3) Search specifically for any published measurement of EV-surface contact-patch area vs. vesicle diameter.** Unchanged from v4 — none found so far, search wasn't exhaustive.

**9. Write-up, figures, poster.**

Timeline still fits. Aim 2's external-validation leg — expected to be a significant, open-ended research task — turned out to be tractable within one session once the right (separate) model was identified, though it required real back-and-forth (a paused-and-asked moment about scope, and iterative parameter search once building began) rather than a quick lookup. The remaining big-effort items are Aim 3 calibration and ranking.

---

## Known risks

- **Rankings might barely change.** Mitigation: the headline deliverable is the *relationship* between distribution shape and achievable selectivity, plus where the standard approach breaks down — informative in both directions.
- **The raw-count-vs-areal-density distinction complicates the story, with real numbers.** Unchanged from v4 — see full detail preserved. Mitigation in place: the physical argument for the default was worked out independently of which answer it would produce, and the size-mismatch-dependent flip gives the write-up a clean, defensible way to present this as a genuine, structured sensitivity finding.
- **The physical argument for `patch_exponent=2` rests on general soft-matter contact-mechanics, not an EV-specific measurement.** Unchanged from v4.
- **Three free parameters (`K_INTRA`, `z`, `k_geom`) now need calibrating against only two measured targets (Mizenko's 10x/5x).** Unchanged from v4 — this is Immediate Next Steps item 1 now.
- **NEW — the superselectivity mechanism (α>1) is not automatically present in this model just because binding is multivalent; it requires operating in a sufficiently dilute regime (small `z`).** Discovered this session while finally exercising `selectivity_alpha()`. At z=1.0 (used throughout every existing demo), the model's α never exceeds 1 — it only exceeds 1 at much smaller z. This is not a flaw in the underlying theory (it matches the founding papers' own stated dilute-regime requirement, and was confirmed independently via the nanostar reproduction), but it IS a real risk for Aim 3: if calibration against Mizenko's numbers lands on a `z` that keeps the model non-superselective, that would be an important, reportable finding requiring careful framing (the nonlinearity is still real and still explains fold-change-ranking errors, even if the model doesn't cross the formal α>1 threshold) — not something to quietly work around by picking a different z than calibration supports.
- **NEW — the nanostar valency-3 reproduction is mechanism-faithful, not digit-faithful, and its result is parameter-sensitive.** The k=3 optimum was reproduced using representative (not the paper's own exact fitted) K_A/K_intra values, and was shown NOT to be robust to K_A (a ten-fold weaker K_A shifts the winner to k=6). This is reported and tested honestly, and is consistent with the paper's own point that valency, binding strength, and concentration jointly determine the optimum — but a reviewer or judge could reasonably ask "how do you know this isn't just curve-fitting to get the answer you wanted?" The honest answer, and the one to give: the representative parameters were chosen to match the paper's *stated qualitative trend* (independently reported by them, not reverse-engineered from wanting k=3 to win), and getting ANY interior optimum at all required real iteration and was not guaranteed by construction (the first naive parameter guess produced no optimum, and even the working parameterization does not make every multivalent k superselective, which is itself informative and reported rather than hidden).
- **The Lennon sign inversion.** Unchanged from v4 — status fully resolved as of v4, unaffected by this session's work.
- **Do not import unmeasured "density category" numbers for any target, GPC1 or otherwise.** Unchanged from v4.
- **Detected vs. absolute counts.** Unchanged from v4.
- **Digitization introduces its own uncertainty.** Unchanged from v4.
- **Liu et al. paywalled.** Unchanged from v4.
- **Distributions underdetermined from summary statistics** for every system other than the now-digitized PANC-1/HPDEC/plasma one. This is Aim 4, not a bug.
- **"This is just theory."** Answer: validated against published experimental signatures not used in fitting (now including the nanostar reproduction, at a mechanism-faithful level), and it makes falsifiable predictions.
- **Genuine worst case:** rankings don't move AND the trade-off curve has no interesting structure AND multi-marker analysis says nothing non-obvious. Unchanged from v4 — seems unlikely but not impossible.

---

## Practical notes

- Fully computational on public data → **no Qualified Scientist, no BSL-2**. An ISEF Adult Sponsor is still required; confirm current form requirements with the fair's Scientific Review Committee.
- ISEF requires disclosing AI use. Keep a running log: date, tool, what was asked, what was done with the answer. The research plan, abstract, poster, and citations must be in my own words. **NEW — this session (Aim 2 leg 3: nanostar validation + exercising selectivity_alpha) should be logged**: date, that Claude was used to (a) research the Linne et al. 2024 paper's actual equations and findings, (b) identify and articulate a real physics mismatch between that paper's system and this project's own model before writing any code, (c) implement a separate, faithful reproduction of the published equations, (d) iteratively search for parameters that reproduce the reported qualitative finding, honestly documenting when they didn't work and why, and (e) run `selectivity_alpha()` for the first time, discovering and documenting the z-dependence finding. The reasoning behind treating this as "mechanism-faithful, not digit-faithful" validation, and the honest limitations found, should be described in your own words in the methods/limitations section, not just cited to this session.
- Confirm the affiliate fair's project-year window; judged data collection must fall within 12 continuous months.
- **Data provenance for the write-up:** if the digitized Figure 2/3 point clouds are used in any figure or analysis presented at the fair, the methods section should describe the digitization approach and cite Lennon et al. 2019 as the source of the original figures. **NEW: if the nanostar validation figure/analysis is used, the methods section should similarly describe that this is an independent reproduction of Linne et al.'s published equations using representative (not their own exact fitted) parameters, and cite Linne et al. 2024.** Judges may ask how confident the nanostar match is — the honest answer ("the qualitative mechanism and optimum are reproduced from the paper's own equations and reported trends; the exact published fit values were not accessible this session") is a defensible one.
- Both prior sessions' AI-use logging notes (Lennon digitization; n-vs-areal-density resolution) remain valid and unchanged from v4 — see that document if starting the log fresh.
- A practical note that WORKED this time, worth repeating going forward: **attach both `ev_capture_model/` and `ev_capture_data/` from the start of a new chat.** This session did so (following v4's own closing advice) and lost no round-trips to a missing folder, unlike the v3→v4 session.
