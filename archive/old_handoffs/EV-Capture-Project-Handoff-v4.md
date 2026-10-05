# Project Handoff v4: Computational Model of Multivalent EV Capture

**Supersedes:** `EV-CAP~1.MD` / `EV-Capture-Project-Handoff-v3.md` (v3, August 2026). Where the two disagree, this document wins. Key changes are listed in "What changed since v3" below.

**Purpose:** hand a new chat everything it needs to continue this project without re-deriving anything or repeating dead ends. Read "What NOT to revisit" carefully — this project reached its current form after three prior conversations of failed shapes, a fourth (v2→v3) that finished the figure-digitization step v2 had left blocked, and a fifth (v3→v4) that resolved v3's top-priority open question: whether the model's receptor-count variable should represent total count or areal density.

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

## What changed since v3

v3's blocking item — "Immediate next steps" item 1, deciding whether the model's `n` should represent total receptor count or areal receptor density — is now **resolved**, along with a concrete implementation and a real finding about when the choice matters.

This session worked from a handoff document (v3) that referenced two folders, `ev_capture_model/` (the Aim 1 code) and `ev_capture_data/` (the digitized Lennon figures), neither of which was initially attached. `ev_capture_data/` came through on the first attachment attempt; `ev_capture_model/` was missing until specifically requested and reattached — worth knowing for next time, since a fresh chat starting from this document will need **both** folders attached before it can do anything beyond read the text.

Five substantive updates:

1. **The n-vs-areal-density question is resolved with a physical argument, not a guess.** Checked what soft-matter adhesion physics predicts for how a vesicle's contact-patch area (the part of its membrane actually touching a capture surface) scales with vesicle size. For a very soft, small adhesive body — which an EV is, given a lipid bilayer's low bending rigidity (20–30 kT, per Jahnke et al. 2024) — general contact-mechanics results predict contact *radius* scales linearly with the vesicle's own radius (the "soft/liquid" adhesion limit), rather than the sub-linear scaling that applies to stiffer bodies (the classical JKR result). That means contact-patch *area* scales as diameter², i.e. as a **fixed fraction of total membrane area regardless of vesicle size** — which is exactly the geometric assumption that was already implicit in the original model's use of raw receptor count `n`. This is a real, sourced argument, not a first-principles EV-specific measurement (none exists), and that limitation is stated plainly in the code.
2. **The model now supports both regimes through one new physical parameter, not two disconnected code paths.** Added a "contact-patch scaling exponent" `p` to `binding_model.py`: `p=2` (soft-adhesion default) reproduces the original size-blind model exactly — verified to floating-point precision by unit test — and `p=0` (fixed-patch alternative) makes the model depend on areal density instead. This directly implements the handoff's own suggested option (c), "run both and report whether target rankings flip," as a first-class, reusable model feature rather than a one-off comparison.
3. **Running both regimes on the real digitized data gives a genuine, robustness-checked finding: whether the ranking flips depends on how much the compared populations differ in *size*, not just receptor count.** For the patient-plasma data (Figure 3c, healthy vs. PDAC-gated, ~3× diameter difference), the ranking **flips**: the size-blind model says PDAC-enriched vesicles capture 1.17–2.07× *more* than healthy; the density-driven model says 0.18–0.30× *less*. For the cell-line data (Figure 2, PANC-1 vs. HPDEC, diameters nearly matched), the ranking does **not** flip — both regimes agree tumor vesicles capture 1.24–1.51× more. This makes physical sense (dividing by a near-constant D² barely changes a ranking) and was checked for robustness across a 3×3 grid of placeholder `sigma`/`k_intra` values, holding in all 9 cases for the plasma data.
4. **A units subtlety was caught and fixed before it could produce a misleading result.** The first version of this comparison used one shared geometric constant across both regimes and produced near-zero, meaningless yields under the density-driven regime — a units artifact (the constant has different physical dimensions in each regime: dimensionless for `p=2`, units of length² for `p=0`), not a physical finding. Caught by inspecting the intermediate numbers rather than trusting the first output, and fixed by calibrating the constant separately per regime against a shared reference point (the average healthy vesicle). Documented in `patch_regime_sensitivity.py` so this mistake isn't quietly repeated.
5. **Four new unit tests added (8 total, all passing)**, plus a standalone analysis script (`src/patch_regime_sensitivity.py`) and summary figure (`figures/patch_regime_sensitivity.png`) that reproduce the flip/no-flip finding directly from the digitized CSVs — not just a one-off script output that could silently stop being true. All four original unit tests, both original demo scripts, and all originally reported validation numbers (1.69 SE, 2.2×10⁻¹⁶, 28.1%) were re-run this session and reproduce exactly, confirming the new code didn't silently change anything it wasn't supposed to.

---

## Current state: done vs. not done

### Done

- **Project shape settled** after three failed versions (see "What NOT to revisit").
- **Novelty verified.** Direct search for "superselectivity" + "extracellular vesicle" returns nothing. The application is open.
- **Aim 1 complete.** Analytical binding model, Monte Carlo bond-counting simulation, population-level model integrating over P(n). Code lives in `ev_capture_model/` (built in a separate Claude Code session; `src/binding_model.py`, `src/monte_carlo.py`, `src/population_model.py`, `tests/test_aim1.py`, plus demo scripts and figures).
- **Verified by independent execution**, not just by reading the docs:
  - Analytical vs. Monte Carlo agree within **1.69 standard errors** over 200,000 trials.
  - Population demo produces a **28.1% relative yield difference** between Poisson and negative binomial at *identical means* — the central claim, demonstrated numerically.
  - All four unit tests pass.
- **Aim 2, legs 1 and 2 complete.** Analytical vs. simulation ✓. Limiting cases ✓ — Poisson theorem holds to 2.2×10⁻¹⁶ at the bond-count partition function level; valency=1 recovers Langmuir.
- **Poisson theorem refined.** It is *exactly* true at the level of the bond-count partition function `q_b`, and only approximately true (1–5%) one step later after the Langmuir conversion to a probability. This refinement is mine and is more precise than the original claim. The tests deliberately check the specific claim that is exact.
- **Monte Carlo is non-circular.** It simulates the primitive Poisson bond-forming process and checks resulting frequencies against the closed form, rather than drawing against the formula's own output.
- **Figures produced:** `population_shape_effect.png` (near poster-ready), `analytical_vs_mc.png`.
- **Mentor/BSL-2 no longer blocking.** Fully computational on public data requires no Qualified Scientist and no lab access. An ISEF Adult Sponsor is still required; confirm current forms with the fair's Scientific Review Committee.
- **NEW — Lennon Figures 2 & 3 fully digitized.** Per-vesicle point clouds for PANC-1/HPDEC × EGFR/CA19-9 (Figure 2, cell lines) and healthy/PDAC patient plasma × EGFR/CA19-9 (Figure 3c). Saved as CSVs with source-image crops and full digitization scripts. See "Data situation" below.
- **Sign-inversion decomposition re-run on real digitized data, not just text.** Result refined (see the v3→v4 predecessor of this document, and "The decomposition" below) — preserved from v3, still current.
- **NEW — Immediate next steps item 1 (n vs. areal density) resolved.** `binding_model.py` now has a `patch_exponent` parameter with two named regimes (`p=2` size-blind default, mathematically identical to the original model; `p=0` density-driven alternative), a physical justification for the default (soft-adhesion contact-mechanics scaling), and a working analysis (`src/patch_regime_sensitivity.py`) run against the real digitized data. See "What changed since v3" above for the full result.
- **NEW — Finding: the ranking flip depends on population size-similarity, not just marker choice.** Plasma data (large size difference between healthy/PDAC) flips ranking between regimes; cell-line data (matched sizes) does not. Robustness-checked across a parameter grid.
- **NEW — Test suite expanded to 8 tests, all passing** (4 original + 4 new, covering backward compatibility, the density-only reduction, joint-yield mechanics, and the headline Lennon-data flip).

### Not done

- **Aim 2, leg 3** — reproducing *external* experimental signatures (nanostar optimal valency of 3; Mizenko's 10×/5× bias). Until this is done, every check is internal.
- **`K_INTRA` and `z` are placeholders**, correctly labeled. Current figures show behavior driven by invented constants. These get calibrated in Aim 3.
- **`selectivity_alpha()` is written but never exercised** — no test, no demo. Check it against a published α before relying on it.
- **Aims 3 and 4 entirely** — no target ranking, no multi-marker joint distribution, no independence test, no sensitivity analysis. (Partial exception: the raw material for Aim 3's PANC-1/HPDEC and healthy/PDAC cases now exists as digitized point clouds — see Done above — but no model has been run against them yet.)
- **Corona/accessibility not in the model.**
- **Incubation-time axis not implemented.**
- **Cancer type / cell line not finally decided** — and the basis for the earlier lean toward ovarian changed with the Lennon finding.
- **ISEF forms, Adult Sponsor, AI-use log** not started.
- ~~model does not yet distinguish raw receptor count `n` from areal receptor density~~ **RESOLVED this session** — see "What changed since v3" above. `binding_model.py` now has both regimes as a single `patch_exponent` parameter, with a physically-justified default and a documented, robustness-checked finding about when the choice matters.
- **NEW — the soft-adhesion contact-patch scaling (`p=2` default) is a physically-motivated approximation, not a direct EV measurement.** No published study was found measuring actual EV-surface contact-patch area as a function of vesicle diameter. If one turns up before Aim 3 calibration, it should replace this approximation. Worth a specific literature check, not just a general assumption to carry forward.
- **NEW — `k_geom` (the contact-patch geometric constant) is uncalibrated, same status as `K_INTRA` and `z`.** The value used in `patch_regime_sensitivity.py` was chosen only to normalize the two regimes against each other for a fair comparison (see "What changed since v3," item 4) — it is not yet fit to any real capture data. Gets calibrated alongside `K_INTRA`/`z` in Aim 3, likely needing a third measured target since there are now three free parameters instead of two.
- **Figure 3c gate polygon was hand-approximated, not pixel-digitized.** The dashed negative-gate boundary in Figure 3c was read visually as a diagonal trapezoid rather than digitized from the dashed line's actual pixels. It reproduces the paper's "<5% of healthy EVs fall in gate" criterion closely (0–0.4% in the reconstruction) and the headline finding is insensitive to this (checked against exactly-resolved points only), but a more precise gate digitization would still be a worthwhile, low-effort improvement before publication-quality figures are made.

---

## The four aims

### Aim 1 — Build the model *(complete, including the n-vs-areal-density resolution)*
Multivalent binding model for a vesicle carrying `n` receptors contacting a surface with ligand density `σ`. Analytical closed form plus independent Monte Carlo. Population-level version integrating over a receptor distribution P(n) rather than a fixed `n` — this population step is the core novelty and is what the field's fold-change approach skips.

**The n-vs-areal-density question (open in v3) is now resolved.** `binding_model.py` gained a `patch_exponent` parameter representing how contact-patch area scales with vesicle diameter: `p=2` (soft-adhesion limit, the physically-motivated default) makes contact-patch receptor count proportional to total `n` regardless of diameter — mathematically identical to the original model — while `p=0` (fixed-patch alternative) makes it proportional to areal density instead. Both are now real, tested, working code paths (`q_bound_patch()`, `capture_probability_patch()`, `joint_population_yield()`), not a decision deferred to Aim 3. See "What changed since v3" above for the physical justification and the resulting finding (ranking flips for size-mismatched populations, doesn't for size-matched ones).

### Aim 2 — Validate three ways *(two of three done)*
1. Analytical vs. Monte Carlo. ✓
2. Limiting cases: Poisson → mean-field exact; valency 1 → Langmuir. ✓
3. **Reproduce known experimental signatures.** ✗ — see Validation Targets.

### Aim 3 — Re-rank real targets, single and multi-marker *(not started, but raw data now exists for the primary case)*
Run the model on published receptor distributions; compute achievable yield and selectivity; compare to fold-change ranking. **The deliverable is the disagreement.**

For the PANC-1/HPDEC and healthy/PDAC-plasma cases specifically, real digitized per-vesicle joint (diameter, copy-number) distributions are now available (`ev_capture_data/`) rather than needing reconstruction from summary statistics alone — a stronger starting point than v2 anticipated for this system. Every other system in the literature table still needs summary-statistic reconstruction.

**Multi-marker capture is a core sub-aim, not a footnote.** If a surface requires two proteins before a vesicle sticks, selectivity multiplies (3× and 3× ≈ 9×) at a yield cost. Three reasons it matters:
- Insurance if single-marker ceilings are low (likely — expect 2–4×, not 20×).
- The field already recommends combinations without quantifying them (Mizenko et al. say "use a combination," never which or at what cost).
- **The independence assumption is itself a finding.** The 9× math assumes the proteins appear independently across vesicles. If they cluster, the multiplication fails. Lennon et al. observed non-uniform EGFR distribution on EV membranes — regions with no detectable EGFR — which *is* clustering. That makes the independence question empirical, not hypothetical.

### Aim 4 — Sensitivity analysis: tell the field what to measure *(not started, but one axis has a preliminary answer)*
Which unknown parameter dominates uncertainty in the selectivity prediction? If distribution *shape* dominates the *mean* by a large factor, the recommendation is concrete: stop reporting mean copy number, report the distribution. Changing what a field measures is a bigger contribution than any single device. The accessibility/corona factor is a natural additional axis here.

**The raw-count-vs-areal-density axis (flagged as a candidate in v3) now has a preliminary answer, run on real data:** the ranking flips for the patient-plasma comparison (large size difference between populations) and does not flip for the cell-line comparison (matched sizes). This is itself a strong, reportable Aim 4 finding — the sensitivity of the model's conclusion to this modeling choice depends systematically on how much the compared populations differ in vesicle size, not on the marker or the absolute receptor counts. This still needs `k_geom`/`K_INTRA`/`z` calibration (Aim 3) before the *exact* fold-numbers can be trusted, but the qualitative structure of the finding — and that it's a real, non-arbitrary sensitivity, not a coin flip — is already established.

---

## Three physical mechanisms the model must include

A naive calculation says capture shouldn't work: a vesicle is a sphere touching a flat surface, so only ~3% of its surface is in the contact patch, giving well under one receptor in contact. But capture demonstrably works. Three mechanisms resolve this, and **all three are real and they interact**:

**(a) Lateral diffusion / receptor recruitment.** EV membranes are fluid, so receptors drift into the contact zone after the vesicle settles. Jahnke et al. 2024 (*J Extracell Vesicles*, DOI 10.1002/jev2.12436) calculate this takes about a second or less for a 100 nm vesicle. Real capture runs 15 minutes to overnight — 10³ to 10⁵ times longer than needed. The static 3% snapshot is the wrong quantity. Independent confirmation from virology (Brownian Adhesive Dynamics, PMC1305223): binding is much greater when receptors can diffuse into the contact zone than when held fixed. **Models that hold receptors still systematically underpredict capture.**

**(b) Vesicle deformation.** Vesicles are soft and flatten on contact. Jahnke et al.: a vesicle engaging 5–50 bonds at ~4 k_BT each, against a bending rigidity of 20–30 k_BT, means a few bonds already deform it. **This is a feedback loop** — bonds form → vesicle flattens → patch grows → more receptors fit → more bonds. **This mechanism is now the explicit basis of the resolved n-vs-areal-density question** (see Aim 1 and "What changed since v3" above): the same low bending rigidity that makes this feedback loop possible is also the physical reason the model's default now treats contact-patch area as scaling with diameter² (soft-adhesion/"liquid" contact-mechanics limit for small, floppy bodies) — i.e., a fixed fraction of membrane area regardless of vesicle size, which makes contact-patch receptor count track total `n` directly. This isn't a coincidence: it's the same physical softness driving both the feedback-loop argument here and the `patch_exponent=2` default in the code.

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
| DNA nanostar optimal valency | 3 | Linne et al., *Soft Matter* DOI 10.1039/d4sm00704b | ✗ |
| Single-tetraspanin capture underperformance, EpCAM | ~10× fewer vs. non-specific | Mizenko et al. 2021 | ✗ |
| Same, Her2 | ~5× fewer | Mizenko et al. 2021 | ✗ |

The last three are the external validation. `K_INTRA` and `z` get calibrated against the Mizenko pair — two free parameters, two measured targets. If physically sensible values can't fit both simultaneously, that's informative: a mechanism is missing.

---

## Key literature

**The anchor paper — Mizenko et al. 2021**, *J Nanobiotechnology* 19:250, DOI 10.1186/s12951-021-00987-1. Open access, ~171 citations.

They measured experimentally that capturing EVs by a single surface protein biases which tumor markers you find — roughly tenfold fewer EpCAM-positive and fivefold fewer Her2-positive vesicles via single-tetraspanin capture versus non-specific capture. Their interpretation is my thesis, stated by them: a large fraction of non-tetraspanin-expressing vesicles may carry the bulk of these markers. Their recommendation is my multi-marker idea: study multiplexing efficiency before choosing a single capture protein.

**Cite prominently as motivation.** What they did NOT do, and what remains mine:
- No theory — they report 10×, can't explain *why* 10× rather than 3× or 50×, can't predict for untested proteins
- No prediction — every new protein/cell line/surface needs a new experiment
- No quantified trade-off — "use a combination" without saying which or at what yield cost
- Only three proteins on one cell line

Framing: *the bias is known and measured — here is the physical theory that explains its magnitude, predicts it for untested proteins, and quantifies the yield cost of the fix everyone recommends but nobody has sized.* **Explaining a known phenomenon is safer than discovering a new one** — it hands you real numbers to validate against.

**The PDAC data source: Lennon et al. 2019**, *J Extracell Vesicles* 8:1685634, DOI 10.1080/20013078.2019.1685634, PMC6844376. Open access.

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

**The decomposition, now run on the actual digitized per-vesicle point cloud from Figure 3c (supersedes the text-only version in v2):**

| Axis | What the digitized data shows | Direction |
|---|---|---|
| Fraction/count of marker-positive vesicles | Whole plasma population: 5× (EGFR) / 15× (CA19-9) more positive EVs in PDAC per the paper's text. Gated PDAC-enriched region specifically: 35× (EGFR) / 130× (CA19-9) per the paper's text. | **Strongly toward tumor** |
| Raw copies per vesicle, gated population vs. all healthy | Digitized: gated PDAC-enriched mean 49.8 copies/EV (EGFR) and 58.5 copies/EV (CA19-9), vs. all-healthy mean 42.7 and 30.0 respectively. **This is 1.2–2× MORE, not fewer.** | **Toward tumor** (reversed from v2's reading) |
| Surface density (copies per nm² of membrane), gated population vs. all healthy | Digitized: gated PDAC-enriched surface density is ~0.079× (EGFR) and ~0.119× (CA19-9) of the healthy mean — an 8- to 13-fold **decrease**. Driven by gated-population mean diameter (176–192 nm) being roughly 3.5× the healthy mean diameter (48–55 nm). | **Away from tumor** |

**What this means, corrected from v2:** the paper's own phrase "fewer target molecules" is a **density** statement, not a **raw-count** statement — the digitized data shows the gated PDAC-enriched population actually carries as many or slightly more absolute target molecules per vesicle than a typical healthy EV, but spread across so much more membrane (larger diameter) that the areal density is far lower. The "opposite axes" structure of the original v2 finding survives — abundance (fraction-positive) and something-about-density both matter and move in different directions relative to tumor status — but which "copy number" quantity moves which way depends on whether you count total receptors or receptors per unit area. **This was the single most important open question for Aim 3's model decision — it is now resolved**, see Aim 1 above and "What changed since v3."

Cross-checked using only exactly-resolved (pixel-exact, non-density-estimated) digitized points, independent of the density-reconstruction methodology used for the saturated ink regions — gives the same qualitative and closely the same quantitative answer (surface-density ratio 0.070–0.080× vs. 0.079–0.119× for the full reconstruction). Not an artifact of the digitization method.

**What full pixel-level access would still add:** the gate polygon itself was hand-approximated rather than pixel-digitized (see "Not done" above) — a more precise digitization of the dashed gate line would tighten the exact fold-numbers, though the qualitative/order-of-magnitude conclusion is not expected to change given the cross-check above.

**Theory foundation (cite, don't reinvent):**
- Martinez-Veracoechea & Frenkel, *PNAS* 2011;108(27):10963 — founding superselectivity paper. Cites αvβ3 targeting as motivation, never connects to EV capture.
- Dubacheva, Curk, Frenkel & Richter, *JACS* 2019;141:2577 — defines α = d ln θ / d ln n_R; α > 1 is the superselectivity signature.
- Curk, Dobnikar & Frenkel — the physics the current code was built from.
- Linne et al., *Soft Matter* 2024, DOI 10.1039/d4sm00704b — nanostar optimal valency.
- *Accounts Chem Res* DOI 10.1021/acs.accounts.2c00672 — superselectivity as a design principle for biomarker surfaces; generalizes beyond nanostars; EVs conspicuously absent.

**Other single-vesicle data sources:**
- Liu et al. 2026, *Advanced Science*, DOI 10.1002/advs.202515129 — only single-EV integrin copy-number study. 12-plex DNA-PAINT at single-copy resolution. Integrins heterogeneously expressed at low copy number; most exosomes carry mis-paired or unpaired integrins. Measured *combinations*, so directly relevant to multi-marker. **Paywalled and bot-blocked** — needs library access.
- Tetraspanin copy numbers by super-resolution on seminal exosomes: CD63 ~12.6, CD81 ~1.6, CD9 ~16.6 average copies. Reviewed in *Chem Biomed Imaging* DOI 10.1021/cbmi.3c00095 — **use this peer-reviewed review rather than the US Patent 12,313,519 cited in v1.**
- Same review reports the Lennon EGFR/CA19-9 numbers (~24 and ~31 in healthy plasma EVs, dropping to ~16 in PDAC) — useful cross-check.
- Deconvolution work: **53.5% of membrane-stained vesicles showed no significant expression of CD63, CD81, or CD9.** The zero-inflation fact.
- *Cytotherapy* 2025, DOI S1465-3249(25)00128-8 — calibrated nano-flow cytometry, copies per particle across three stem-cell sources: 50–70% CD81+, 40–60% CD63+, 15–35% CD73+; ~30–40 CD81 and ~41–53 CD63 copies per positive particle. Detection limit 4–10 molecules/particle.
- PMC9369185 — ExoView across HeLa, HepG2, H1650, SH-SY5Y.
- Mizenko et al. 2021 — SK-OV-3 (ovarian) with co-occurrence data.

**Contested — do not cite as established, and do not rank as a target:** GPC1 as a pancreatic EV biomarker. Melo et al. 2015 *Nature* reported near-perfect discrimination; roughly ten subsequent studies produced contradictory results. Specifically checked: one study found GPC1 enriched in PDAC circulating exosomes relative to source tumors, but could not distinguish PDAC from benign pancreatic disease using exosomal GPC1, and found no significant GPC1 difference between normal and PDAC tissue even in matched pairs (Oncotarget, PMC5922373). Another, using direct mass-spec quantification, found exosomal GPC1 "not diagnostic for PDAC" with levels relatively similar across normal, PDAC, and chronic pancreatitis (PMC5386003).

**No single-vesicle copy-number distribution exists for GPC1** — searched specifically, found nothing. Every GPC1 measurement in the literature is bulk (ELISA, mass spec, flow cytometry panels), never per-vesicle molecule counts of the kind Lennon et al. did for EGFR/CA19-9. This means GPC1 **cannot be entered into Aim 3's target ranking** — there is no distribution to rank it with, only a bulk-level controversy.

A "benign disease causes a moderate density bump, cancer causes an extreme one" mechanism was proposed externally to explain the controversy, with invented copy-number categories (healthy 0, benign 1–4, cancer 15–50+) presented as if measured. **This is not supported by the data** — the replication studies above found GPC1 essentially flat across conditions in several cohorts, not an intermediate benign value. Do not use this framing or these numbers anywhere in the write-up. GPC1's only legitimate role here is as a **discussion-section case study**: a high-profile bulk biomarker controversy that single-vesicle distributional data could in principle help resolve, but hasn't, because that data doesn't exist — a clean illustration of the field-wide gap this project is about.

---

## Data situation (honest)

**Ruled out:** the Mizenko Zenodo repository (DOI 10.5281/zenodo.5189322) is public and CC0 but contains 96.5 GB of raw `.pgm` microscope images, not extracted per-vesicle numbers. Figure 7.zip (7.2 kB) was downloaded and contains nothing usable. **Do not send me back to this dataset.**

**NEW — Lennon Figures 2 & 3 digitized and saved.** Location: `ev_capture_data/` in the workspace folder (persists across sessions; the source PDF pages do not need to be re-uploaded). Contents:
- Per-vesicle CSVs for Figure 2 panels a, d, e, f (PANC-1/HPDEC × EGFR/CA19-9) — exact vesicle counts matching the paper (1651, 197, 736, 104), with a rescaled version whose mean/CV exactly match the paper's reported values.
- Per-vesicle CSVs for Figure 3c (patient plasma, healthy vs. PDAC, pooled across all subjects, both markers) — shape-accurate reconstruction, no summary-statistic rescaling available since the paper doesn't publish per-point stats for this pooled panel.
- Source image crops for every digitized panel, plus full 400-DPI page renders of both PDF pages, so any panel can be re-examined or re-digitized differently without re-uploading anything.
- Full Python digitization pipeline (axis calibration, blob detection, density-surface reconstruction for the unavoidably-overlapping dense corner of each plot, rescaling, and the decomposition analysis) — reproducible and auditable.
- `README.md` in that folder documents the method, its limits, validation checks, and the full decomposition result in detail. **Read that file before using the CSVs for anything quantitative** — in particular, it explains which points are pixel-exact (`source = resolved`) versus density-reconstructed (`source = density_estimated`), and why the dense-corner problem is a fundamental limit of digitizing a static scatter image, not a fixable bug.

**Primary method for every OTHER system in the literature table (no per-vesicle figure available):** digitize published per-vesicle scatter plots and histograms where they exist (WebPlotDigitizer, free, browser-based, for any future figure access), and reconstruct distribution families from summary statistics (percent positive + mean copy number + co-occurrence fractions) where only summaries are published.

**Why reconstruction is acceptable:** the claim was never "EV receptor counts follow distribution X." The claim is "capture depends on distribution shape, here is how strongly, and here is the range of answers across all plausible shapes." Testing sensitivity to the shape assumption *is* Aim 4. This makes the shape-sensitivity analysis mandatory rather than optional.

**Prepared answer to the circularity objection.** A judge may ask: *"you assumed the shape, isn't that circular?"* No — assuming a family and then reporting sensitivity to that assumption is more honest than fitting one shape and presenting it as truth. I need this in my own words before February. (This objection applies less forcefully now to the PANC-1/HPDEC and healthy/PDAC systems specifically, since those now have digitized real distributions rather than reconstructed ones — but it still applies to every other system in the comparison table.)

**Physical sanity bound on any fitted distribution:** a 100 nm EV has roughly 31,000 nm² of outer surface; a large receptor ectodomain footprint is tens of nm². If a fitted distribution puts meaningful mass above ~500 copies, it is physically wrong.

**Open lead, not yet pursued:** the Carney lab (UC Davis) published 96 GB under CC0, signalling openness. Rachel Mizenko (first author, ORCID 0000-0002-9991-8799) and Randy Carney (senior, ORCID 0000-0001-8193-1664) could be emailed for tabular per-vesicle data. Also a plausible mentor lead. Deferred by choice, not ruled out. The Lennon/Jenkins group is a second such lead — now somewhat less urgent for the PANC-1/HPDEC/plasma system specifically, since usable per-vesicle data for that system has been extracted from the published figures directly, though authentic (non-digitized) tabular data would still be strictly better if obtainable.

---

## Cancer type — decided by data, not by prior commitment

**Do not anchor on any one cell line.** The framework is cell-line-agnostic: distributions in, yield and selectivity out. An earlier conversation drifted into anchoring on MDA-MB-231 out of inertia; I caught it. Don't repeat that.

Revised data-availability picture:

| System | Data quality | Notes |
|---|---|---|
| **PANC-1 vs. HPDEC (pancreatic)** | **Best matched pair — now with digitized per-vesicle data** | Per-vesicle counts, EGFR + CA19-9, tumor AND normal, plus patient plasma. Open access. Figures 2 and 3 now fully digitized as point clouds (see Data situation). |
| SK-OV-3 (ovarian) | Strong | Co-occurrence + three tumor markers; no matched normal line |
| HeLa, HepG2, H1650, SH-SY5Y | Good | Standard tetraspanins only |
| Stem cells (3 sources) | Good | Actual calibrated copy numbers |
| MDA-MB-231 (breast) | Thin | Scattered percentages |

**Current plan:** the Lennon data makes pancreatic a data-supported primary case rather than only a forward-looking one, and this is now stronger than it was in v2 — real digitized per-vesicle joint distributions exist for this system, not just summary statistics. Validate across everything with adequate data and report the comparison — that is more informative than picking one and removes any suspicion of cherry-picking.

**Personal note:** my grandfather died of pancreatic cancer. It belongs in the motivation and conclusion, stated plainly once. The project should not structurally depend on it — but the data no longer forces a choice between the story and the evidence, and now has an even stronger evidentiary base than it did in v2.

---

## Scope: what this project does and doesn't claim

The project shifted from engineering (*build a thing that works*) to science (*determine what's possible and how to choose well*).

- **The model predicts both yield and selectivity** — same calculation. But **selectivity is the more trustworthy output** because it's a ratio: uncertainties in ligand density, contact area, incubation conditions, and detection efficiency affect tumor and normal predictions the same way and partially cancel.
- **The real deliverable is the trade-off curve**, not a single optimum. Few ligands → low yield, good selectivity. Many ligands → high yield, poor selectivity. The curve also shows which combinations are *physically impossible* regardless of engineering.
- **Relevance is asymmetric, state it honestly:** high specificity matters most for *diagnostics* (tumor EVs are a tiny minority in blood; low yield is fixable by running more sample, poor specificity isn't). For *therapeutic* filtering, yield matters much more — this project is relevant but incomplete there.
- **Hard ceiling worth remembering:** if only 20% of vesicles carry a target, yield can never exceed 20%. Bulk measurements cannot see this — they can't distinguish "every vesicle has 10 copies" from "one in ten has 100."
- **A transferable idea worth developing:** a *computational* limit of detection. Given measured distributions, what is the minimum tumor-EV fraction that any affinity surface could detect? That is a legitimate and strong output, and it is the computational analogue of the LOD an experimental project would report.

**The decomposition that matters most.** Published copy numbers average over all vesicles including zeros, blending two different quantities:

> **copies per carrier = overall average ÷ fraction carrying**

Worked example (positivity fractions estimated — illustrative, not publishable): CD9 17.0 ÷ ~0.33 ≈ 52. CD81 1.6 ÷ ~0.10 ≈ 16. A tenfold apparent gap collapses to threefold.

**This decomposition has now been run on the real digitized Lennon matched pair (superseding the "step 2" item from v2's next-steps list, which is now DONE) — see "The decomposition" above for the full result.** The refined finding: raw copies-per-carrier in the PDAC-enriched population is *higher*, not lower, than in healthy EVs; it's the *areal density* that drops sharply, because those carriers are much larger vesicles.

---

## What NOT to revisit

Settled. Reopening these wastes time.

- ❌ **Hydrogel formulation optimization / DOE / GPR / Bayesian optimization.** Circular — no measured objective function exists. Killed two prior attempts.
- ❌ **Calibrating against published capture efficiencies.** Surveyed; metrics incompatible, ligand density essentially never reported. Dead end.
- ❌ **The Mizenko Zenodo dataset.** 96 GB of images. Checked.
- ❌ **Pure literature analysis with no model.** Explicitly rejected — I need to build something.
- ❌ **Anchoring on a cell line for non-data reasons.**
- ❌ **Modelling receptor counts as Poisson or Gaussian.** Poisson is the exact distribution under which my central finding vanishes by construction (mean-field becomes exact, so shape stops mattering). Using it would produce a clean-looking null result that is an artifact of the assumption. Gaussian is wrong for count data — continuous and unbounded below. **Negative binomial** is correct and is what the existing code uses. This was recommended by an external chatbot; it is wrong.
- ❌ **Defining success as "the model captures cancer EVs and ignores healthy ones."** That is a predetermined outcome, not a hypothesis test. A null or modest result is scientifically valid when framed as hypothesis-testing, and that framing was settled deliberately.
- ❌ **Using GPC1 as a ranked target, or citing any "density category" numbers for it.** No single-vesicle GPC1 data exists. An external chatbot invented copy-number categories to explain the GPC1 replication controversy via a density mechanism; several careful replication studies found GPC1 levels flat across normal/benign/cancer tissue, contradicting that mechanism. GPC1 may be used only as a discussion-section case study of the bulk-vs-single-vesicle gap, never as a modeled/ranked target.
- ❌ **NEW — Re-deriving from scratch whether contact-patch area scales with vesicle diameter.** This session checked general soft-matter/JKR contact-mechanics literature and settled on the soft-adhesion ("liquid limit") scaling as the physically motivated default (`patch_exponent=2`), with the fixed-patch alternative (`patch_exponent=0`) available and tested. Don't restart this literature check from zero in a future session — the reasoning and sources are in `binding_model.py`'s module docstring and in "What changed since v3" above. DO revisit it specifically if an EV-specific contact-patch measurement turns up (see Immediate next steps item 8) — that would be new evidence, not a re-litigation.
- ❌ **Re-litigating whether the Figure 3c sign-inversion finding is "real."** It has now been checked twice independently: once from the paper's prose (v2), once from the actual digitized per-vesicle point cloud including a cross-check restricted to pixel-exact points only (v3). Both agree on the qualitative structure. What changed between v2 and v3 was *which* copy-number quantity moves in which direction — raw count vs. areal density — not whether the phenomenon exists. That question is now further resolved in the model itself (see "What changed since v3" above) — don't restart this check from scratch.
- ⚠️ **"Capture efficiency" across papers is not comparable.** Denominators vary wildly (% of input vs. % vs. ultracentrifugation vs. fluorescence a.u.), and ultracentrifugation itself recovers only 5–25%. Always check the denominator. Belongs in the limitations section.
- ⚠️ **Superselectivity is not a threshold or an on/off switch.** It is a steep but continuous power-law response (α > 1 over an intermediate range, then saturation). Monovalent binding is Langmuir — it rises and saturates, it is not a straight line. The correct figure is two curves on log-log axes with different slopes, not "diagonal vs. hockey stick." Do not describe it as "unbreakable velcro" or a density threshold.
- ⚠️ **Classical superselectivity theory is not directly transferable.** The original theory describes a multivalent *probe* discriminating between *surfaces* of differing receptor density. This project inverts the geometry: a ligand-coated surface discriminating between *vesicles* of differing receptor count, where the number in the contact patch is small and finite. The continuum assumptions are strained. **That gap is the contribution** — do not paper over it.

---

## Immediate next steps, in order

Ordering principle: do the thing that could most change the project's direction first. v3's step 1 (decide n vs. areal density) is now **done** — resolved with a physical argument, implemented, tested, and run against real data (see "What changed since v3" above). Numbering below restarts again to reflect the new state.

**1. Finish Aim 2, leg 3 — NOT STARTED, highest priority now.** Reproduce the nanostar valency-3 optimum. Test `selectivity_alpha()` against a published α value. This converts validation from internal to external, and should happen before further model features are added on top of an internally-validated-only foundation.

**2. Calibrate `K_INTRA`, `z`, and now `k_geom` against Mizenko's 10× and 5×.** Three free parameters now instead of two (the new `patch_exponent` machinery added a geometric constant, `k_geom`, alongside the original `K_INTRA` and `z`) — two measured targets. If physically sensible values can't fit all three simultaneously, or if the fit is underdetermined with only two targets for three parameters, that's informative and worth flagging rather than forcing a fit. Consider whether a third measured target is needed, or whether `patch_exponent` itself should be fixed at a physically-motivated value (2.0, the default) rather than treated as a free calibration parameter — probably the latter, since it represents a geometric regime choice, not a fitted bond-strength number.

**3. Add `f_accessible` to the model.** One multiplier on `n`. Don't model corona formation mechanistically; carry the uncertainty and test sensitivity in Aim 4.

**4. Aim 3 — single-marker ranking**, starting with the PANC-1/HPDEC and healthy/PDAC-plasma cases using the digitized real distributions in `ev_capture_data/` and the now-working `joint_population_yield()` function, then extending across every other system with adequate (reconstructed) data, compared against fold-change ranking. **Run every ranking under both `patch_exponent` regimes** (or at minimum, check whether the two populations being compared differ substantially in mean diameter — if not, the choice likely won't matter, per this session's finding; if so, report both numbers) rather than silently picking one.

**5. Aim 3b — multi-marker joint distributions and the independence test.** Lennon's observed EGFR clustering makes independence an empirical question.

**6. Aim 4 — sensitivity analysis**, including distribution-family choice, `f_accessible`, and a fuller version of the raw-count-vs-areal-density sensitivity check (this session covered two datasets at placeholder parameters; Aim 4 should re-run it with calibrated parameters and across the full literature table, not just the Lennon pair).

**7. (Optional, low priority) Tighten the Figure 3c gate-polygon digitization.** Currently hand-approximated as a diagonal trapezoid; a pixel-precise digitization of the dashed gate line would sharpen the exact fold-numbers reported in "The decomposition," though the qualitative conclusion is not expected to change (cross-checked against exactly-resolved points already).

**8. (Optional, worth a quick check before Aim 3) Search specifically for any published measurement of EV-surface contact-patch area vs. vesicle diameter.** The `patch_exponent=2` default rests on a general soft-matter contact-mechanics argument, not an EV-specific measurement — none was found this session, but the search wasn't exhaustive. If one exists, it should replace or validate the current approximation.

**9. Write-up, figures, poster.**

Timeline still fits. Parameter extraction is essentially complete for the primary pancreatic system, and the n-vs-density modeling question (originally expected to need Aim 3-level effort) turned out to be resolvable directly from existing digitized data and general physics literature — faster than anticipated. Remaining September time can go toward Aim 2 leg 3 and Aim 3 calibration/ranking.

---

## Known risks

- **Rankings might barely change.** Mitigation: the headline deliverable is the *relationship* between distribution shape and achievable selectivity, plus where the standard approach breaks down — informative in both directions. "The mean is a good proxy under realistic EV distributions" is a real finding, and the Poisson theorem explains why.
- **The raw-count-vs-areal-density distinction DOES complicate the story, exactly as this risk predicted — now with real numbers, not a hypothetical.** This session confirmed the model gives two disagreeing selectivity numbers for the patient-plasma comparison (PDAC captures 1.17-2.07x more under one regime, 0.18-0.30x less under the other), and this needs to be reported honestly rather than resolved by picking whichever number is more convenient. **Mitigation now in place, not just a plan:** the physical argument for the default (`patch_exponent=2`, soft-adhesion scaling) was worked out independently of which answer it would produce, and the finding that the flip is specific to size-mismatched populations (not a universal ambiguity) gives the write-up a clean, defensible way to present this as a genuine, structured sensitivity finding rather than an embarrassing inconsistency. This should be presented as a *feature* of the analysis (the model is honest about when a modeling choice matters) rather than hidden.
- **NEW — the physical argument for `patch_exponent=2` rests on general soft-matter contact-mechanics, not an EV-specific measurement.** No published study measuring EV-surface contact-patch area vs. vesicle diameter was found (search wasn't exhaustive — worth another look before Aim 3 calibration, see Immediate next steps item 8). If a future finding contradicts the soft-adhesion assumption, the `p=0` alternative and the general `patch_exponent` parameter are already built and tested, so revising the default would not require rebuilding the model — just recalibrating with a different `p`.
- **NEW — three free parameters (`K_INTRA`, `z`, `k_geom`) now need calibrating against only two measured targets (Mizenko's 10x/5x).** This is a real complication for Aim 3 that didn't exist in v3. See Immediate next steps item 2 for how to handle it (likely: fix `patch_exponent` at its physically-motivated value rather than treating it as a free fit parameter, and treat `k_geom` calibration as bundled with `K_INTRA` rather than fully independent).
- **The Lennon sign inversion — status upgraded from v2 in v3, and now built into the model in v4.** The v2 risk note said the "fewer molecules" claim was directional-only, unquantified, and dependent on the gate definition. By v3 it was quantified from digitized data, with the qualitative "opposite axes" structure surviving a cross-check restricted to exactly-resolved points. What remained uncertain in v3 — whether "raw count" or "areal density" is the physically correct quantity for the binding model — is now resolved with a physical argument and implemented as a real, tested model parameter (see Aim 1 and "What changed since v3" above). What remains uncertain now: (a) the hand-approximated gate polygon could still shift the exact fold-numbers somewhat if pixel-digitized precisely (unchanged from v3); (b) the `patch_exponent=2` default is a general physics argument, not an EV-specific measurement (see the two Known risks items directly above).
- **Do not import unmeasured "density category" numbers for any target, GPC1 or otherwise.** An external chatbot invented specific copy-number categories (e.g. "benign 1–4 copies, cancer 15–50+") to fit a story about a bulk biomarker controversy; no such single-vesicle data exists for the target in question. Treat any copy-number figure from any source as untrustworthy until traced to a specific single-vesicle measurement paper.
- **Detected vs. absolute counts.** Lennon reports detected densities; the underestimation factor is unknown. Partly cancels in ratios; state it.
- **Digitization introduces its own uncertainty.** The dense low-diameter/low-copy corner of every scatter plot is fundamentally unresolvable into individual points from a static image (ink saturates into solid blobs) — this was reconstructed via density estimation, not exact recovery, for roughly 70-90% of points in each panel depending on the panel. Validated by exact-N matching, mean/CV rescaling, visual shape comparison, and (for the headline Figure 3c finding specifically) a cross-check against only the exactly-resolved subset. State this limitation explicitly in the methods/limitations section — it is honest and defensible, not a weakness to hide.
- **Liu et al. paywalled.** Build target-agnostic; use whatever distributions are obtainable.
- **Distributions underdetermined from summary statistics** for every system other than the now-digitized PANC-1/HPDEC/plasma one. This is Aim 4, not a bug.
- **"This is just theory."** Answer: validated against published experimental signatures not used in fitting, and it makes falsifiable predictions.
- **Genuine worst case:** rankings don't move AND the trade-off curve has no interesting structure AND multi-marker analysis says nothing non-obvious. All three at once seems unlikely given how large and clearly non-random the zero-inflation is, but it isn't impossible.

---

## Practical notes

- Fully computational on public data → **no Qualified Scientist, no BSL-2**. An ISEF Adult Sponsor is still required; confirm current form requirements with the fair's Scientific Review Committee.
- ISEF requires disclosing AI use. Keep a running log: date, tool, what was asked, what was done with the answer. The research plan, abstract, poster, and citations must be in my own words. **This session (Lennon figure digitization) should be logged**: date, that Claude was used to extract per-vesicle data from uploaded PDF figures via an image-analysis pipeline (blob detection + density reconstruction), and that the resulting sign-inversion refinement was derived from that extracted data. **Start the log now if not already started** — reconstructing it in February is much harder.
- Confirm the affiliate fair's project-year window; judged data collection must fall within 12 continuous months.
- **Data provenance for the write-up:** if the digitized Figure 2/3 point clouds are used in any figure or analysis presented at the fair, the methods section should describe the digitization approach (referencing `ev_capture_data/README.md`'s method summary) and cite Lennon et al. 2019 as the source of the original figures. Judges may ask how the data was extracted from a PDF — the honest answer ("pixel-level image analysis with a documented, cross-validated reconstruction method for regions where points visually overlap") is a defensible one, not something to obscure. Preserved from v3, still current.
- **NEW — this session (n-vs-areal-density resolution) should be logged for the AI-use disclosure**, same as the digitization session: date, that Claude was used to (a) research general soft-matter contact-mechanics literature to justify a modeling choice, (b) implement and test the resulting code change, and (c) run it against the digitized data to produce the flip/no-flip finding. The physical reasoning and its limitations (no direct EV contact-patch measurement exists) should be described in your own words in the methods/limitations section, not just cited to this session.
- **NEW — a practical note for whoever starts the next chat from this document:** this session began with only `ev_capture_data/` attached; `ev_capture_model/` (the actual Aim 1 code) had to be specifically requested before it appeared. **Attach both folders from the start** next time to avoid losing a round-trip to this same issue.

