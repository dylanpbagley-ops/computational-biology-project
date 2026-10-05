# Aim 1 Walkthrough: Building the Multivalent EV Capture Model

This document explains, in plain language, what was built during this session, why each piece exists, the equations behind it, and one important thing that got corrected along the way. Read this alongside the code files in `src/` -- the code has detailed comments too, but this is the "big picture" version.

## Where the equations actually came from

Before writing any code, I checked the two founding papers on "superselectivity" (the phenomenon your project is built on):

- Martinez-Veracoechea & Frenkel, *PNAS* 2011
- Dubacheva, Curk, Frenkel & Richter, *JACS* 2019
- The clearest derivation of both is in a book chapter: Curk, Dobnikar & Frenkel, arXiv:1608.01222 ("Design principles for super selectivity using multivalent interactions")

My first instinct for the binding formula was `theta = 1 - (1-p)^n` (the chance that at least one of n independent coin flips comes up heads). **This turned out to be wrong for this theory** -- it's explicitly called out in the paper as the "monovalent-like" limit that does NOT show superselectivity. I caught this by checking the primary source before writing code, rather than after. The actual theory builds the binding probability from a *combinatorial* count of how many ways bonds can form, not from independent probabilities. I'm flagging this because it's exactly the kind of "sounds right but isn't" mistake you asked me to catch and admit, not bury.

## The physical picture and how it maps to your project

The papers describe a multivalent **particle** with `k` ligands binding to a **surface** with `n_R` receptors. Your project is flipped: the **surface** has the ligands (density `sigma`), and the **vesicle** has `n` receptors. The math doesn't care which side is called "particle" vs. "surface" -- what matters is picking the right *disorder regime*, since the paper gives three different versions depending on what's allowed to move around:

1. Fixed ligand and receptor counts, all can reach all -- exact combinatorial sum
2. Receptor count is Poisson-distributed (mobile receptors), ligand count fixed
3. **Both** receptor and ligand counts are Poisson-distributed (both mobile) -- this is the "large colloid" case in the paper's own words

I chose regime 3, because your handoff's mechanism (a) -- lateral diffusion lets receptors reshuffle into the contact patch during the ~15 minutes to overnight incubation -- means your EV's receptors really are effectively randomly repositioning, and a "large colloid" (an EV, 30-150 nm, versus nanometer-scale ligands) meeting a surface with many independently-positioned ligands is exactly the physical picture this regime describes.

## The equations, translated into our symbols

**Bond-count partition function** (how many effective ways bonds can form):

```
q_b(n, sigma) = exp(n * sigma * K_INTRA) - 1
```

- `n` = number of receptors on this specific vesicle
- `sigma` = ligand density on the capture surface (your design variable)
- `K_INTRA` = a single lumped number standing in for bond strength, contact patch area, and other geometric factors. This is the "knob" we'll eventually calibrate against real data (Aim 3) -- for now it's just a placeholder number.

**Converting that into an actual probability** (the Langmuir wrapper):

```
theta(n, sigma) = z * q_b / (1 + z * q_b)
```

- `z` = a separate dimensionless number capturing everything about *getting the vesicle to the surface in the first place* -- concentration, incubation time, diffusion -- as opposed to the multivalent bond-counting itself.

The important thing to notice: `n` sits **inside an exponential**. Doubling the receptor count doesn't just double your chances -- it can multiply your capture probability by a large factor, because more receptors means combinatorially more chances for a bond. That's the whole nonlinearity this project is about.

## What was built and validated

| File | What it does |
|---|---|
| `src/binding_model.py` | The analytical (closed-form) single-vesicle model: `q_bound()`, `capture_probability()`, `langmuir_probability()`, `selectivity_alpha()` |
| `src/monte_carlo.py` | A brute-force simulation of the same physics (see below for why this counts as an independent check) |
| `src/compare_analytical_vs_mc.py` | Runs both side by side across many receptor counts and confirms they agree |
| `src/population_model.py` | Extends the single-vesicle formula to a whole population with a distribution of receptor counts -- this is Aim 1's core novelty |
| `src/demo_population_effect.py` | A demo showing two populations with the *same average* receptor count but different shapes give meaningfully different capture yields |
| `tests/test_aim1.py` | The formal unit tests for Aim 2's required limiting cases |

### Validation 1: analytical vs. Monte Carlo agree

The Monte Carlo simulation doesn't use the `q_bound()` formula at all -- it simulates the underlying random process directly (drawing a random number of bonds from a Poisson distribution with mean `n * sigma * K_INTRA`, since that's the physical assumption the formula encodes) and checks whether at least one bond formed, repeated 200,000 times per data point. Comparing this to the closed-form answer, the two matched everywhere to well within statistical noise (worst case was about 1.7 standard errors out of thousands of trials, which is unremarkable). See `figures/analytical_vs_mc.png`.

### Validation 2: the two required limiting cases

**Poisson mean-field theorem.** You proved earlier that mean-field theory is exact for Poisson-distributed receptors. While building the test for this, I found the exactness lives at a more specific level than I first assumed:

- `E[q_b(n)]` (averaging the *bond-count partition function* over a Poisson population) matches its closed form `exp(mean_n * (exp(sigma*K_INTRA) - 1)) - 1` to machine precision -- this is genuinely, provably exact, confirmed by the Poisson moment-generating function.
- But converting that into an actual capture probability requires one more step (the Langmuir wrapper), and that step is *also* nonlinear. Averaging the full capture probability over the population is only approximately equal to plugging the averaged `q_b` into the Langmuir formula -- off by roughly 1-5% in the ranges tested, more at extreme values of `z`.

I flagged this distinction to you directly mid-build rather than writing a test that would silently pass on a looser claim than what's actually true. You asked me to test at the level where it's genuinely exact (the `q_b` level) and document the second, smaller approximation clearly -- both are now in the code comments and in `test_aim1.py`.

**Valency = 1 recovers Langmuir binding.** Similarly, I found this only holds *exactly* in the dilute/weak-binding limit (low `sigma * K_INTRA`), not at every ligand density -- because our exponential formula `exp(x) - 1` only equals the linear Langmuir form `x` approximately, for small `x`. The gap is about 5% at moderate density and shrinks proportionally as density drops, exactly as expected from the math (this is the first term of a Taylor series). You confirmed this is a known, acceptable property of the chosen physical regime rather than a bug, so the test checks for convergence in the dilute limit rather than exact equality everywhere -- and documents why.

Both corrections are the kind of thing you said you wanted acknowledged directly rather than smoothed over, so I want to make sure they're visible here, not just buried in code comments.

### The core novelty demonstration

`demo_population_effect.py` builds two receptor-count populations with the **exact same average** (20 receptors), but different shapes: one Poisson (the "textbook" shape, no zero-inflation), one negative-binomial (heavily zero-inflated -- about 16% of vesicles have zero receptors, closer to what real EV data looks like per your handoff). Feeding both into the population-average yield calculation across a range of ligand densities shows the two populations produce yields that differ by up to **28% at the same mean receptor count**. This is a working, numeric demonstration of the project's central claim: the mean alone is not sufficient information, and fold-change ranking (which only ever looks at the mean) can be misled by distribution shape. See `figures/population_shape_effect.png`.

## What's NOT done yet (next steps, not gaps)

This session covered Aim 1 and the first half of Aim 2 (limiting-case validation). Still ahead, per your timeline:

- Reproducing the specific experimental signatures in your Validation Targets table (Mizenko et al.'s ~10x/~5x underperformance numbers, the DNA nanostar optimal valency of 3) -- this needs literature parameter extraction, which you've scheduled for November.
- Calibrating `K_INTRA` and `z` against real or reconstructed distributions (Aim 3).
- Multi-marker capture and the independence assumption (Aim 3's core sub-aim).
- Sensitivity analysis (Aim 4).

## A note on how I worked

I deliberately verified the source equations against the primary literature *before* writing code, rather than coding from memory and hoping it was right. That surfaced the `1-(1-p)^n` mistake early, when it was cheap to fix, instead of after building several more aims on top of a wrong foundation. I also stopped and asked you directly at both points where the math turned out to be subtler than the handoff implied, rather than picking a resolution myself and mentioning it as a footnote. Both of those pauses are documented above and in the code.

---

## Session addendum: resolving n (total count) vs. areal receptor density

This addendum documents work done in a later session, to resolve "Immediate next steps" item 1 from the v3 project handoff: whether the model's receptor-count variable `n` should represent total receptors on the vesicle, or receptors per unit membrane area.

### The physical question

The original model (everything above this addendum) feeds raw receptor count `n` straight into `q_bound()`, with no reference to vesicle size. This was fine until the digitized Lennon Figure 3c data showed that PDAC-enriched plasma vesicles are ~3x larger in diameter than typical healthy vesicles, with total receptor count moving in the OPPOSITE direction from areal density between the two populations. The model needed a way to represent this distinction.

### The resolution: a single new geometric parameter, not two competing models

Rather than picking either "total count" or "areal density" and hard-coding it, `binding_model.py` gained a `contact_patch_receptors()` function that computes the number of receptors actually available for bonding inside the vesicle's contact patch, as a function of total count `n`, diameter `D`, and a new **contact-patch scaling exponent `p`**:

```
n_contact = k_geom * n * D^(p - 2)
```

Two named special cases:
- **`p = 2`** (soft-adhesion limit, the default): diameter drops out entirely, `n_contact` is proportional to total count `n` alone. This is mathematically IDENTICAL to the original model (verified by unit test to floating-point precision) -- it just makes explicit a geometric assumption that was always implicit.
- **`p = 0`** (fixed-patch limit): `n_contact` becomes proportional to areal density alone. Diameter enters explicitly.

The `p = 2` default is not arbitrary: it follows from soft-matter adhesion physics for very soft, small bodies (contact radius scales linearly with the body's own radius in this "liquid" regime, rather than the sub-linear scaling that applies to stiffer bodies), which is the appropriate regime for a floppy lipid EV (bending rigidity 20-30 kT, per Jahnke et al. 2024) meeting a much stiffer functionalized capture surface. This is a real, sourced physical argument, but it is still an approximation -- no direct measurement of EV contact-patch size vs. diameter was found in the literature, and this is flagged honestly in the code.

### The result: it depends on how much the compared populations differ in size

Running both regimes on the real digitized data (`src/patch_regime_sensitivity.py`) gives a genuinely important finding:

- **Patient plasma (Figure 3c, healthy vs. PDAC-gated, ~3x diameter difference):** the ranking FLIPS. Size-blind (p=2) says PDAC-enriched vesicles capture 1.17-2.07x MORE than healthy; density-driven (p=0) says they capture 0.18-0.30x LESS.
- **Cell lines (Figure 2, PANC-1 vs. HPDEC, diameters nearly matched, ~1x ratio):** the ranking does NOT flip. Both regimes agree PANC-1 (tumor) captures 1.24-1.51x more than HPDEC (normal).

This makes sense physically and isn't a coincidence: when two populations have similar diameters, dividing receptor count by a nearly-constant D² barely changes anything, so the choice of `p` barely matters. The choice only matters when vesicle SIZE itself differs substantially between the populations being compared -- exactly the situation in the gated plasma data, and not the situation in the cell-line data. Checked for robustness across a 3x3 grid of `sigma`/`k_intra` placeholder values; the plasma-data flip holds in all 9 combinations tested, so it is not an artifact of one parameter choice.

### A units subtlety worth flagging

The lumped geometric constant `k_geom` has different physical units under `p=2` (dimensionless, a fraction of membrane area) versus `p=0` (units of length², an absolute patch area). Using the same numeric value for both regimes produces a meaningless comparison (this was tried first and produced near-zero yields under `p=0` -- a units artifact, not a physics result). The fix: calibrate `k_geom` separately per regime so both agree on a shared reference point (the average normal/healthy vesicle gets exactly 1 contact receptor under both regimes), which is documented in `patch_regime_sensitivity.py`.

### What this means going forward

This is now flagged as a live open question for Aim 3/Aim 4 rather than a settled implementation detail: the model can be run in either regime (or values of `p` between them), and target rankings should be reported under both, with the disagreement itself reported as a finding when it occurs. New functions added: `binding_model.contact_patch_receptors()`, `binding_model.q_bound_patch()`, `binding_model.capture_probability_patch()`, `population_model.joint_population_yield()`, `population_model.compare_patch_regimes()`. Four new unit tests added to `tests/test_aim1.py` (8 total, all passing), plus `src/patch_regime_sensitivity.py` as the standalone analysis script and `figures/patch_regime_sensitivity.png` as the summary figure.

---

## Session addendum: Aim 2 leg 3 -- external validation (nanostar valency-3, and finally exercising `selectivity_alpha()`)

This addendum documents a later session's work on the v4 handoff's top-priority "not started" item: reproducing the two external experimental signatures in the Validation Targets table, specifically the DNA nanostar optimal-valency-3 finding (Linne et al. 2024) and exercising `binding_model.selectivity_alpha()` (written since Aim 1 but never tested) against a real theoretical signature.

### An honest complication surfaced immediately, and was brought to you before proceeding

The nanostar valency-3 optimum, on close reading of the actual paper (Linne, Heemskerk, Zwanikken, Kraft & Laan, *Soft Matter* 2024, 20:8515-8523, DOI 10.1039/d4sm00704b, open access via PMC11484159), turns out to come from mechanics our EV model does not have: a fixed, small number of discrete binding arms (valency `k`, a real design parameter) meeting VARIABLE receptor density, where the optimum arises because higher `k` means faster binding from solution, which saturates the surface before the interesting superselective transition becomes visible at experimentally accessible receptor densities. Our EV model has the geometry flipped (vesicle with variable receptor count meeting FIXED ligand density) and has zero surface-saturation/depletion mechanics -- every vesicle in our model captures independently, with no cap on total capture. Reusing `binding_model.py` under relabeled variables would not have tested anything real about this mechanism. This was flagged directly before any code was written, and the decision (build the real mechanics as a separate module) was yours, not assumed.

### What was built: `src/nanostar_validation.py`, a genuinely separate model

This new file implements Linne et al.'s own published equations directly (their Eq. 2-4): a monovalent Langmuir isotherm for `k=1`, and an avidity extension for `k>1` where `K_A` (first-bond strength) is held roughly constant across valencies (matching what the paper itself reports finding) while `K_intra` (subsequent-bond strength) decreases with valency (also matching the paper's reported trend) -- the trade-off between the combinatorial `(k-1)` boost from more arms and this per-arm weakening is exactly the mechanism the paper's Discussion section identifies as the source of the k=3 optimum.

**An honest limitation, encountered and handled in real time, not smoothed over:** the paper's own Figure 2(e) table -- containing their actual fitted `K_A`/`K_intra` numbers for k=1,3,6,10 -- did not come through as machine-readable text from either the PMC full-text API or the RSC-hosted PDF (both were checked; the table appears to be an embedded graphic in both sources). Rather than inventing numbers and presenting them as the paper's own, or giving up on the validation, the approach taken was: use the paper's own *qualitative* statement of the relationship (K_A roughly constant, K_intra decreasing with k) to choose representative parameters, and report the result honestly as "mechanism-faithful, not digit-faithful" -- a real, reproducible confirmation that the described trade-off mechanism CAN produce an interior valency optimum, not a claim of an exact quantitative match to their published numbers.

**Getting a working parameterization took real iteration, documented rather than hidden:** the first parameter guess (K_A=1.0, mild K_intra decrease) produced NO optimum at all -- max selectivity was highest at k=1, which was a sign the model wasn't even in the genuinely dilute regime the theory requires (the same "must be dilute enough to see the effect" condition that matters for `selectivity_alpha()` itself, see below). Fixing this required first recognizing that K_A needed to be small enough that the system starts near-unbound at low receptor density, then finding a K_intra(k) decay steep enough to overcome the combinatorial boost larger k automatically gets. The final working parameterization (K_A=0.01, K_intra(k) ~ k^-3.5) does reproduce the k=3 optimum -- but a further honest check found this is **not robust to K_A**: a ten-fold weaker K_A shifts the optimum to k=6 instead. This is reported directly in the script's own output (not just in a comment), and is treated as consistent with, not contradicting, the paper's own point that valency, binding strength, and concentration jointly determine the optimum.

**A second honest limit, found while writing the test suite:** at the representative baseline parameters, k=10 doesn't exceed alpha=1 at all -- it actually falls below the plain Langmuir ceiling. Checked systematically: no tested K_intra(k) decay rate makes k=3 win over k=6 and k=10 while ALSO keeping k=10 above alpha=1; the two requirements trade off. This is plausible physically (an over-decorated particle whose extra arms barely engage can behave worse than a clean single-site binder) but is a property of this specific representative parameterization, not a claim that real k=10 nanostars are sub-Langmuir (Linne et al.'s own data shows they aren't, just less selective than k=3 or k=6). The test suite checks only that k=3 -- the paper's actual reported winner -- exceeds 1, not that every multivalent k does.

### A genuinely new finding while finally exercising `selectivity_alpha()`

`binding_model.selectivity_alpha()` had existed since the first session but had never been tested or run against real numbers -- flagged explicitly in the v4 handoff's "Not done" list. Running it for the first time surfaced something worth knowing: **at z=1.0 -- the activity/concentration value used in every existing demo script in this codebase** (`compare_analytical_vs_mc.py`, `demo_population_effect.py`, `patch_regime_sensitivity.py`) **-- alpha never exceeds 1 anywhere**, meaning none of the existing demo figures are technically showing superselective behavior by the strict Martinez-Veracoechea/Frenkel definition (alpha > 1), even though the underlying exponential nonlinearity is still doing real work in those figures. This is not a bug: it traces to the same "must be in a sufficiently dilute regime" condition that the founding superselectivity papers state explicitly, and that showed up independently while building the nanostar reproduction above. Using a much smaller z (e.g. 1e-4) with the same sigma and k_intra recovers the textbook superselective shape exactly: alpha starts near 1 at low n, rises to a real peak (~6, in one tested case) at intermediate n, then decays toward 0 as capture probability saturates at large n. This is now documented directly in `selectivity_alpha()`'s own docstring, and is worth remembering for Aim 3/4: whichever way `z` ends up being calibrated (dilute vs. concentrated conditions) will determine whether the superselectivity mechanism is actually operating in the final calibrated model, or just present as unused machinery.

### What was added

- `src/nanostar_validation.py` -- new, separate model implementing Linne et al.'s own equations; produces `figures/nanostar_valency_validation.png` when run directly.
- A new comment block added to `binding_model.selectivity_alpha()`'s docstring documenting the z-dependence finding.
- `tests/test_aim2_leg3.py` -- seven new tests: three exercising `selectivity_alpha()` for the first time (dilute-regime requirement, low-n sanity bound, numerical stability near n=0), and four validating the nanostar reproduction (monovalent ceiling holds generally, k=3 achieves genuine superselectivity, k=3 is the winner at representative parameters, and the K_A-sensitivity/k=10-sub-Langmuir limitations are locked in as documented, checked facts rather than left as one-off script output).
- All 15 tests (8 original + 7 new) pass via both `python3 test_*.py` and `pytest`.

### What this does and doesn't accomplish for Aim 2

This should be read as a genuine, defensible qualitative/mechanism-level external validation -- the trade-off mechanism the paper describes really can produce an interior valency optimum matching their reported k=3, using our own independent implementation of their equations -- but not a precise quantitative match to their published fit values, which were not recoverable from the available sources this session. The Mizenko et al. 10x/5x numbers (the other two rows in the Validation Targets table) remain correctly assigned to Aim 3's `K_INTRA`/`z`/`k_geom` calibration step, not to this validation leg -- Mizenko's paper reports fold-underperformance numbers, not a superselectivity alpha, so it was never the right target for exercising `selectivity_alpha()` against.
