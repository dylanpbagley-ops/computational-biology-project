# Limitations: Why K_INTRA, z, and k_geom Are Not Calibrated, and Why That's OK

This document exists so the write-up's limitations section doesn't have to be improvised at the last minute. It answers three questions plainly: what isn't calibrated, why six full sessions of searching never found a way to calibrate it, and what that does (and doesn't) cost the project. It also states, formally, a mathematical fact about the model that explains *why* the search kept failing in the specific way it did -- this second part came out of an independent external review of the project (see below) and sharpens something this project already suspected into something now actually proven.

## What "uncalibrated" means here, precisely

The model has three lumped parameters that never got a single fitted numeric value:

- **K_INTRA** -- how strongly a single receptor-ligand bond forms, folded together with the contact patch area and some geometric constants (see `binding_model.py`'s module docstring for the full derivation).
- **z** -- a catch-all "activity" number standing in for vesicle concentration, incubation time, and how easily a vesicle diffuses to the capture surface in the first place.
- **k_geom** -- a geometric constant relating total receptor count to how many of those receptors actually sit inside the vesicle's contact patch with the surface.

None of these were ever fit to a single real number. Instead, `k_intra_z_plausibility_bounds.py` builds a physically-grounded *range* that each parameter should plausibly fall in, anchored to a real, measured antibody surface density (from the Lopez Baltazar dissertation -- see "Key literature" in the handoff). Every result this project has produced (the ranking-flip findings, the superselectivity-reachability finding, the stem-cell marker crossover) was computed as a sweep across this range, not at one fitted point.

## Why calibration was never found, despite six full search sessions

Six separate sessions (spanning roughly two months of this project's timeline) searched the literature for a data source that could pin these parameters down. Every search used a different strategy: engineered EV platforms with tunable receptor counts, virus-like particles, synthetic liposomes with programmable ligand density, PAMAM dendrimer/DNA systems, and a controlled-surface-density EV biosensor. All six came back empty, but they came back empty in the *same specific way*, which is itself worth stating plainly because it's informative:

**Every real candidate source splits the two things a calibration needs -- a genuinely controlled/known particle-side receptor (or ligand) count, and a genuine surface-immobilized (not solution-phase, not cell-suspension) binding assay -- across different experiments, different papers, or different analytes.** Papers that control and measure receptor copy number precisely (e.g. HaloTag-engineered EVs) don't run a surface-capture binding assay at all. Papers that run a real surface-immobilized binding assay with independently verified surface density (e.g. the IRIS biosensor paper checked in the sixth search session) explicitly treat the EV-side receptor count as unknown and build their whole analysis around not needing it. Nobody appears to have published a study that puts a particle of known valency onto a receptor-density-known surface and measures the resulting binding curve, for any EV-relevant system.

This could mean the field genuinely hasn't done this experiment yet (EV engineering and biosensor/SPR work seem to be pursued by largely separate research communities with different goals), or it could mean the right paper exists but uses different terminology than what was searched. Both possibilities are stated here rather than picking one, because neither can be ruled out.

Given this, the decision (made explicitly by the project's author, not assumed) was to close the search permanently rather than keep spending project time on a seventh, eighth, or ninth attempt, and to make the range-sweep methodology -- already demonstrated on three independent real-data-grounded systems -- the project's actual answer to the calibration gap, rather than treating the gap as a placeholder waiting to be filled later.

## The formal reason calibration was never going to work from a single experiment, stated precisely

An independent external review of this project (a red-team-style assessment from an AI research tool called Biomni, dated 26 August 2026) pointed out something this project's own sessions had described informally but never proven: `sigma`, `K_INTRA`, and `k_geom` are not merely *uncalibrated* -- they are **not separately identifiable from any single-condition experiment, even in principle.** This was checked directly, algebraically and numerically, rather than taken on faith from the report.

Here's what that means and why it's true. Look at the core equation:

```
q_bound(n, sigma, K_INTRA) = exp(n * sigma * K_INTRA) - 1
```

Only the *product* `n * sigma * K_INTRA` ever enters this formula. If you double `sigma` and simultaneously halve `K_INTRA`, the product `sigma * K_INTRA` is exactly unchanged, so `q_bound` -- and therefore every downstream quantity (capture probability, yield, selectivity) -- comes out bit-for-bit identical. This was checked directly in the actual code this session: holding `sigma * K_INTRA` fixed while rescaling the two in opposite directions by a factor of 37.5 produces `q_bound` values that agree to the last bit (difference exactly 0.0, not just "very small").

The same thing happens with `k_geom` and `K_INTRA`, and this turns out to hold more broadly than first expected. The general contact-patch formula is:

```
n_contact = k_geom * n * diameter^(patch_exponent - 2)
q_bound_patch = exp(n_contact * sigma * K_INTRA) - 1
```

`k_geom` enters this formula as a simple linear multiplier *regardless of `patch_exponent`* -- the `diameter^(patch_exponent-2)` term changes how vesicle diameter enters, but never changes how `k_geom` and `K_INTRA` combine. So the combined product `k_geom * K_INTRA` is the only thing that matters at *any* `patch_exponent`, not just the `patch_exponent=2` default. This was checked directly this session at both `patch_exponent=2` (factor of 5, matching to floating-point precision) and `patch_exponent=0` (also checked, also matching) -- confirming the confound is `patch_exponent`-*independent*, a real, slightly non-obvious property of the model worth stating precisely rather than assuming it only applies to the default regime.

**A numerical note worth carrying forward, found while writing the regression tests for this section:** naively checking these rescalings with an exact bit-for-bit equality check (`q_bound` computed directly, compared for exact equality) gives a spuriously large "failure" at big rescaling factors -- not because the identifiability claim is wrong, but because `q_bound = exp(x) - 1` is exponentially sensitive to `x`, so a floating-point rounding error in `x` far too small to see when printing it (order 1e-16, from the rescaling multiplication itself) becomes a huge absolute difference in `exp(x)` once `x` is large (at `x≈133`, `exp(x)≈7×10^57`, so even a relative error of 1e-16 in `x` becomes an absolute error of order 1e41 in `exp(x)`). The correct way to check this numerically is in log-space (comparing `x` itself, or equivalently `log(q_bound+1)`), which is linear and doesn't amplify the floating-point noise -- this is the same general "large dynamic range needs relative, not absolute, comparison" lesson this project's `aim3_aim4_range_sweep.py` already documents for a related reason (`MIN_RELIABLE_EXPONENT_AT_N1`). See `tests/test_identifiability.py` for the actual corrected check.

**What this means in plain terms:** no single experiment -- no matter how well-designed, no matter how precisely measured -- could ever separately tell you the individual values of `sigma`, `K_INTRA`, and `k_geom`. It could only ever tell you the *combined product* of some of them. Separating the individual factors would require multiple, independent experiments that vary each one while holding the others fixed -- exactly the "known particle-side valency AND known surface density, in the same experiment" combination that six sessions of searching never found published anywhere. This isn't a limitation of this project's search skill; it's a structural fact about what kind of data the model actually needs, which is a stronger and more useful statement than "we couldn't find it."

This also means the negative search result and this identifiability argument reinforce each other rather than being two separate problems: the reason the search kept failing in exactly the same way (controlled valency and controlled surface density always show up in different papers) is *because* only an experiment that controls both simultaneously could ever separate these parameters -- an experiment that controls only one of them, however good, was never going to be enough even if it had been found.

## What this costs the project, and what it doesn't

The impact is narrow and specific, not diffuse:

**Stays available:** any claim of the form "nonlinear capture reorders these targets relative to fold-change ranking, across roughly X% of the physically plausible parameter space." This is exactly the kind of claim every current finding already makes -- the ranking-flip results, the superselectivity-reachability result, and the stem-cell marker crossover were all computed as range-sweeps from the start, specifically anticipating that a single calibrated number might never arrive. Nothing already built needs to be redone, retracted, or reframed because of this limitation.

**Stays unavailable:** any claim of the form "for a real anti-CD81 capture assay under these specific conditions, selectivity will be exactly N-fold." Single-number point predictions remain conditional on "somewhere within this several-order-of-magnitude range" rather than collapsing to one answer.

For a solo, computational, ISEF-level project, reporting a well-justified range instead of an arbitrary or over-confident point estimate is a legitimate scientific choice, not a consolation prize -- and it is a stronger position to be in than being asked hard calibration questions by a judge and improvising an answer live, since the honest answer (here is exactly why this can't be calibrated from any single experiment, here is the range instead, here is why the range methodology is the deliberate response to that fact) is already worked out.

## What to actually say in the write-up

State plainly, in the limitations section:

1. `K_INTRA`, `z`, and `k_geom` are physically bounded (via a real measured antibody surface density) but not fitted to a single value.
2. Explain why in the identifiability terms above: these parameters are not separately identifiable from any single-condition experiment, because only certain *products* of them (`sigma * K_INTRA`, `k_geom * K_INTRA`) ever enter the model's predictions -- this is a provable property of the equations, confirmed algebraically, not just an informal caveat about missing data.
3. State that six full search sessions looked for a data source that varies both particle-side valency and surface density independently in the same experiment, and found none -- with the informative structural pattern (every real candidate splits the two ingredients across different papers/experiments/analytes) stated as a real, if inconclusive, finding about the current state of the field.
4. Frame the range-sweep methodology as the deliberate, designed response to this identifiability gap, not an apology for a missing number -- and note it has now been demonstrated on three independent, real-data-grounded systems (the two Lennon-based comparisons and the stem-cell marker-ranking extension).
5. Audit prediction language throughout the rest of the write-up to make sure every claim is phrased as a range statement, not a point prediction, consistent with what the underlying computations actually support.
