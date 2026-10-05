# Aim 4: Distribution-Family and Plausibility-Bounds Assumption Sensitivity — Build Log

**Session:** the session immediately following the v16 restructuring session's same-day follow-ups (still August 26, 2026). Picked up `EV-Capture-CURRENT.md`'s "Immediate next steps" item 1: "Aim 4 — the rest of the sensitivity analysis. Distribution-family choice, a fuller sensitivity analysis of the plausibility bounds' own assumptions beyond the single widened check already done."

## What "the single widened check already done" refers to, precisely

This phrase does not appear literally in any source file — it refers to the `target_exponents` diagnostic widening run across the v9→v10 and v10→v11 handoff transitions (archived in `EV-Capture-Project-Handoff-v15.md`, "What changed since v10" / "since v9"). That check widened `dimensionless_exponent_bounds()`'s `target_exponents` ceiling from `(0.01, 0.1, 1.0)` to include `5.0, 10.0, 20.0`, found the qualitative superselectivity finding unchanged, and the user explicitly decided to leave the official bound as-is and treat the widened check as a write-up footnote only — a closed decision (see "What NOT to revisit"). This session's task was the two sensitivity questions that check was never meant to answer: distribution-family choice, and the plausibility bounds' *other* unexamined assumptions.

## Scoping: where distribution-family choice is (and isn't) a live question

Before writing any code, the two families of headline findings were checked directly to see which one distribution-family sensitivity could even apply to:

- **The Lennon-based real-data findings** (CA19-9 flip, retired-EGFR flip) never fit a distribution at all. `population_model.joint_population_yield()`'s own docstring is explicit: it takes real per-vesicle `(n, diameter)` pairs directly from a digitized CSV and averages over them — "no explicit P(n) needs to be fit first." Distribution-family choice is therefore **not a live sensitivity axis** for these findings.
- **The stem-cell CD81/CD63/CD73 marker-ranking crossover** (`stem_cell_marker_ranking_extension.py`) is the *only* headline finding built from a reconstructed population: a negative binomial fitted to two summary statistics (fraction positive, mean copies if positive) per marker/tissue, via a 1D root-find on the dispersion parameter. That file's own module docstring already states the project's philosophy on this ("assuming a distribution family and reporting sensitivity to it is more honest than fitting one shape and presenting it as truth") but had never actually delivered the sensitivity check. This is where the real work was.

The same file also has a second, related, never-tested assumption: every marker/tissue is assigned the *same* generic lognormal diameter distribution (120nm mean, CV=0.3), explicitly flagged in its own docstring as an assumption rather than a measurement.

Separately, `k_intra_z_plausibility_bounds.py`'s `k_intra_bounds_from_real_sigma()` bakes in two more assumptions besides `target_exponents`: the assumed receptor-count grid (`n_values=(1,10,50,100)`), and the sigma anchor itself (the Lopez Baltazar "40% coverage" operating point, 4.44×10¹¹ Abs/cm², vs. the same dissertation's own 100%-close-packed ceiling, `SIGMA_MAX_LOPEZ_BALTAZAR_ABS_PER_CM2` = 1.11×10¹² Abs/cm² — a real, already-defined value that no prior session had ever actually swept against).

## New module: `src/aim4_distribution_and_assumption_sensitivity.py`

Four checks, each run to completion and independently re-verified by inspecting raw numbers (not just trusting a printed summary), per this project's standing discipline:

**1. Zero-inflated Poisson as an alternative to negative binomial.** A zero-inflated Poisson (ZIP) is a structurally different generative story from negative binomial — a hard Bernoulli on/off switch times a Poisson draw, vs. negative binomial's continuous Gamma-mixture-of-Poissons — even though both can be matched to the same two summary statistics. Because ZIP has exactly two free parameters and exactly two targets, the fit is a direct algebraic construction (`zero_inflation_prob = 1 - fraction_positive`, `poisson_rate = mean_copies_if_positive`), not a root-find. This was checked empirically, not just algebraically: a 200,000-sample draw at every one of the 9 real marker/tissue combinations recovers both target statistics to within 0.01 (fraction) and 0.5 (mean), and the two families' overall `mean_n` values agree to ~0.008% (the tiny residual traced to the negative-binomial fit's own root-find not landing exactly on the naive product — a real, small property of the existing negative-binomial code, not a bug in the new ZIP code; confirmed directly by checking `scipy.stats.nbinom.pmf(0, r, p)` against the target zero-fraction, which matches to ~1e-9).

**Result: the crossover survives.** CD81 leads the marker-level ranking at the largest plausibility-range K_INTRA under the ZIP reconstruction, at all 4 of the same seed offsets (0, 100, 500, 1000) already used to establish the negative-binomial version's seed-robustness. This is a genuine, new robustness finding — the crossover is not an artifact of the specific choice of negative binomial.

**2. Diameter-model sensitivity.** Four variants tested against the same 4 seeds: lognormal at CV=0.3 (existing default), CV=0.6, CV=0.15, and a gamma distribution (a genuinely different tail shape) matched to the same mean and CV=0.3. **Result: the crossover survives all four**, at every seed.

**3. Sensitivity to the assumed receptor-count grid.** Three grids tested for `k_intra_bounds_from_real_sigma()`'s `n_values` parameter: the default `(1,10,50,100)`, a low-shifted `(1,5,25,75)`, and a high-shifted `(5,20,100,300)` (deliberately extended past the documented 1–100 range to test a genuinely wider assumption). **Result: `any_superselective=True` under all three**; max_alpha ranges 10.44–10.46, fraction of the grid superselective ranges 0.577–0.679 (some variation in the *fraction*, but the qualitative "superselectivity is reachable" finding does not flip under any tested grid).

**4. Sensitivity to the sigma anchor.** The default 40%-surface-coverage operating point vs. the same dissertation's own 100%-close-packed ceiling (`SIGMA_MAX_LOPEZ_BALTAZAR_ABS_PER_CM2`). **Result: max_alpha is identical to machine precision (10.461142... both ways)** — exactly as predicted by the identifiability proof in `Limitations_and_Calibration_Gap.md` (only the product sigma×K_INTRA enters the model, and `k_intra_bounds_from_real_sigma()` re-derives K_INTRA consistently against whichever sigma is passed in, so the swept product grid is unchanged by construction). This is a nice independent confirmation, via a genuinely new code path, of a fact the identifiability proof already guaranteed algebraically — worth noting in the write-up as a consistency check that passed, not just a redundant re-statement.

## Tests

`tests/test_aim4_distribution_and_assumption_sensitivity.py`, 6 new tests, all passing:
1. ZIP fit hits both targets empirically for every marker/tissue combination.
2. ZIP and negative-binomial `mean_n` agree to within the small tolerance the negative-binomial root-find itself introduces (documented above).
3. Crossover survives the ZIP family (all 4 seeds).
4. Superselectivity reachability stable across the 3 tested `n_values` grids.
5. Superselectivity reachability stable across the 2 tested sigma choices, plus an exact-match check on max_alpha between them.
6. Crossover survives all 4 diameter-model variants (all 4 seeds each).

Full suite: **87 tests, all passing** (81 carried forward + 6 new). Every affected script's own `__main__` entry point (`aim3_aim4_range_sweep.py`, `stem_cell_marker_ranking_extension.py`, `k_intra_z_plausibility_bounds.py`, `binding_model.py`, and the new module itself) was run directly (`python3 <file>.py`), not just via pytest — all exit cleanly with no errors, per this project's standing "a passing test suite is not sufficient" discipline.

## What this closes, and what it doesn't

This closes "Immediate next steps" item 1 from `EV-Capture-CURRENT.md` as originally scoped: distribution-family choice and the plausibility bounds' other assumptions (beyond the already-closed `target_exponents` check) have now been checked, and every headline finding tested survives every alternative tried.

**What this does NOT do:** it does not test every conceivable alternative distribution family (e.g., a beta-binomial, or a genuinely different diameter-copy-number correlation structure) — the two families chosen (zero-inflated Poisson vs. negative binomial; lognormal vs. gamma) were picked specifically for being structurally different generative stories, not just different-shaped curves, so a "survives both" result is a meaningful robustness statement, but it is not an exhaustive proof that no distributional assumption could ever change the ranking. This should be stated in the write-up's limitations section in the same spirit as the other honestly-scoped findings in this project (a real, checked robustness result, not an unbounded claim).
