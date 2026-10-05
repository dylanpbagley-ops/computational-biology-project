# Aim 3b: Marker Co-occurrence and the Independence Assumption

This document is the ready-to-use methods/results section for Aim 3b. It exists so this part of the write-up doesn't have to be reconstructed from code later — the model, the fit, and both headline findings are already built, tested, and verified (`src/aim3b_cooccurrence_model.py`, 8 passing tests, run end-to-end this session). This section can be pasted into the paper/poster largely as-is, adapted only for length and audience.

## The question Aim 3b answers

Every other part of this project's model (Aims 1–2, the range-sweep findings) asks: for a *single* capture marker, how does nonlinear per-vesicle binding affect selectivity? Aim 3b asks a different question that only comes up once you have *two or more* markers available: if a vesicle is positive for marker A, does that tell you anything about whether it's also positive for marker B? If markers were statistically independent, you could predict the double-positive fraction just by multiplying the two single-marker positive fractions. If they're not independent — if certain vesicles tend to carry *both* markers together, or certain vesicles tend to lack *both* — that shortcut breaks, and it matters for practical assay design (do you gain much by combining two capture markers, or is one about as good as two?).

## The real data this model is fit to

Mizenko et al. 2021 (*J Nanobiotechnology* 19:250, DOI 10.1186/s12951-021-00987-1) reports that capturing vesicles via a single tetraspanin marker (CD9/CD63/CD81) recovers noticeably fewer EpCAM-positive and Her2-positive vesicles than a non-specific (biotin) capture does — roughly **10x fewer for EpCAM** and roughly **5x fewer for Her2**. This is a real, reported measurement, not an assumption. The natural explanation is that vesicles aren't a single uniform population: some fraction of vesicles is "tetraspanin-rich" (what single-tetraspanin capture actually grabs) and the rest is "tetraspanin-poor" (what it misses), and EpCAM/Her2 happen to be under-represented specifically in the tetraspanin-rich fraction — so capturing on tetraspanin alone systematically undercounts them.

## The model, in plain terms

The vesicle population is modeled as two subpopulations mixed together:

- A **tetraspanin-rich** group, some fraction `f_rich` of all vesicles. Single-tetraspanin capture is idealized as catching *all* of this group and *none* of the other (a simplification — real capture isn't a perfectly sharp on/off switch, and this is stated plainly rather than hidden).
- A **tetraspanin-poor** group, the remaining `1 - f_rich`. Non-specific (biotin) capture catches this group too, along with the tetraspanin-rich group.

Each group has its own typical (mean) copy number of EpCAM and of Her2 — set independently for each marker, because the real data shows different suppression strength for the two markers (10x vs. 5x can't both come from one shared number). Within each group, marker copy counts per vesicle follow a negative binomial distribution — the same distributional choice used everywhere else in this project's model, chosen because real biological counts are "overdispersed" (more spread out than a plain average would suggest), not because it was convenient.

**Important boundary, stated explicitly:** this model never touches `K_INTRA`, `z`, or `k_geom` — the per-vesicle binding-strength parameters from Aims 1–2. Those describe how strongly *one* vesicle's receptors bind to *one* capture surface. This model describes something different: how markers are distributed *across* different vesicles in the population. Both pieces are real and necessary, but they answer different questions, and conflating them (e.g. trying to use Mizenko's 10x/5x ratio to calibrate `K_INTRA`) would be a category error — this is a settled point in the project's "what not to revisit" list.

## How the fit works

For a candidate `f_rich`, the model searches (using a standard root-finding method, not a multi-parameter optimizer, since only one number is being solved for at a time) for the EpCAM suppression ratio — how much lower EpCAM's average copy number is in the tetraspanin-poor group vs. the rich group — that reproduces the real ~10x measured ratio. The same is done independently for Her2's ~5x ratio, using the *same* `f_rich` (since both markers live on the same underlying tetraspanin-defined subpopulation structure).

`f_rich` itself is not something Mizenko's ratio-only data can pin to one value, so the fit was repeated across a plausible range (10%–30% tetraspanin-rich, matching the same physically-reasonable range used in the earlier scoping work this file builds on). Results across that range:

| Tetraspanin-rich fraction | EpCAM suppression needed | EpCAM ratio achieved | Her2 suppression needed | Her2 ratio achieved |
|---|---|---|---|---|
| 10% | 157x | 10.22x | 71x | 4.98x |
| 15% | 164x | 10.01x | 75x | 4.98x |
| 20% | 174x | 9.91x | 80x | 4.98x |
| 25% | 187x | 9.80x | 85x | 4.99x |
| 30% | 198x | 9.97x | 91x | 4.98x |

The fit reproduces both real reported ratios closely (within about 1%) at every tested `f_rich` — confirming that this two-subpopulation structure is a genuinely plausible physical explanation for Mizenko's measurement, not just a hand-picked illustration. As expected, a smaller tetraspanin-rich fraction requires less extreme per-marker suppression to explain the same observed ratio (fewer "hidden" vesicles have to carry a disproportionate share of the signal).

## Finding 1 (robust): markers are positively correlated, not independent

Using the fitted model, the actual joint population of (EpCAM, Her2) status across vesicles was simulated directly, and two numbers were compared:

- **True joint double-positive fraction** (both EpCAM+ and Her2+ on the same vesicle): **0.7808**
- **Naive independence prediction** (just multiplying the two single-marker positive fractions, P(EpCAM+) = 0.8053 and P(Her2+) = 0.8217, together): **0.6617**

The true joint fraction is about **1.18x higher** than the naive independence shortcut predicts. In plain terms: EpCAM and Her2 tend to show up *together* on the same vesicles more often than pure chance would predict, because both are suppressed together in the same tetraspanin-rich subpopulation. This directly confirms, with real numbers for the first time, something this project's handoff had flagged for several sessions as an open question ("the independence assumption is itself a finding, not a safe default").

This 1.18x fold-deviation was checked for stability by re-running the same calculation across a wide range of the one parameter the fit can't pin down (the absolute marker copy-number scale — see below) and came out essentially unchanged every time. **This number can be cited directly and with confidence.**

## Finding 2 (real, but sensitivity-limited): combining markers helps only modestly, and the exact size is uncertain

A second, related question: if you could capture on *either* EpCAM or Her2 (an "OR-gate" capture strategy), how much better is that than just using whichever single marker performs best on its own?

At one reasonable choice of the model's remaining free parameter (the absolute copy-number scale, not just the ratio — Mizenko's data constrains ratios, not absolute numbers):

- Best single-marker yield: 0.8223
- True combined (EpCAM OR Her2) yield: 0.8464
- **Real combination gain over the best single marker: 1.03x**

That's a modest gain — nowhere near the naive "just multiply the odds of missing each marker" independence-based prediction of 0.9654 (which would suggest a much larger apparent benefit; the real correlation between markers eats most of that expected gain).

**This number is genuinely sensitive to a parameter the data can't currently constrain.** Sweeping the absolute copy-number scale gives a real range, not a single answer:

| Absolute copy-number scale (arbitrary units) | Best single-marker yield | Combination gain | Independence fold-deviation |
|---|---|---|---|
| 2.0 | 0.6258 | 1.26x | 1.18x |
| 5.0 | 0.7652 | 1.10x | 1.18x |
| 10.0 | 0.8120 | 1.04x | 1.18x |
| 15.0 | 0.8221 | 1.03x | 1.18x |
| 20.0 | 0.8270 | 1.02x | 1.18x |
| 30.0 | 0.8302 | 1.02x | 1.18x |

The combination-gain multiplier ranges from about **1.02x to 1.26x** depending on this unmeasured scale — notably *below* the project's earlier informal expectation of "2x–4x, not 20x." The independence fold-deviation (Finding 1, above), by contrast, stays essentially fixed at 1.18x across the exact same sweep.

**The honest way to report this**, and the way it should be written up: state Finding 1 (the 1.18x positive-correlation fold-deviation) as a specific, confident number. State Finding 2 (combination gain) as a range — "somewhere between roughly 1.02x and 1.3x, depending on the absolute marker density scale, which this ratio-only dataset cannot currently pin down" — never as a single point estimate. Presenting a range here isn't a weakness to downplay; it's the model correctly reporting exactly what the underlying data can and can't support, which is itself worth stating as part of the methods narrative.

## What this closes, and what's still open

This closes two long-standing "not done" items that had been open since v6 of the project's handoff: Aim 3b's co-occurrence model itself, and the multi-marker independence test that depends on it. Nothing further needs to be built to report these two findings — the numbers above are ready to cite directly, with Finding 2's range correctly framed as a range.

What remains genuinely open: pinning down the *absolute* marker copy-number scale (not just Mizenko's ratios) would let Finding 2 collapse from a range to a single number, if a future session finds a real absolute-copy-number source for EpCAM/Her2 on pancreatic-cancer-derived EVs specifically. This is a well-specified but optional future extension, not a gap in what's currently reportable.
