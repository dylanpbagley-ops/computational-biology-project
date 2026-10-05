# Session build log: Gillespie real units, Aim 3b write-up, f_accessible real-scenario application

(The session immediately following the v15→v16 restructuring, August 26, 2026. This is an archive file per `PROCESS.md`'s consolidation instructions — `EV-Capture-CURRENT.md` points here for the full detail; its own "What changed this session" section carries only a short summary now that this session is no longer the most recent one.)

## Context

This session completed items 1, 2, and 3 from the "Immediate next steps" list left by the v15→v16 restructuring session: pin down the Gillespie cooperative model's real physical units, write the Aim 3b real-data write-up, and apply `f_accessible` to a real scenario. All three were genuinely open, well-specified items — no discrepancy requiring a user check-in this time (unlike the prior session's step-4 mismatch).

## 1. The Gillespie cooperative model's real physical units

Retrieved Linne et al. 2024's actual reported Methods/Results numbers via the paper's RSC open-access HTML full text (fetched directly, not summarized secondhand): `rho_n = 1e-8 M` (their real fixed nanostar concentration for the valency sweep) and `sigma_R` swept over `1,000–600,000 um^-2` (their real reported range), plus `A = 22 mm^2` (the flow-cell area used to make their K_A dimensionless-compatible).

The more important discovery: **Table 2's `K_A`/`K_intra` are NOT the values that actually produce the paper's k=3 optimum.** Table 2 is a separate per-valency experimental fit (shared K_A/K_intra forced across k=3,6,10 specifically to isolate Delta_G_c's effect). The values that DO produce the k=3 optimum, per the paper's own Fig. 4c–e caption, are three illustrative (K_A, K_intra) pairs, with panel (e) — K_A=3.0e-7 M⁻¹, K_intra=4.6e-12 — explicitly flagged by the caption as "corresponds best to the experimental data."

Re-ran the sweep with these real, paper-sourced numbers (new function `run_cooperative_valency_sweep_real_units()` in `src/gillespie_nanostar_cooperative.py`). Result: **k=3 still does not win with genuine separation**, but the reason is now mechanistically precise rather than an open "wrong regime" question:

- At `rho_n=1e-8 M` across the paper's own reported σR range (1,000–600,000 μm⁻²), all three tested (K_A, K_intra) pairs from Fig. 4c–e give essentially identical max-alpha values across k=3/6/10 (checked to 4 decimal places).
- Direct inspection of the stationary distribution P(b) confirmed why: across that entire σR range, P(b≥2) is 0.0000 (to four decimals) for every tested k — the population is pinned almost entirely at b=0 or b=1. The first-bond rate (0.65–391, comparable to k_off=1) is fast enough to reach b=1 readily, but the second-bond rate (proportional to σR·K_intra, only 4.6e-9 to 2.76e-6) is 6–9 orders of magnitude below the unbinding rate — so a second bond essentially never forms before the first one breaks, and Delta_G_c (which only ever affects the second bond and beyond) never gets a chance to act.
- A wider exploratory sweep (rho_n from 1e-14 to 1e-2 M, σR up to 1e14 μm⁻², both far outside the paper's own reported experimental ranges) does eventually show P(b≥2) becoming non-negligible and k=3 > k=10 there (0.115 vs 0.028 at σR=1e10, rho_n=1e-6 M) — but by that point Y(σR) is already close to saturating toward 1.0 for all three valencies simultaneously, so the alpha metric (a log-derivative) flattens for all three together rather than showing a real selectivity-optimum separation.

Conclusion: this is not evidence the paper's own finding is wrong. Their k=3 optimum comes from a full joint fit of Delta_G_c together with K_A and K_intra to their own raw intensity data — a joint optimization this project's approach (testing already-published numbers as fixed inputs) cannot recover by construction. Fully closing this gap would require re-fitting all three parameters jointly to the paper's own raw Θ-vs-σR data (available at 4TU Research Data per their Data Availability statement) — judged out of scope for what this project needs, since the methodological point (multivalent capture is nonlinear, cooperativity can shift an optimum) stands on its own without a digit-exact k=3 reproduction.

2 new tests added to `tests/test_gillespie_nanostar_cooperative.py`: one locking in the real-units parameter values against the paper's own reported numbers (so a future edit can't silently drift them), one locking in that the three valencies still tie within 1e-2 of each other at the real reported units.

## 2. Aim 3b real-data write-up

New file: `docs/Aim3b_Cooccurrence_and_Independence_Writeup.md`. A ready-to-paste methods/results section, written directly from `aim3b_cooccurrence_model.py`'s own `__main__` output (re-run this session to confirm the exact numbers, not copied from memory): the fit to Mizenko's real ~10x (EpCAM) / ~5x (Her2) ratios across a 10–30% tetraspanin-rich-fraction range, the robust 1.18x independence-fold-deviation finding (positive correlation between EpCAM and Her2 positivity), and the honestly-ranged 1.02x–1.26x combination-gain sensitivity finding (explicitly NOT presented as a single point value, since it depends on an absolute marker copy-number scale Mizenko's ratio-only data cannot constrain). No new code was needed for this item — it was purely a writing task built on already-tested code.

## 3. f_accessible applied to the CA19-9 real headline finding

Searched PubMed (via a dedicated PubMed MCP tool, not general web search) for real literature on EV protein-corona accessibility effects. Found and read in full: Försönits, Tóth, Jezsoviczky et al. 2025, "Improved Accessibility of Extracellular Vesicle Surface Molecules Upon Partial Removal of the Protein Corona by High Ionic Strength," *Journal of Extracellular Vesicles* 14(7):e70124, DOI 10.1002/jev2.70124, PMID 40660090, PMC12259388.

This paper directly demonstrates the mechanism `f_accessible` represents: a protein corona forms around real EVs (human blood plasma, HEK293T-PalmGFP and THP1 cell lines) and measurably masks surface markers from antibody/annexin V detection; stripping the corona (high-salt washing) significantly improves detection of 15 of 37 tested plasma-EV surface markers (Wilcoxon test, p<0.05 to p<0.001 depending on marker) and 9 of 37 THP1-EV surface markers.

**Critically, this paper does not compare tumor vs. normal EVs** — it compares corona-intact vs. corona-stripped EVs within a single population (healthy donor plasma, or one non-cancer cell line). A separate, targeted PubMed search for a pancreatic-cancer-specific or tumor-vs-normal EV corona/accessibility comparison (search terms: pancreatic ductal adenocarcinoma, EV surface protein, corona, accessibility, tumor, normal, comparison, MUC1/EGFR/TROP2) returned zero results — a genuine, checked negative finding. Additionally, even within the Försönits paper, the actual fold-change magnitudes per marker are shown only in bar-chart figures (Figure 4C/4F), not given as numeric values in the article's own text or tables — confirmed by fetching the full PDF text and grepping for numeric fold-change values, none found outside p-value reporting.

Given this, rather than inventing a plausible-sounding tumor/normal accessibility split with no real source, the new function `ca199_selectivity_sensitivity_to_accessibility()` (in `src/aim3_aim4_range_sweep.py`) runs three explicitly-labeled scenarios against the real digitized CA19-9 data:

1. `no_effect` (f_accessible=1.0 both) — the existing headline finding, ratio range 1.227–1.449.
2. `illustrative_equal_direction` (f_accessible=0.7 both) — confirmed IDENTICAL to `no_effect` at every tested K_INTRA value, verifying `binding_model.py`'s own equal-accessibility-cancels-out proof holds against real digitized data, not just synthetic examples.
3. `illustrative_unequal_direction` (tumor f_accessible=0.49, normal=0.7 — an invented, clearly-labeled 30% additional tumor-specific reduction) — ratio drops to 1.014–1.028: still above 1 (the flip survives in this specific tested scenario) but dramatically weakened from the no-effect baseline.

The function's returned dict includes a `"note"` field stating this honest data-gap finding in plain language, and 3 new tests lock in: the equal-case-matches-baseline guarantee, the unequal-case-weakens-but-doesn't-reverse result, and that the note field actually documents the DOI and the "no comparison found" finding (so a future edit can't silently turn it into generic boilerplate).

## Verification

All 81 tests pass: 76 at the start of this session + 2 new (Gillespie real-units tests, in `tests/test_gillespie_nanostar_cooperative.py`) + 3 new (accessibility-sensitivity tests, in `tests/test_aim3_aim4_range_sweep.py`) = 81. Every module in `src/` (14 total, no new files this session) was run standalone (`python3 <file>.py`) and confirmed to execute cleanly end to end, per the project's standing rule about running actual entry points after touching anything they depend on.
