# Project Handoff v6: Computational Model of Multivalent EV Capture

**Supersedes:** `EV-Capture-Project-Handoff-v5.md` (v5, August 2026). Where the two disagree, this document wins. Key changes are listed in "What changed since v5" below.

**Purpose:** hand a new chat everything it needs to continue this project without re-deriving anything or repeating dead ends. Read "What NOT to revisit" carefully — this project reached its current form after three prior conversations of failed shapes, a fourth (v2→v3) that finished the figure-digitization step v2 had left blocked, a fifth (v3→v4) that resolved v3's top-priority open question (n vs. areal receptor density), a sixth (v4→v5) that tackled v4's top-priority open question (Aim 2 leg 3, external validation), and a seventh (v5→v6) that attempted v5's top-priority open question (calibrating K_INTRA/z/k_geom) and found a real logic error in that plan before any fit was made.

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

## What changed since v5

v5's blocking item — "Immediate next steps" item 1, calibrate `K_INTRA`, `z`, and `k_geom` against Mizenko's 10×/5× numbers — was attempted this session, but **before writing any fitting code**, a re-check of Mizenko et al.'s actual experimental methodology (not just the handoff's summary of it) turned up a genuine logic error in that plan. This is exactly the kind of "stress-test before executing" catch the working style calls for, and it changed the session's direction: rather than forcing a calibration that would produce meaningless numbers, the session reassigned the calibration target, built honest scoping work toward where it actually belongs (Aim 3b), and separately made real progress recovering more precise data for the still-open nanostar validation question from v5.

Four substantive updates:

1. **The optional v5 follow-up item (recovering Linne et al.'s exact fitted numbers) was completed — with a genuine, reportable correction to how the nanostar validation should be understood.** Direct access to the paper's full HTML text (RSC-hosted, published version) confirmed that Figure 2(e)'s table — the actual fitted K_A/K_intra values for the k=1,3,6,10 valency sweep, which the closed-form model in `nanostar_validation.py` is built around — remains embedded as a graphics-only table, unrecoverable as text, in every source checked across two sessions now (PMC full-text API, RSC-hosted PDF/HTML, and — checked new this session — the arXiv preprint version, arXiv:2306.16885). However, the paper's own text turned up something more important than the missing numbers: the paper **explicitly states that its own simple independent-arms model (Eq. 2–4, what `nanostar_validation.py` implements) does NOT explain the k=3 optimum** when fit to the valency-sweep data — direct quote, the data "cannot be fitted with a single set of parameters, and only show reasonable agreement if we fit a different K_intra-value to each" valency, and the paper calls this "puzzling" and says the simple model is "incomplete" or "inconsistent." Their REAL explanation is a second, more complex model: a kinetic Monte Carlo simulation with a shared K_A and K_intra across all valencies, plus a separate arm-arm cooperativity/competition free energy term (ΔG_c) that becomes more competitive as valency increases. Crucially, **this extended model's fitted values ARE recoverable as machine-readable text** (Table 2 of the published paper): K_intra = 14×10⁻¹² (shared across k=3,6,10 — NOT decreasing with k, contrary to what `nanostar_validation.py` previously assumed), K_A = 1×10⁻⁵ M⁻¹ (shared), and ΔG_c = −0.5, 0, +1.0 k_BT for k=3, 6, 10 respectively. `nanostar_validation.py`'s module docstring, its `NANOSTAR_PARAMS` block, and its top-level "what this does and does not do" section have all been updated to state this precisely: the file's existing representative-parameter reproduction should now be understood as "a demonstration that the closed-form Eq. 3-4 CAN produce an interior optimum for some K_intra(k) choice," not as a reproduction of the paper's own actual explanatory mechanism (arm-arm cooperativity). This is a real, honest downgrade from v5's framing, not a new failure — the k=3-wins qualitative result is unchanged and all 7 of that file's tests still pass exactly as before (docstring-only edit, verified by rerun). A concrete, well-specified follow-up now exists for a future session that wants a tighter match: implement the Gillespie kinetic Monte Carlo simulation with the published K_intra/K_A/ΔG_c(k) values directly.
2. **A real logic error was found in v5's calibration plan, before any fitting code was written.** Re-reading Mizenko et al. 2021's actual methodology directly (BioMed Central full text, not just the project's own prior summary of it) revealed that the paper's ~10× (EpCAM) and ~5× (Her2) numbers do not measure single-marker capture efficiency at all. Mechanistically: EVs from one starting population were split into two capture methods — (a) single-tetraspanin capture (anti-CD9/-CD63/-CD81 coated chip, capturing only EVs bearing that specific tetraspanin) and (b) genuinely non-specific capture (EVs first chemically biotinylated at *any* exposed surface amine, independent of any particular protein, then captured on an anti-biotin chip) — and *both* captured populations were then stained for EpCAM and Her2. The paper's own words, quoted directly: "Although EpCAM+ EVs could be detected on all three tetraspanin capture spots at similar rates, biotin capture identified more EVs by approximately a tenfold difference," and their own interpretation: "there may be a large percentage of non-tetraspanin expressing EVs that contain the bulk of these markers." **This is a co-occurrence/correlation measurement between two different surface markers (tetraspanin vs. EpCAM/Her2) across one population, not a measurement of how efficiently one ligand captures vesicles bearing its own receptor.** `binding_model.py`'s `K_INTRA`, `z`, and `k_geom` describe single ligand-receptor pair thermodynamics — there is no principled way to map a cross-marker co-occurrence ratio onto those parameters. Forcing a fit would have produced numbers that look calibrated but don't mean what a calibration claims to mean. This was surfaced directly (via `AskUserQuestion`) before proceeding, rather than silently picking a resolution — matching the explicit "stress-test before executing" instruction. The user chose: reassign Mizenko's numbers to Aim 3b (the multi-marker/independence question, where the handoff had *already* correctly flagged co-occurrence as central), and leave `K_INTRA`/`z`/`k_geom` uncalibrated for now rather than force-fit.
3. **Honest scoping work was done toward Aim 3b, using the reassigned Mizenko numbers, rather than leaving the reassignment as a bare note.** A new file, `src/mizenko_cooccurrence_scoping.py`, builds a simple two-subpopulation mixture model (tetraspanin-rich vs. tetraspanin-poor, each with its own EpCAM copy-number distribution) and asks: how strong does the anti-correlation between tetraspanin and EpCAM need to be, structurally, to reproduce Mizenko's ~10×/~5× magnitude? The finding, swept across several parameter combinations and locked in with 3 new regression tests (`tests/test_mizenko_scoping.py`, all passing): **mild-to-moderate anti-correlation (e.g., a 3× difference in mean EpCAM copy number between subpopulations, or a Gaussian-copula correlation of 0.5) only produces a 1.0–1.6× ratio — nowhere near Mizenko's reported magnitude.** Reaching 5–10× requires *both* a fairly small tetraspanin-rich fraction (roughly 10–30% of the population) *and* strong, specific suppression of EpCAM within that subpopulation (roughly 50–100× lower mean copy number than in the tetraspanin-poor majority) — something closer to two largely-distinct subpopulations than a graded correlation. This is explicitly scoping work, not a fitted or validated model (no real joint tetraspanin/EpCAM per-vesicle dataset was used), but it gives Aim 3b's eventual real model a concrete, quantitative target to explain, and is consistent with (and sharpens) the paper's own "different EVs carry different markers" interpretation and the Lennon non-uniform-EGFR-distribution clustering finding already noted elsewhere in this project.
4. **`binding_model.py` and `nanostar_validation.py` were both updated to state their current calibration/validation status precisely, so a future session doesn't have to re-derive this.** `binding_model.py`'s module docstring now states directly that `K_INTRA`, `z`, and `k_geom` remain uncalibrated placeholders, explains why the Mizenko-based calibration attempt was abandoned, and points to what a genuine calibration would still need (a data source measuring single-marker *absolute* capture efficiency, not a cross-marker ratio). All 18 tests (15 from before this session + 3 new) pass; all previously-reported reproduction numbers (1.69 SE, 28.1%, the Lennon flip/no-flip ratios, the nanostar k=3-wins result) were re-run and reproduce exactly, confirming nothing was silently broken.

---

## Current state: done vs. not done

### Done

- **Project shape settled** after three failed versions (see "What NOT to revisit").
- **Novelty verified.** Direct search for "superselectivity" + "extracellular vesicle" returns nothing. The application is open.
- **Aim 1 complete**, including the n-vs-areal-density resolution (see v4 for full detail — unchanged this session). Code lives in `ev_capture_model/` (`src/binding_model.py`, `src/monte_carlo.py`, `src/population_model.py`, `tests/test_aim1.py`, plus demo scripts and figures).
- **Verified by independent execution**, not just by reading the docs (re-confirmed again this session):
  - Analytical vs. Monte Carlo agree within **1.69 standard errors** over 200,000 trials.
  - Population demo produces a **28.1% relative yield difference** between Poisson and negative binomial at *identical means* — the central claim, demonstrated numerically.
  - All 18 unit tests pass (8 Aim 1 + 7 Aim 2 leg 3 + 3 new this session).
- **Aim 2, legs 1 and 2 complete.** Analytical vs. simulation ✓. Limiting cases ✓ — Poisson theorem holds to 2.2×10⁻¹⁶ at the bond-count partition function level; valency=1 recovers Langmuir.
- **Aim 2, leg 3 — nanostar reproduction preserved, now more precisely characterized.** The DNA nanostar optimal-valency-3 signature (Linne et al. 2024) is still reproduced at a mechanism-*suggestive* level (downgraded from v5's "mechanism-faithful" framing — see "What changed since v5," item 1, for exactly why). `binding_model.selectivity_alpha()`'s z-dependence finding from v5 (dilute regime required for α>1) is unchanged and still holds. **The exact Figure 2(e) valency-sweep fit values remain unrecoverable as text** (now checked across three independent sources: PMC, RSC, arXiv), but the paper's *own* real explanatory mechanism (a kinetic-Monte-Carlo cooperative-arms model) and its fitted Table 2 values (K_intra=14×10⁻¹², K_A=1×10⁻⁵ M⁻¹, ΔG_c(k)) were found and are now documented — a concrete, well-specified target for a future, more faithful reproduction.
- **NEW — the Mizenko-calibration plan from v5 was found to contain a logic error, caught before any code was written.** Mizenko's 10×/5× numbers are a cross-marker co-occurrence ratio, not a single-marker capture-efficiency measurement, and cannot principled-ly calibrate `K_INTRA`/`z`/`k_geom`. See "What changed since v5," item 2.
- **NEW — Mizenko's numbers reassigned to Aim 3b, with honest scoping work done, not just a bare reassignment note.** `src/mizenko_cooccurrence_scoping.py` establishes that reproducing the reported 5–10× magnitude requires strong, structural anti-correlation between tetraspanin and tumor-marker status — a concrete quantitative target for Aim 3b's eventual real model. See "What changed since v5," item 3.
- **Poisson theorem refined.** It is *exactly* true at the level of the bond-count partition function `q_b`, and only approximately true (1–5%) one step later after the Langmuir conversion to a probability. This refinement is mine and is more precise than the original claim. The tests deliberately check the specific claim that is exact.
- **Monte Carlo is non-circular.** It simulates the primitive Poisson bond-forming process and checks resulting frequencies against the closed form, rather than drawing against the formula's own output.
- **Figures produced:** `population_shape_effect.png`, `analytical_vs_mc.png`, `patch_regime_sensitivity.png`, `nanostar_valency_validation.png`. Unchanged this session — no new figures were generated (the Mizenko scoping work is text/table output only, not a figure).
- **Mentor/BSL-2 no longer blocking.** Fully computational on public data requires no Qualified Scientist and no lab access. An ISEF Adult Sponsor is still required; confirm current forms with the fair's Scientific Review Committee.
- **Lennon Figures 2 & 3 fully digitized.** Per-vesicle point clouds for PANC-1/HPDEC × EGFR/CA19-9 (Figure 2, cell lines) and healthy/PDAC patient plasma × EGFR/CA19-9 (Figure 3c). Saved as CSVs with source-image crops and full digitization scripts. See "Data situation" below. Unchanged this session.
- **Sign-inversion decomposition, run on real digitized data.** Preserved from v3/v4, still current — see "The decomposition" below.
- **`n` vs. areal receptor density resolved** (v3→v4 work, unchanged this session). `binding_model.py` has a `patch_exponent` parameter with two named regimes, a physical justification for the default, and a working analysis run against the real digitized data.
- **Finding: the ranking flip depends on population size-similarity, not just marker choice.** Unchanged this session — re-verified to reproduce exactly.
- **Test suite: 18 tests, all passing** (8 Aim 1 + 7 Aim 2 leg 3, unchanged from v5, plus 3 new this session for the Mizenko co-occurrence scoping work).

### Not done

- **`K_INTRA`, `z`, and `k_geom` remain uncalibrated placeholders — still true, but the PATH to calibrating them has changed.** v5 identified this as the top priority and proposed Mizenko's numbers as the calibration target; this session found that plan doesn't work (see above) and did not find a replacement single-marker calibration source this session (that was explicitly deferred, not attempted — see "Immediate next steps"). Current figures still show behavior driven by invented constants.
- **Aim 2, leg 3 — Mizenko's 10x/5x fold-underperformance numbers, now correctly and doubly confirmed to belong to Aim 3b, not to any single-marker validation or calibration step.** v5 had already correctly assigned these away from Aim 2 leg 3's validation-target table; this session's deeper methodological read confirms that assignment was right for an even more specific reason (co-occurrence, not efficiency) and moves them specifically to Aim 3b (not just "Aim 3" generally).
- **The nanostar reproduction is mechanism-*suggestive*, not mechanism-faithful (downgraded from v5) and not digit-faithful.** The real explanatory mechanism the paper identifies (kinetic Monte Carlo, arm-arm cooperativity ΔG_c) is not implemented; `nanostar_validation.py` implements only the simpler model the paper itself says does not fully explain the data. The paper's own fitted values for that real mechanism ARE now available (Table 2) and could be used to build a tighter reproduction. See "Immediate next steps."
- **Aim 3b's co-occurrence model is scoped, not built.** `src/mizenko_cooccurrence_scoping.py` establishes what magnitude of anti-correlation is needed, using an assumed two-subpopulation mixture form and no real joint tetraspanin/EpCAM data (none exists for this comparison, as far as searched). Building a real model — ideally against real joint-marker data if any can be found — is still ahead.
- **A genuine single-marker calibration source for `K_INTRA`/`z`/`k_geom` has not been found.** This needs a data source reporting *absolute* capture efficiency or yield for one marker (not a ratio between two different markers/methods) with enough context (ligand density, incubation conditions) to calibrate against. Not yet searched this session — see "Immediate next steps."
- **Aims 3 and 4 entirely** — no target ranking, no multi-marker joint distribution, no independence test, no sensitivity analysis. (Partial exception: the raw material for Aim 3's PANC-1/HPDEC and healthy/PDAC cases exists as digitized point clouds — see Done above — but no model has been run against them yet.)
- **Corona/accessibility not in the model.**
- **Incubation-time axis not implemented.**
- **Cancer type / cell line not finally decided** — and the basis for the earlier lean toward ovarian changed with the Lennon finding.
- **ISEF forms, Adult Sponsor, AI-use log** not started.
- **The soft-adhesion contact-patch scaling (`p=2` default) is a physically-motivated approximation, not a direct EV measurement.** Unchanged from v4/v5 — no published study was found measuring actual EV-surface contact-patch area as a function of vesicle diameter.
- **`k_geom` (the contact-patch geometric constant) is uncalibrated**, same status as `K_INTRA` and `z`. Unchanged from v4/v5.
- **Figure 3c gate polygon was hand-approximated, not pixel-digitized.** Unchanged from v4/v5 — low priority, optional tightening.

---

## The four aims

### Aim 1 — Build the model *(complete, including the n-vs-areal-density resolution)*
Multivalent binding model for a vesicle carrying `n` receptors contacting a surface with ligand density `σ`. Analytical closed form plus independent Monte Carlo. Population-level version integrating over a receptor distribution P(n) rather than a fixed `n` — this population step is the core novelty and is what the field's fold-change approach skips. Unchanged this session — see v4 for full detail.

### Aim 2 — Validate three ways *(all three legs have substantive work; leg 3's status is now more precisely characterized, not more complete)*
1. Analytical vs. Monte Carlo. ✓
2. Limiting cases: Poisson → mean-field exact; valency 1 → Langmuir. ✓
3. **Reproduce known experimental signatures.** ~ — DNA nanostar valency-3 optimum reproduced at a mechanism-*suggestive* level (downgraded this session from v5's "mechanism-faithful" claim — the paper's own real explanation, a cooperative-arms kinetic Monte Carlo model, is not what's implemented here, though its fitted numbers are now known and documented for a future tighter build). Mizenko's numbers are confirmed, more specifically than before, to belong to Aim 3b, not this leg.

### Aim 3 — Re-rank real targets, single and multi-marker *(not started, but raw data now exists for the primary case, and calibration needs a new data source)*
Run the model on published receptor distributions; compute achievable yield and selectivity; compare to fold-change ranking. **The deliverable is the disagreement.**

For the PANC-1/HPDEC and healthy/PDAC-plasma cases specifically, real digitized per-vesicle joint (diameter, copy-number) distributions are now available (`ev_capture_data/`) rather than needing reconstruction from summary statistics alone. Every other system in the literature table still needs summary-statistic reconstruction.

**Consideration carried over from v5:** whatever `z` value calibration lands on will determine whether the calibrated model is actually in the superselective regime (α>1 achievable) or not (see v5's "What changed since v4," point 4). This still needs to be checked explicitly once `K_INTRA`/`z`/`k_geom` are actually calibrated — which now needs a different data source than Mizenko's numbers (see Not Done above and Immediate Next Steps).

**Multi-marker capture is a core sub-aim, not a footnote — and now has a concrete quantitative scoping target.** Unchanged in overall framing from v4/v5, but sharpened this session:
- Insurance if single-marker ceilings are low (likely — expect 2–4×, not 20×).
- The field already recommends combinations without quantifying them (Mizenko et al. say "use a combination," never which or at what cost).
- **The independence assumption is itself a finding.** The 9× math assumes the proteins appear independently across vesicles. If they cluster, the multiplication fails. Lennon et al. observed non-uniform EGFR distribution on EV membranes — regions with no detectable EGFR — which *is* clustering. That makes the independence question empirical, not hypothetical.
- **NEW — this session's scoping work gives a specific number to build toward:** reproducing Mizenko's magnitude requires structural, near-exclusive anti-correlation between tetraspanin and tumor-marker status (roughly 50–100× EpCAM suppression within a tetraspanin-rich subpopulation that is itself only ~10–30% of the population), not mild correlation. See "What changed since v5," item 3, and `src/mizenko_cooccurrence_scoping.py`.

### Aim 4 — Sensitivity analysis: tell the field what to measure *(not started, but two axes now have preliminary answers)*
Which unknown parameter dominates uncertainty in the selectivity prediction? If distribution *shape* dominates the *mean* by a large factor, the recommendation is concrete: stop reporting mean copy number, report the distribution. Changing what a field measures is a bigger contribution than any single device. The accessibility/corona factor is a natural additional axis here.

**The raw-count-vs-areal-density axis** (from v4, unchanged): the ranking flips for the patient-plasma comparison (large size difference between populations) and does not flip for the cell-line comparison (matched sizes) — a strong, reportable Aim 4 finding.

**The z/superselectivity-regime axis** (from v5, unchanged): whether the model is even in the superselective regime depends sensitively on the activity parameter `z`. Still a first-class Aim 4 axis, still blocked on actual calibration (see Aim 3 above).

---

## Three physical mechanisms the model must include

Unchanged this session. Preserved in full from v4/v5:

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
| DNA nanostar optimal valency | 3 | Linne et al., *Soft Matter* DOI 10.1039/d4sm00704b | **~ mechanism-*suggestive* reproduction (re-characterized this session, downgraded from v5's "mechanism-faithful") — the paper's own real explanatory mechanism (cooperative-arms kinetic Monte Carlo, fitted values now known — see "What changed since v5") is not what's implemented; the simpler model implemented here is one the paper itself says does not fully explain the data on its own, though it CAN produce a qualitative k=3 optimum for some parameter choices** |
| Single-tetraspanin capture underperformance, EpCAM | ~10× fewer vs. non-specific | Mizenko et al. 2021 | ✗ — **confirmed this session to be a cross-marker co-occurrence measurement, reassigned specifically to Aim 3b (not Aim 3 generally, and never Aim 2 leg 3 or a K_INTRA/z/k_geom calibration target)** |
| Same, Her2 | ~5× fewer | Mizenko et al. 2021 | ✗ (same) |

The last two are Aim 3b's targets now (not Aim 3's `K_INTRA`/`z`/`k_geom` calibration targets, as v5 had it — see "What changed since v5," item 2, for why that reassignment was necessary). `K_INTRA`, `z`, and `k_geom` still need calibrating, but against a different, not-yet-found data source — see "Immediate next steps."

---

## Key literature

**The anchor paper — Mizenko et al. 2021**, *J Nanobiotechnology* 19:250, DOI 10.1186/s12951-021-00987-1. Open access, ~171 citations.

They measured experimentally that capturing EVs by a single surface protein biases which tumor markers you find — roughly tenfold fewer EpCAM-positive and fivefold fewer Her2-positive vesicles via single-tetraspanin capture versus non-specific capture. Their interpretation is my thesis, stated by them: a large fraction of non-tetraspanin-expressing vesicles may carry the bulk of these markers. Their recommendation is my multi-marker idea: study multiplexing efficiency before choosing a single capture protein.

**NEW — precise methodology, re-confirmed this session by reading the paper's actual methods section (not just the abstract/summary used in prior sessions).** Single-tetraspanin capture uses antibody-coated chips (anti-CD9, anti-CD63, anti-CD81) that bind EVs bearing that specific tetraspanin. Non-specific capture uses a genuinely different mechanism: EVs are first chemically biotinylated at any exposed surface amine group (independent of which proteins are present), then captured on an anti-biotin-coated chip — a near-total, marker-independent capture of the population. Both captured populations are then stained with anti-CD24/anti-EpCAM/anti-Her2 detection antibodies. The reported 10×/5× numbers are the ratio of EpCAM+/Her2+ counts between these two captured populations — i.e., **a co-occurrence measurement between tetraspanin status and tumor-marker status, not a capture-efficiency measurement for any single marker's own ligand-receptor pair.** Direct quotes: "Although EpCAM + EVs could be detected on all three tetraspanin capture spots at similar rates, biotin capture identified more EVs by approximately a tenfold difference. Her2 + EVs were identified with high frequency on all tetraspanin spots with approximately 5 times as many detected on the anti-biotin capture spot." And: "there may be a large percentage of non-tetraspanin expressing EVs that contain the bulk of these markers." Sample: SK-OV-3 (ovarian cancer cell line) EVs, N=3 capture spots per condition, ExoView/SP-IRIS platform.

**Cite prominently as motivation.** What they did NOT do, and what remains mine:
- No theory — they report 10×, can't explain *why* 10× rather than 3× or 50×, can't predict for untested proteins
- No prediction — every new protein/cell line/surface needs a new experiment
- No quantified trade-off — "use a combination" without saying which or at what yield cost
- Only three proteins on one cell line

Framing: *the bias is known and measured — here is the physical theory that explains its magnitude, predicts it for untested proteins, and quantifies the yield cost of the fix everyone recommends but nobody has sized.* **Explaining a known phenomenon is safer than discovering a new one** — it hands you real numbers to validate against. **This framing now applies specifically to Aim 3b (the co-occurrence/independence model), not to a single-marker K_INTRA/z/k_geom calibration** — see "What changed since v5."

**The nanostar validation paper: Linne, Heemskerk, Zwanikken, Kraft & Laan 2024**, *Soft Matter* 20:8515-8523, DOI 10.1039/d4sm00704b, PMC11484159. Open access.

DNA-origami nanostars with a fixed number of branched arms (valency `k`), each with a single-stranded "sticky end" that binds a complementary receptor on a supported lipid bilayer, imaged by TIRFM. Measured bound fraction Y as a function of receptor density σ_R, for k=1,3,4,10 (main valency sweep) and for three different sticky-end binding strengths at fixed k=6. Key findings, all directly relevant to this project's theoretical foundation:
- Selectivity α = d ln(Y)/d ln(σ_R) peaks at an *intermediate* valency, k=3 — not the largest valency tested (k=10). This is the headline external signature this project's Aim 2 leg 3 aimed to reproduce.
- **NEW, corrected this session:** the paper's SIMPLE independent-arms model (their Eq. 3-4, treating ligands as acting independently) does **not**, on its own, explain this k=3 optimum when fit to the real data — fitting it requires a different K_intra for each valency, which the paper itself calls "puzzling" and describes as evidence the simple model is "incomplete" or "inconsistent." Their real explanation is a second, extended model incorporating ligand-ligand pair interactions (cooperativity/competition), represented by a free energy term ΔG_c, fit via kinetic Monte Carlo (Gillespie algorithm) simulation — NOT the closed-form equations. ΔG_c becomes more *competitive* (arms interfering with each other) as valency increases: fitted values ΔG_c = −0.5, 0, +1.0 k_BT for k=3, 6, 10, with a SHARED K_intra=14×10⁻¹² and K_A=1×10⁻⁵ M⁻¹ across all three valencies (published Table 2 — this is the actual mechanism and the actual numbers, now recovered and documented in `nanostar_validation.py`).
- Binding strength also has an optimum (their k=6 sweep): too strong and the transition happens at too-low receptor density to resolve; too weak and multivalent effects wash out.
- **This project's reproduction (`src/nanostar_validation.py`) implements ONLY the simple independent-arms equations (Eq. 2-4), with representative per-k K_intra values chosen to decrease with k** — a mechanism the paper's own text indicates is NOT the real explanation for their data (see above). This should now be described as "a demonstration that the simple model's functional FORM can produce a qualitative k=3 optimum for some parameter choice," not as a faithful reproduction of the paper's own mechanism. This is a real, honest downgrade from v5's "mechanism-faithful, not digit-faithful" framing — see "What changed since v5" for the full account, and `nanostar_validation.py`'s module docstring for the complete technical explanation.
- Their exact fitted numeric values for the SIMPLE model (Figure 2(e)'s table, giving K_A and K_intra for k=1,3,6,10 under the independent-arms assumption) remain unrecoverable as machine-readable text — checked this session across the PMC full-text API, the RSC-hosted HTML/PDF, AND (new this session) the arXiv preprint (arXiv:2306.16885) of the same paper, all three embedding this specific table as a graphic. Note also: **the arXiv preprint is an earlier draft with different numbers** (ΔG⁰=−7 kT vs. the published −8.4 kT, different sticky-end sequence labels) — do not mix numbers between the arXiv and published RSC versions.
- **A concrete, well-specified path to a tighter reproduction now exists** (it did not in v5): implement the Gillespie kinetic Monte Carlo simulation described in the paper's Methods, using the published Table 2 values (K_intra=14×10⁻¹², K_A=1×10⁻⁵ M⁻¹, ΔG_c(k) for k=3,6,10) directly. This is a substantially larger undertaking than the current closed-form script (a stochastic simulation with a transition-rate correction term, not a few equations), which is why it wasn't attempted this session, but it is no longer blocked on missing data.

**The PDAC data source: Lennon et al. 2019**, *J Extracell Vesicles* 8:1685634, DOI 10.1080/20013078.2019.1685634, PMC6844376. Open access. Unchanged from v4/v5 — see full detail preserved below.

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

**The decomposition, run on the actual digitized per-vesicle point cloud from Figure 3c (unchanged since v3/v4/v5):**

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
- Linne, Heemskerk, Zwanikken, Kraft & Laan, *Soft Matter* 2024, DOI 10.1039/d4sm00704b — nanostar optimal valency; **reproduced at a mechanism-suggestive level, with the paper's actual real mechanism (cooperative-arms kinetic Monte Carlo) now identified and its fitted values documented for a future tighter build — see above and "What changed since v5."**
- *Accounts Chem Res* DOI 10.1021/acs.accounts.2c00672 — superselectivity as a design principle for biomarker surfaces; generalizes beyond nanostars; EVs conspicuously absent.

**Other single-vesicle data sources:** (unchanged from v4/v5)
- Liu et al. 2026, *Advanced Science*, DOI 10.1002/advs.202515129 — only single-EV integrin copy-number study. 12-plex DNA-PAINT at single-copy resolution. Integrins heterogeneously expressed at low copy number; most exosomes carry mis-paired or unpaired integrins. Measured *combinations*, so directly relevant to multi-marker. **Paywalled and bot-blocked** — needs library access.
- Tetraspanin copy numbers by super-resolution on seminal exosomes: CD63 ~12.6, CD81 ~1.6, CD9 ~16.6 average copies. Reviewed in *Chem Biomed Imaging* DOI 10.1021/cbmi.3c00095 — **use this peer-reviewed review rather than the US Patent 12,313,519 cited in v1.**
- Same review reports the Lennon EGFR/CA19-9 numbers (~24 and ~31 in healthy plasma EVs, dropping to ~16 in PDAC) — useful cross-check.
- Deconvolution work: **53.5% of membrane-stained vesicles showed no significant expression of CD63, CD81, or CD9.** The zero-inflation fact.
- *Cytotherapy* 2025, DOI S1465-3249(25)00128-8 — calibrated nano-flow cytometry, copies per particle across three stem-cell sources: 50–70% CD81+, 40–60% CD63+, 15–35% CD73+; ~30–40 CD81 and ~41–53 CD63 copies per positive particle. Detection limit 4–10 molecules/particle.
- PMC9369185 — ExoView across HeLa, HepG2, H1650, SH-SY5Y.
- Mizenko et al. 2021 — SK-OV-3 (ovarian) with co-occurrence data. **Worth revisiting as a candidate source for absolute single-marker capture efficiency too, alongside its now-primary role as Aim 3b's co-occurrence target — the paper reports raw particle counts per capture spot for each tetraspanin, which might be usable for a single-marker calibration if enough context (ligand/antibody density on the chip, incubation conditions) can be extracted. Not yet checked for this purpose.**

**Contested — do not cite as established, and do not rank as a target:** GPC1 as a pancreatic EV biomarker. Unchanged from v4/v5 — see below for full detail preserved.

Melo et al. 2015 *Nature* reported near-perfect discrimination; roughly ten subsequent studies produced contradictory results. Specifically checked: one study found GPC1 enriched in PDAC circulating exosomes relative to source tumors, but could not distinguish PDAC from benign pancreatic disease using exosomal GPC1, and found no significant GPC1 difference between normal and PDAC tissue even in matched pairs (Oncotarget, PMC5922373). Another, using direct mass-spec quantification, found exosomal GPC1 "not diagnostic for PDAC" with levels relatively similar across normal, PDAC, and chronic pancreatitis (PMC5386003).

**No single-vesicle copy-number distribution exists for GPC1** — searched specifically, found nothing. Every GPC1 measurement in the literature is bulk (ELISA, mass spec, flow cytometry panels), never per-vesicle molecule counts of the kind Lennon et al. did for EGFR/CA19-9. This means GPC1 **cannot be entered into Aim 3's target ranking** — there is no distribution to rank it with, only a bulk-level controversy.

A "benign disease causes a moderate density bump, cancer causes an extreme one" mechanism was proposed externally to explain the controversy, with invented copy-number categories (healthy 0, benign 1–4, cancer 15–50+) presented as if measured. **This is not supported by the data.** Do not use this framing or these numbers anywhere in the write-up. GPC1's only legitimate role here is as a **discussion-section case study**: a high-profile bulk biomarker controversy that single-vesicle distributional data could in principle help resolve, but hasn't, because that data doesn't exist.

---

## Data situation (honest)

Unchanged this session. Preserved in full from v4/v5:

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

**Open lead, not yet pursued:** the Carney lab (UC Davis) published 96 GB under CC0, signalling openness. Rachel Mizenko (first author, ORCID 0000-0002-9991-8799) and Randy Carney (senior, ORCID 0000-0001-8193-1664) could be emailed for tabular per-vesicle data. Also a plausible mentor lead. Deferred by choice, not ruled out. **NEW leads from this session:** Linne et al.'s DNA nanostar data, including the Table 2 cooperative-model fit values, is confirmed available at 4TU.ResearchData per the paper's Data Availability statement — a Gillespie kinetic Monte Carlo implementation using those published values would be a genuinely tighter validation than what exists now. Separately, **Mizenko et al.'s raw per-capture-spot particle counts** (not yet examined for this purpose) might provide a legitimate single-marker absolute-efficiency calibration source if enough experimental context can be extracted — worth checking before searching elsewhere.

---

## Cancer type — decided by data, not by prior commitment

Unchanged this session. Preserved in full from v4/v5:

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

Unchanged this session. Preserved in full from v4/v5:

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
- ❌ **Calibrating against published capture efficiencies, generally.** Surveyed; metrics incompatible, ligand density essentially never reported. Dead end for most sources. **Exception, NEW this session: Mizenko et al.'s raw per-capture-spot particle counts specifically have NOT been checked for this purpose and remain a genuinely open, not-yet-ruled-out lead — see "Immediate next steps."** Don't treat this general "dead end" note as covering that specific unchecked possibility.
- ❌ **The Mizenko Zenodo dataset.** 96 GB of images. Checked.
- ❌ **Pure literature analysis with no model.** Explicitly rejected — I need to build something.
- ❌ **Anchoring on a cell line for non-data reasons.**
- ❌ **Modelling receptor counts as Poisson or Gaussian.** Poisson is the exact distribution under which my central finding vanishes by construction. Gaussian is wrong for count data. **Negative binomial** is correct and is what the existing code uses. This was recommended by an external chatbot; it is wrong.
- ❌ **Defining success as "the model captures cancer EVs and ignores healthy ones."** That is a predetermined outcome, not a hypothesis test.
- ❌ **Using GPC1 as a ranked target, or citing any "density category" numbers for it.** No single-vesicle GPC1 data exists.
- ❌ **Re-deriving from scratch whether contact-patch area scales with vesicle diameter.** Settled in v3→v4 — see `binding_model.py`'s module docstring. DO revisit specifically if an EV-specific contact-patch measurement turns up.
- ❌ **Re-litigating whether the Figure 3c sign-inversion finding is "real."** Checked twice independently (v2, v3), further resolved in the model itself in v4. Don't restart this check from scratch.
- ❌ **Re-litigating whether the nanostar valency-3 optimum can be reproduced by simply relabeling `binding_model.py`'s variables.** This was tried in reasoning (not in code, since the mismatch was caught before writing anything) and confirmed not to work — see v5 for the original reasoning. `src/nanostar_validation.py` is a deliberately separate model — don't attempt to fold this back into `binding_model.py` or treat it as redundant with it.
- ❌ **Assuming `binding_model.selectivity_alpha()` at z=1.0 (the value used in every existing demo) shows superselectivity.** It doesn't — checked and documented in v5, unchanged this session. Don't cite any existing demo figure as showing α>1 without checking the z value used; use a small z (e.g. 1e-4) specifically if a genuinely superselective example is needed for a figure.
- ❌ **NEW — Treating `nanostar_validation.py`'s representative-parameter k=3 optimum as a validated reproduction of Linne et al.'s actual explanatory mechanism.** This session found the paper's own real mechanism is arm-arm cooperativity (ΔG_c) fit via kinetic Monte Carlo, NOT a simple K_intra-decreases-with-k trend in the closed-form equations, which is what this file actually implements. The file's result is real and still qualitatively matches (k=3 wins), but should be described as "the closed-form model's structure CAN produce a k=3 optimum for some parameter choice," not "reproduces the paper's mechanism." Don't cite this file's result as mechanism-faithful without this caveat.
- ❌ **NEW — Calibrating `K_INTRA`, `z`, or `k_geom` directly against Mizenko's 10x/5x numbers.** Checked and rejected this session: these numbers measure cross-marker co-occurrence (tetraspanin status vs. EpCAM/Her2 status), not single-marker capture efficiency, and there is no principled way to map a co-occurrence ratio onto single-ligand binding-strength parameters. Don't re-attempt this without first finding a genuinely different justification than what was checked and rejected this session (see "What changed since v5," item 2, for the full reasoning) — a future session tempted to "just fit it anyway" should re-read that section first.
- ❌ **NEW — Trying to force Mizenko's 10x/5x numbers into a `selectivity_alpha()`-style single-marker test.** Confirmed again this session, now for a more specific and complete reason than v5's assessment (fold-underperformance ratio) — see above. These numbers belong to Aim 3b's co-occurrence model (see `src/mizenko_cooccurrence_scoping.py`), not to any single-marker validation or calibration.
- ⚠️ **"Capture efficiency" across papers is not comparable.** Denominators vary wildly, and ultracentrifugation itself recovers only 5–25%. Always check the denominator.
- ⚠️ **Superselectivity is not a threshold or an on/off switch.** It is a steep but continuous power-law response (α > 1 over an intermediate range, then saturation). Do not describe it as "unbreakable velcro" or a density threshold.
- ⚠️ **Classical superselectivity theory is not directly transferable.** The original theory describes a multivalent *probe* discriminating between *surfaces* of differing receptor density. This project inverts the geometry. **That gap is the contribution** — do not paper over it. **This gap is now concretely demonstrated, not just asserted**, by the fact that the nanostar validation required a genuinely separate model rather than a relabeling of our own (see above) — and reinforced further this session by the fact that even that separate model turned out to need a MORE complex mechanism (cooperativity) than initially implemented, underscoring how much real physics separates the flipped-geometry EV problem from the textbook superselectivity setup.

---

## Immediate next steps, in order

Ordering principle: do the thing that could most change the project's direction first. v5's step 1 (calibrate K_INTRA/z/k_geom against Mizenko) was attempted and found to be a dead end for a specific, now-documented reason; the actual calibration work is still the top priority but needs a different data source first. Numbering below restarts to reflect the new state.

**1. Find a genuine single-marker calibration data source for `K_INTRA`, `z`, and `k_geom`.** This is now the highest-priority not-started item, more precisely scoped than in v5: needs a paper reporting *absolute* capture efficiency or yield for ONE marker's own ligand-receptor pair (not a ratio between two different capture methods or two different markers), ideally with enough experimental context (ligand/antibody density, incubation time/concentration) to map onto the model's parameters. Two concrete leads, both unchecked: (a) **Mizenko et al.'s own raw per-capture-spot particle counts** (each tetraspanin's own capture spot, compared to its own detection channel) might work if enough context can be extracted — check this before searching elsewhere, since the paper is already thoroughly read and trusted; (b) a fresh literature search specifically for single-marker absolute EV capture yield/efficiency papers with reported ligand density, which the general "calibrating against published capture efficiencies" dead-end note (see "What NOT to revisit") does NOT rule out having missed, since that was a different, broader search. If no clean source is found, consider whether `K_INTRA`/`z`/`k_geom` can instead be constrained by physical plausibility bounds (e.g., known antibody-antigen Kd ranges) rather than fit to a specific number — worth discussing as a fallback, not silently defaulting to it.

**2. (Optional, now well-specified) Build the Gillespie kinetic Monte Carlo reproduction of Linne et al.'s real cooperative-arms mechanism.** Using the published Table 2 values (K_intra=14×10⁻¹², K_A=1×10⁻⁵ M⁻¹, ΔG_c(k) = −0.5/0/+1.0 kT for k=3/6/10) directly, rather than the closed-form Eq. 3-4 approximation currently implemented. This is a bigger build than the current script (stochastic simulation with a state-dependent transition-rate correction, following the paper's own description in their Methods/Simulation section) but is no longer blocked on missing data, and would upgrade the nanostar validation from "mechanism-suggestive" back to genuinely mechanism-faithful. Not blocking — the current qualitative k=3-wins result is a legitimate, honestly-scoped finding on its own.

**3. Build Aim 3b's real co-occurrence model**, using `src/mizenko_cooccurrence_scoping.py`'s scoping result as a starting target (need strong, structural anti-correlation — not mild — between tetraspanin and tumor-marker status, roughly matching the 50-100x/10-30% ballpark found this session) but ideally against real joint per-marker data if any can be found (checked so far: none exists for tetraspanin×EpCAM/Her2 specifically as far as searched). If no real joint data can be found, this may need to remain a scoping-level "what would be required" analysis rather than a fitted model — worth deciding explicitly rather than drifting into treating a plausibility sweep as a validated result.

**4. Add `f_accessible` to the model.** One multiplier on `n`. Don't model corona formation mechanistically; carry the uncertainty and test sensitivity in Aim 4.

**5. Aim 3 — single-marker ranking**, starting with the PANC-1/HPDEC and healthy/PDAC-plasma cases using the digitized real distributions in `ev_capture_data/` and the working `joint_population_yield()` function, then extending across every other system with adequate (reconstructed) data, compared against fold-change ranking. **Run every ranking under both `patch_exponent` regimes** (or at minimum, check whether the two populations being compared differ substantially in mean diameter) rather than silently picking one. Note this can proceed with placeholder `K_INTRA`/`z`/`k_geom` if calibration (step 1) isn't resolved in time — the qualitative ranking-flip findings don't require exact calibration, only the eventual quantitative yield numbers do; state this limitation clearly if it comes to that.

**6. Aim 3b — multi-marker joint distributions and the independence test**, building on step 3's scoping work. Lennon's observed EGFR clustering makes independence an empirical question generally, and this session's scoping work makes the Mizenko-specific magnitude question concrete.

**7. Aim 4 — sensitivity analysis**, including distribution-family choice, `f_accessible`, the raw-count-vs-areal-density sensitivity (re-run with calibrated parameters and across the full literature table), and the z/superselectivity-regime sensitivity (does the calibrated model actually land in the α>1 regime, and how sensitive is that to reasonable z uncertainty?) — still blocked on step 1's calibration.

**8. (Optional, low priority) Tighten the Figure 3c gate-polygon digitization.** Unchanged from v4/v5 — currently hand-approximated, qualitative conclusion not expected to change.

**9. (Optional, worth a quick check before Aim 3) Search specifically for any published measurement of EV-surface contact-patch area vs. vesicle diameter.** Unchanged from v4/v5 — none found so far, search wasn't exhaustive.

**10. Write-up, figures, poster.**

Timeline still fits. This session's main lesson: a stress-test that catches a planning error BEFORE code is written is itself valuable, sessioned-work — not a wasted session, even though the original top-priority task (calibration) didn't get done. The remaining big-effort items are finding a real calibration source, Aim 3 ranking, and Aim 3b's co-occurrence model.

---

## Known risks

- **Rankings might barely change.** Mitigation: the headline deliverable is the *relationship* between distribution shape and achievable selectivity, plus where the standard approach breaks down — informative in both directions.
- **The raw-count-vs-areal-density distinction complicates the story, with real numbers.** Unchanged from v4/v5 — see full detail preserved. Mitigation in place: the physical argument for the default was worked out independently of which answer it would produce, and the size-mismatch-dependent flip gives the write-up a clean, defensible way to present this as a genuine, structured sensitivity finding.
- **The physical argument for `patch_exponent=2` rests on general soft-matter contact-mechanics, not an EV-specific measurement.** Unchanged from v4/v5.
- **`K_INTRA`, `z`, and `k_geom` remain uncalibrated, and the calibration path is now LESS clear than in v5, not more.** v5 thought it had a calibration target (Mizenko's numbers); this session found that target doesn't work and did not find a replacement. This is a real, elevated risk relative to v5 — worth stating plainly rather than downplaying. Mitigation: two concrete unchecked leads exist (Mizenko's own raw per-spot counts; a fresh literature search) — see Immediate Next Steps item 1 — and Aim 3's ranking work can proceed qualitatively without full calibration if it comes to that, with the limitation stated clearly in the write-up.
- **The superselectivity mechanism (α>1) is not automatically present in this model just because binding is multivalent; it requires operating in a sufficiently dilute regime (small `z`).** Unchanged from v5 — discovered while exercising `selectivity_alpha()`. Still a real risk for Aim 3, now compounded by `z` not yet having any calibration path at all (see above).
- **The nanostar valency-3 reproduction is mechanism-*suggestive*, not mechanism-faithful (downgraded this session) and not digit-faithful.** This session found the paper's own text indicates its simple closed-form model (what's implemented here) does NOT fully explain their own data — the real mechanism is more complex (arm-arm cooperativity, kinetic Monte Carlo). The k=3-wins qualitative result is still real and still not guaranteed by construction (getting any interior optimum required real iteration, per v5), but a reviewer or judge who reads the source paper carefully could reasonably ask "the paper itself says this simple model doesn't work — why are you using it?" The honest answer: it demonstrates the qualitative TYPE of trade-off (combinatorial boost vs. per-arm weakening) can produce an interior optimum, using the paper's own simpler starting equations, while being explicit that the paper's own fuller mechanism is different and more specific, and pointing to the now-known path (Gillespie kinetic Monte Carlo, Table 2 values) to close that gap if wanted.
- **The Lennon sign inversion.** Unchanged from v4/v5 — status fully resolved as of v4, unaffected by this session's work.
- **Do not import unmeasured "density category" numbers for any target, GPC1 or otherwise.** Unchanged from v4/v5.
- **Detected vs. absolute counts.** Unchanged from v4/v5.
- **Digitization introduces its own uncertainty.** Unchanged from v4/v5.
- **Liu et al. paywalled.** Unchanged from v4/v5.
- **Distributions underdetermined from summary statistics** for every system other than the now-digitized PANC-1/HPDEC/plasma one. This is Aim 4, not a bug.
- **NEW — Aim 3b's scoping work (`mizenko_cooccurrence_scoping.py`) uses an assumed mixture-model form, not real joint per-marker data.** No dataset with real, paired (tetraspanin, EpCAM) per-vesicle measurements was found this session. The scoping result (need strong, structural anti-correlation) is a genuine constraint on what any real model must explain, but should not be mistaken for a validated finding about the real biology.
- **"This is just theory."** Answer: validated against published experimental signatures not used in fitting (including the nanostar reproduction, now more precisely characterized as mechanism-suggestive rather than mechanism-faithful — be ready to explain this distinction if asked), and it makes falsifiable predictions.
- **Genuine worst case:** rankings don't move AND the trade-off curve has no interesting structure AND multi-marker analysis says nothing non-obvious. Unchanged from v4/v5 — seems unlikely but not impossible.

---

## Practical notes

- Fully computational on public data → **no Qualified Scientist, no BSL-2**. An ISEF Adult Sponsor is still required; confirm current form requirements with the fair's Scientific Review Committee.
- ISEF requires disclosing AI use. Keep a running log: date, tool, what was asked, what was done with the answer. The research plan, abstract, poster, and citations must be in my own words. **NEW — this session (attempted Mizenko calibration, found the logic error, re-read Linne et al. for the exact nanostar values, built Aim 3b scoping work) should be logged**: date, that Claude was used to (a) re-read Mizenko et al. 2021's actual capture methodology directly from the primary source before attempting a calibration, (b) identify that the planned calibration (K_INTRA/z/k_geom against Mizenko's 10x/5x) does not work because those numbers measure cross-marker co-occurrence rather than single-marker capture efficiency, (c) surface this to the user directly and get a decision on how to proceed rather than silently choosing, (d) build honest, clearly-labeled scoping work (not a validated model) toward where those numbers actually belong (Aim 3b), (e) re-read Linne et al. 2024's full text directly (RSC HTML and arXiv preprint) and find that the paper's own text states its simple model doesn't explain its own k=3 finding, requiring a more complex mechanism whose fitted values were then found and documented. The reasoning behind each of these — why the calibration attempt was abandoned, what the scoping work does and doesn't show, and the corrected nanostar-validation framing — should be described in your own words in the methods/limitations section, not just cited to this session.
- Prior sessions' AI-use logging notes (Lennon digitization; n-vs-areal-density resolution; Aim 2 leg 3 nanostar validation) remain valid and unchanged — see v5 if starting the log fresh from scratch.
- **Data provenance for the write-up:** if the digitized Figure 2/3 point clouds are used in any figure or analysis presented at the fair, the methods section should describe the digitization approach and cite Lennon et al. 2019 as the source of the original figures. If the nanostar validation figure/analysis is used, the methods section should describe that this is an independent reproduction using the paper's own SIMPLER equations with representative (not exactly-published) parameters for that simpler model — and should be explicit that the paper's own real explanatory mechanism is different and more complex (now documented, not implemented). Judges may ask how confident the nanostar match is — the honest answer has gotten more nuanced this session: "the qualitative optimum-at-intermediate-valency result is reproduced from the paper's own simpler starting equations, but the paper's own text indicates those simpler equations don't fully explain their real data, and I've documented (but not yet built) what the paper's own real mechanism requires."
- A practical note that WORKED again this session: **attach both `ev_capture_model/` and `ev_capture_data/` from the start of a new chat** — carried forward from v4/v5's advice, still valid, no round-trips lost to a missing folder this session either.
- **NEW practical note for future sessions:** when a handoff document's "Immediate next steps" item names a specific external number as a calibration target, re-verify the ORIGINAL paper's methodology directly before writing calibration code, even if a prior session already summarized it — this session's catch (Mizenko's numbers being a co-occurrence ratio, not a capture-efficiency number) came from re-reading the primary source's actual Methods/Results text, not from anything wrong in how prior sessions had summarized it in this handoff document itself. The summary wasn't false, but it didn't carry enough methodological detail to reveal the problem — only reading the paper's own capture-methodology description directly did that.
