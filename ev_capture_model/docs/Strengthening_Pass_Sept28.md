# Strengthening Pass — September 28, 2026

Goal: make each expansion angle *stronger*, not just broader. For each angle: the weakest link, what was
done about it, what was found, and what to say (and not say) in the write-up.

Everything below is reproducible: scripts in `ev_capture_model/src/aim13_scoping/`, tests in
`tests/test_strengthening_pass.py` (23 tests), and `verify_documented_numbers.py` now runs 44 checks.

---

## 0. Downloads

| File | Result |
| --- | --- |
| crisprSQL (`100720.zip`, 4.5 MB) | Downloaded through your Chrome; unpacked. 25,632 sites, 87 guides, 17 studies. |
| NAR 2023 supplementary (`gkad702_Supplemental_Files.zip`, 14 MB) | Downloaded through your Chrome (the site's bot-check cleared on its own after a few seconds). Tables S1–S12. |
| DepMap `CRISPRGeneEffect.csv` | **Not obtained.** The portal shows a "verify you're a person" gate. I do not complete those, so this one is yours if you still want it. |

Location: `ev_capture_data/crispr_offtarget/` (SHA-256 hashes recorded in `angle5_precommitment.md`).
The container and your computer's shell are both blocked from these hosts; only your real Chrome could reach them.

DepMap is now optional: the angle that needed it most (Angle 5 essentiality weighting) was dropped for other reasons (below).

---

## 1. Results-document drift — the list is longer than we thought

The verifier now catches **four** stale numbers in `Aim13_Results_Two_Comparator_Ceiling.md`. All are
old-panel values that were never regenerated after the panel grew:

| Section | Document says | Live code says |
| --- | --- | --- |
| §23.4 specificity (safe vs efficacy-fail) | 0.700, p = 0.043, n = 10 | **0.746**, p = 0.011, n = 12 |
| §23.7 reference agreement | ρ = 0.550, n = 42 | **ρ = 0.592**, n = 44 |
| §27.1–27.2 autopsy replication | confinement 0.750, p = 0.046, 9 v 8 (17 targets) | **0.778, p = 0.025, 9 v 9 (18 targets)**; ceiling 0.556 (was 0.528); pooled 0.543 (was 0.514) |
| §27.3 registry control | 0.495, p = 0.530, 20 targets | **0.448, p ≈ 0.68, 24 targets (11 v 13)** |

§27.4's sentence "17 labelled targets… one relabelled target away from crossing 0.05" also needs rewriting for 18 targets.
Direction of every change is neutral-to-favourable; none reverses a conclusion. Fixing them is your rewrite — §23.4's
interpretive prose needs re-deriving, not find-and-replace.

---

## 2. Angle 3 — confidence intervals (the outcome tests)

**Weakest link:** point estimates alone; and two claims that were never given an uncertainty — the replications and the
claim that *expendability adds something beyond confinement*.

**Done:** `outcome_test_robustness.py`.

| Question | Result |
| --- | --- |
| Autopsy replication precision | 0.778, 95% CI [0.519, 0.963] (9 v 9) |
| Registry control precision | 0.448, 95% CI [0.210, 0.699] (11 v 13) |
| Could the specificity control have failed? (power) | 98% power against an effect of AUC 0.88; **69% at 0.75; 33% at 0.65**; min detectable at 80% power ≈ 0.79 |
| Paired AUC gain from adding expendability | +0.033, CI [0.000, +0.100]; carried by three targets (CD22, CD33, CD38) |
| **Direct test: haematopoietic home vs safety** | 8/19 safe vs 0/11 toxic; **Fisher p = 0.013; permutation p = 0.011** |

**What this means.** Two different lenses give different answers about the second factor, and both are true:

* As a *ranking improvement* it is small (+0.033) and bootstrap-fragile — it rests on three targets.
* As an *association* it is real (p ≈ 0.01): every toxic target has a non-haematopoietic home; 8 of 19 safe ones have a haematopoietic one.

**Say:** "Confinement carries most of the ranking power; expendability is a real, statistically supported association
but adds modestly to AUC, driven by CD22/CD33/CD38." **Don't say:** that the two-factor rule is significantly better
than one-factor as a classifier.

**The control is stronger than assumed, with a limit worth stating:** it would almost certainly have caught confinement
tracking "popularity" at the strength it tracks toxicity (98% power at 0.88), so passing means something. It cannot rule out
a modest popularity effect (≤ 0.65).

*Drift caught in my own work:* the first version reconstructed the autopsy group from the wrong key (pooled instead of
within-patient) and read a stale registry JSON. Same failure mode as before; both fixed and pinned in tests.

---

## 3. Angle 1 — continuous weights and the discovery screen

**Weakest links:** (a) hand-assigned weights; (b) the gate claim; (c) the candidate list itself.

### 3a. Weight uncertainty
Perturbing every cell-type weight by ±0.25 (2,000 draws) and re-ranking: **only 9 of the top 25 stay in the top 25 in ≥ 90%
of draws.** The list is stable at the very top and soft below rank ~10.

### 3b. The gate claim was overstated — corrected
The docstring said the ordering is "stable across the upper half" of the gate range. Measured: the known broad-gene artefacts
(AKAP10, TGFBR1) are demoted at *every* gate (ranks ≥ 203 of 214), and the top 25 is identical — but the top 25 is
insensitive by construction, and top-100 overlap falls to ~81–84 of 100 at the ends. Docstring rewritten; **ranks beyond ~50
should not be read as precise.**

### 3c. The candidate list had a real quality problem — found and fixed
Inspecting the stability output showed the "top 10" contained things that are not cell-surface targets: **immunoglobulin
kappa variable genes** (clonal B-cell-receptor fragments), **GAGE-family genes** (multi-copy intracellular
cancer-testis antigens), and secreted proteins. Cause: the screen accepted any gene with an HPA immunofluorescence
"plasma membrane" observation, which fires on any staining pattern, including antibody cross-reactivity.

`candidate_hygiene.py` adds a filter with rules fixed from HPA's own protein-class fields:
hard-exclude (i) immunoglobulin/TCR genes and (ii) genes whose only surface claim is *contradicted* by HPA's topology prediction.

**It was checked, and the first version failed the check.** v1 also excluded any "predicted intracellular" gene without a
CD-marker label — it falsely excluded **27 of 107** known clinically pursued surface targets (MET, MSLN, DLL3, FOLH1, ROR1,
PTK7 …), because HPA routinely gives real membrane proteins both labels. v1 was abandoned and v2 built so positive membrane
evidence is decisive. v2: **6 of 107 (6%) false exclusions, an upper bound** (FN1, KLK2 and NY-ESO-1's gene aren't surface
proteins at all). All 10 positive-control recoveries survive. Revision history is kept in the module docstring.

**Corrected funnel:** 214 → 208 (widened exclusion) → **134 pass hygiene → 128 also not already pursued.**

New stable, genuinely-new, hygiene-passing shortlist (≥ 90% of weight draws in top 25): **ERVFRD-1, PSG1, DCT, HAPLN3, CD1B, GPR50, CD1E.**
Caveats to carry with them: PSG1 carries an HPA "secreted" annotation (kept because it is also a CD-marker gene — flagged, not
excluded); ERVFRD-1's IF reliability is "Uncertain"; DCT is a melanosomal enzyme, so the accessible-epitope question is live.
**This remains a hypothesis-generating screen; it knows nothing about internalisation, shedding, epitope access or dependence.**

*Process note (for honesty in the paper):* the filter was tightened after inspecting results. That is defensible as pipeline
QC only because (i) rules were fixed from independent HPA fields, (ii) the false-exclusion rate was measured on known targets,
(iii) positive controls were re-run, and (iv) the v1 failure is on the record.

---

## 4. Angle 2 — the clinical-units ceiling: diagnosed

**Weakest link:** tiny, unexplained numbers (0.6–1.9% at a 1% harm budget) that didn't rank real targets sensibly.

The diagnosis took three steps, each with a pre-stated hypothesis:

1. **Noise floor** (`clinical_units_noise_floor.py`). The myeloma panel contains FITC and PE control channels. At 1% harm the
   best control scores **1.3%** — above BCMA (1.0%), CD38 (0.6%) and SLAMF7 (0.8%). Across budgets from 1% to 50%, at most
   5 of 61 markers clear both the control and a permutation-null envelope. **On this dataset the metric is at the noise floor.**
2. **Is it a bug?** No. Published per-compartment AUCs are equally low (BCMA 0.61, CD38 0.29). The metric was faithful to a
   dataset in which tumour and normal barely separate.
3. **Rerun where separation is strong** (`clinical_units_dlbcl.py`). DLBCL CyTOF populations reproduce the published
   AUCs to 0.000. Per compartment at a 1% harm budget: **CD19 kills 89.5% of tumour while sparing 99% of T cells** — the metric
   works — but only 1.2% vs "other leukocytes" and 1.0% vs normal B cells.

**Mechanism, tested directly** (`clinical_units_contamination.py`). A budget b sits at the (1−b) quantile of the normal
compartment; if that compartment contains B-lineage cells, the threshold lands inside the tumour's own range and kill
collapses to ≈ b. Estimating B-like contamination *without* the target marker: the "other leukocytes" compartment is
**12%–95% B-like** depending on threshold; removing them lifts CD19's 1%-harm kill from 1.2% to as much as 76%.
That spread *is* the result: **the low-budget ceiling is not identifiable from these compartments** — it measures gating purity.

**Implication for the project.**
* Demote Angle 2 from a headline to a **methods caution**: "tail-based clinical units require normal compartments that are pure
  at the level of the harm budget; ours are not."
* What survives as a clean, legible statement: the **T-cell compartment** (contamination 0–3%): at 1% T-cell harm, CD19 90%,
  CD20 86%, CD22 50%, CD79B 23%, CD40 44%, CD80 9%, HLA-DR 28%. Approved > failed for 3 of 4 approved, with CD79B (23%) below CD40 (44%).
  n = 4 v 3: descriptive only.
* Side effect worth knowing: the "other leukocytes" compartment used in earlier three-compartment work is heavily B-contaminated.
  That only *lowers* ceilings, so the published upper-bound logic stands, but anything read as an absolute "spared" claim against it should be checked.

---

## 5. Angle 4 — the DCT dossier
No new analysis this pass beyond the hygiene results above (DCT passes hygiene, stays top-3, stable in 100% of weight draws).
Still open: a check for any surface-directed DCT programme that has appeared, and DepMap dependence (blocked).

---

## 6. Angle 5 — Cas9 variants: precommitted, tested, and mostly *not* supported as planned

**What the audit changed.** The NAR 2023 data is 306 guides × ≤ 7 *designed* mismatch sites (not genome-wide natural off-targets), so
the original within-guide-curves plan is infeasible. crisprSQL has only 28% of sites with gene IDs. Essentiality weighting dropped.

**The plan was written down first** (`angle5_precommitment.md`, with data hashes): the ceiling theorem's premise says a variant is a
monotone transform of one latent score, so all sites — on and off target — should lie on one curve.

| Test | Pre-stated prediction | Result |
| --- | --- | --- |
| T1: WT → HiFi/LZ3 on-target, one monotone map | R² 0.5–0.85 | **R² 0.07 / 0.02** — wrong; below the 0.40 "premise fails" line |
| T2: same, off-target | tighter than T1 | 0.36 / 0.39 — as predicted (within-class weaker) |
| T3: on-target lies on the off-target curve | positive residual | **+0.11 / +0.055**, CIs exclude 0 — as predicted (variants "move the curve") |

**Post-hoc check (labelled as such).** HiFi and LZ3 agree with *each other* on on-target activity (ρ = 0.87) though each barely tracks WT
(0.25, 0.16). So the T1 failure is not noise: WT appears to **saturate** (all guides ≈ 30%), making it a poor x-axis; T3's positive residual is
also what saturation would produce. The better-designed HiFi → LZ3 comparison: **T1 R² 0.70; off-target R² 0.82; on-target residual −0.005
(CI [−0.010, −0.002])** — on- and off-target sites lie essentially on one curve.

**Decision, by the rule fixed in advance:** both precommitted tests failed the premise, so **Angle 5 is dropped as a validation section.**
The post-hoc result is a hypothesis, not a finding. If you want it back, the honest route is a *new precommitted* HiFi→LZ3 test on an
independent dataset (Tables S11/S12, the T-TISS genome-wide data, are already in hand) — worth doing only if it stays inside one bounded section.
Otherwise a 2–3 sentence note under *scope of the theorem* covers it ("the ceiling bounds a fixed score; engineered variants may change the score").

---

## 7. Summary of what each angle now is

| Angle | Before | Now |
| --- | --- | --- |
| 1 Continuous weights / screen | 214 "new" candidates, ranks taken at face value | 128 hygiene-passing new candidates; top ~10 stable to weight uncertainty, lower ranks soft; gate claim corrected |
| 2 Clinical units | Unexplained tiny numbers | Diagnosed (noise floor + normal-compartment contamination); demote to a methods caution + T-cell-compartment statements |
| 3 CIs | Point estimates | Intervals on every headline incl. replications; control power measured; expendability shown significant as an association, modest as an AUC gain |
| 4 Dossier | DCT pick | Survives hygiene and weight perturbation; open items unchanged |
| 5 Cas9 | Planned section | Precommitted test failed the premise → dropped as validation; post-hoc lead recorded |

## 8. Decisions for you
1. Rewrite §23.4, §23.7, §27.1–27.4 (four verifier failures until you do). Change "214" → "128 new, hygiene-passing" wherever cited.
2. Angle 2: agree to demote to a methods caution?
3. Angle 5: drop (recommended) or spend a bounded effort on the precommitted HiFi→LZ3 confirmation on T-TISS?
4. DepMap: only worth doing if you want tumour-dependence added to the dossier.
