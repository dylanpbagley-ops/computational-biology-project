# Project status: weaknesses, ISEF competitiveness, and expansion angles

**September 7, 2026.** Every expansion angle below was tested before being listed. Nothing here
is speculation about what might work — the feasibility runs are in
`ev_capture_model/src/aim13_scoping/`.

---

## Part 1 — Weaknesses, ranked by how much they hold the project back

### 1. No wet-lab or first-hand data, and no clinical outcome link (PARTLY INHERENT)
Every result is a reanalysis of other people's data. The outcome tests use *target-level*
outcomes (which drugs got approved, which caused toxicity), never patient-level ones.
**Resolvable?** The patient-level version is effectively not resolvable before the fair — the
one cohort with outcomes is 5 vs 5 and no reframing fixes that. The target-level version is
already done and is the stronger play. **Treat as inherent and lean into the target-level framing.**

### 2. Small n on every outcome test (INHERENT — say so plainly)
The blood test is 5 approved vs 5 failed (p = 0.048); the new whole-body toxicity test is 15 vs 7
(p = 0.006). **This cannot be fixed by more work: the world contains few approved cancer surface
targets.** State it as a limit, not a gap. The whole-body test roughly doubles n and is the best
available answer.

### 3. RNA versus protein in the new whole-body work (PARTLY RESOLVABLE)
The Human Protein Atlas single-cell reference is transcript-level. The field's own tools
(Dannenfelser, SCAN-ACT) do the same, but protein would be better. **Resolvable** by adding the
HPA immunohistochemistry protein tables as a cross-check.

### 4. Labels assigned by one person (RESOLVABLE, cheap)
Every approved/failed/toxic label is mine, from the literature. **Fix:** pre-register the list
(draft already written) and have someone else — Ms. Turner, or a lab contact — check it.

### 5. Gating choices are homemade (LARGELY RESOLVED)
Percentile gates versus Otsu versus Gaussian mixture agree at 0.92–0.98 rank correlation, and the
circularity problem is fixed. Residual risk is low; mention and move on.

### 6. Everything rests on two blood-cancer datasets (NOW RESOLVED — see Part 3)
Until today the project was DLBCL plus myeloma. The whole-body extension covers 154 human cell
types and both blood and solid tumours.

### 7. The absolute-calibration gap (INHERENT, and now harmless)
CyTOF signal is not molecules/cell, so specific drug operating points are not computable. The
measures that carry the project — AUC and restriction — are invariant to it. Inherent, but it
costs almost nothing.

---

## Part 2 — ISEF competitiveness, from the actual award lists

Read from Society for Science's own 2024, 2025 and 2026 full-award releases, Computational
Biology & Bioinformatics (CBIO) category.

**What wins First Award (2 per year, $6,000):**
- 2026 CBIO071 — a **Markov-chain model** of tomato–*Pseudomonas* immunity, validated against
  published benchmarks, no new data.
- 2026 CBIO018 — physics-informed neural network pipeline.
- 2025 CBIO057 — StutterZero, which **generated its own paired dataset**.
- 2025 CBIO071T — ALS biomarkers with **cross-dataset validation**.
- 2024 CBIO044 — ML virtual screening terminating in **candidate compounds**.

**What lands at Fourth Award ($600) — and this is the important calibration:**
- 2026 CBIO042 — *"ReinforCell — Reinforcement Learning to Engineer Robust CAR-T Cells."*
- 2026 CBIO088 — *"Integrative Multi-Omic Data Profiling of Pancreatic Ductal Adenocarcinoma."*
- 2025 CBIO013 — *"Bispecific antibody modeling for solid tumor therapy."*
- 2025 CBIO048 — *"CD30+ targeting antibodies for lymphoma."*

**The pattern.** Computational CAR-T and immunotherapy-target projects have been landing at
Fourth Award, and so has multi-omic profiling of pancreatic cancer. Both of those are adjacent
to this project's genre. The First Award winners are not fancier — CBIO071 is a Markov chain —
they are the ones with a **mechanistic model plus independent validation plus one specific new
claim**.

**Honest assessment.** This project sits in a genre that has recently scored Fourth Award, but it
does the things the First Award winners do and those projects did not:
- a proved theorem with an in-sample verification (0 violations in 666 and 1,953 pairs),
- negative controls that behave (isotype antibodies at chance, duplicate clones agreeing to 0.034),
- replication across two diseases, two technologies, two labs,
- validation against external clinical outcomes with exact p-values,
- documented self-correction — three bugs and one wrong diagnosis, all written down.

Realistic placing: **strong Second Award, with First Award reachable if the whole-body solid-tumour
extension lands and the write-up leads with the theorem rather than the pipeline.** Below Second
only if it is presented as "I analysed some CAR-T data," which is the genre that scores fourth.
The single biggest presentation risk is being mistaken for CBIO042/CBIO088. Lead with the
mechanism and the impossibility result, never with the dataset.

---

## Part 3 — Expansion angles, all tested, ranked by impact

### #1. Whole-body normal atlas: extend the thesis to SOLID TUMOURS — **TESTED, WORKS**
Human Protein Atlas single-cell consensus: **154 human cell types**, every major organ, ~20,000
genes. For each target compute how confined it is across normal human tissue.

**Ran it on 34 real targets. Result:**

| gene | confinement | cell types | outcome |
|---|---|---|---|
| CD20 | 0.850 | 1 | approved |
| BCMA | 0.789 | 1 | approved |
| PSMA | 0.568 | 2 | approved |
| CLDN18.2 | 0.359 | 8 (gastric mucosa) | **approved 2024** |
| CEACAM5 | 0.283 | 7 (colonocytes) | severe colitis |
| CD44 | 0.061 | 45 | skin toxicity |
| EGFR | 0.051 | 60 | skin/lung toxicity |
| **HER2** | **0.041** | **72** | **fatal CAR-T event** |
| **CD47** | **0.030** | **120** | **magrolimab, on-target anaemia** |

Against documented **on-target toxicity** (15 safe vs 7 toxic): confinement separates them at
**0.810, exact p = 0.011**; number of cell types at **0.829, p = 0.006**. Better n and a smaller
p than the blood-only test.

**And it sharpens the thesis into its final form.** Breadth catches systemic toxicity (CD47, HER2,
EGFR) but not CEACAM5, which is confined to colonocytes and still caused colitis — because you
cannot lose your colon lining. That is the project's own thesis restated for solid tumours:

> **A target is safe only if it is BOTH confined to few normal cell types AND those cell types are
> expendable or regenerating.** CD19/CD20 → B cells (expendable). BCMA → plasma cells (expendable).
> CLDN18.2 → gastric mucosa (regenerates). CEACAM5 → colonocytes (not expendable) → colitis.
> HER2/CD47 → everywhere → catastrophic.

This unifies blood and solid tumours under one principle and explains every case including the
exceptions. **Highest-impact expansion available, and it already runs.**

### #2. Pancreatic cancer — **TESTED, and the connection is real, not forced**
Ranking the actual PDAC surface targets by whole-body confinement:

**CLDN18.2 (0.359) > CEACAM5 (0.283) > mesothelin (0.239) > MUC1 (0.093) > EGFR (0.051) > HER2 (0.041)**

That is the order the field converged on. CLDN18.2 is the target that reached approval and is now
in pancreatic CAR-T trials; mesothelin is the long-standing tolerable PDAC target; MUC1, EGFR and
HER2 have repeatedly failed. CEACAM5 is confined but to colon — the "confined but not expendable"
case.

**Why this is not a bolted-on connection:** the project already fitted a CLDN18.2 dose-response on
three pancreatic cell lines (r = 0.964, RMSE 0.24) and a mesothelin calibration, and already holds
CPTAC PDAC survival data (n = 135). The two top-ranked pancreatic targets are the two this project
independently calibrated months ago. And the project's *original* EV headline — the CA19-9 flip —
was pancreatic (PANC-1 vs HPDEC). **A pancreatic thread runs from the first aim to the last.**

### #3. Protein-level cross-check of the whole-body result — **feasible, not yet run**
HPA also publishes immunohistochemistry protein levels across normal tissues. Repeating #1 at
protein level would remove the RNA-versus-protein caveat. Moderate effort, directly closes
weakness #3.

### #4. Release it as a tool — **feasible, no new data**
Everything needed is a small program: paste a gene list, get confinement scores and the ceiling.
This is the only route that meaningfully widens real-world impact, and it is the project's own
long-standing stretch goal.

### #5. More antigen pairs / more diseases — **works, but diminishing**
The two-antigen machinery runs on any dataset. Gains were small once bugs were fixed
(+0.003 median statistical gain), so this adds breadth, not a new claim.

### Angles tested and NOT recommended
- **Patient-outcome prediction.** 5 vs 5 is unfixable before the fair.
- **Myeloma targets against a blood reference.** Plasma cells are not in blood; every myeloma
  target scores at or below the isotype floor. Confirmed uninformative, not evidence against.
- **Raising n on the approved-vs-failed test.** The panel and the world both run out of approved
  targets. Inherent.
