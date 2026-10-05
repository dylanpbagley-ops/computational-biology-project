# CAR-T Antigen-Selection Pivot — Full Implementation Plan

**Status:** Proposed. Awaiting Dylan's approval before execution.
**Basis:** `PIVOT-FOR-NOVELTY.md` (prior session's research), re-verified this session via fresh web search (the "no tool convolves per-cell antigen-density distribution with the nonlinear kill curve for target selection" gap still holds as of August 2026 — a 2026 Nature Communications AI-guided-CAR paper addresses signaling-domain engineering and pathway modulation, a different axis, not distribution-aware target selection by density).

**Decision already made per the pivot document's own Recommendation 2**: keep all existing EV work intact as the validation backbone proving the general principle; the CAR-T application is presented as the higher-impact use of the same machinery, not a replacement.

---

## 1. The headline claim this pivot is built to produce

> **"A distribution-aware, superselectivity-based scoring method — applied to real per-cell antigen-density data — predicts which candidate CAR-T target antigen (or antigen pair) minimizes the escape fraction (the fraction of tumor cells too antigen-low to be reliably killed), and this ranking differs from the ranking produced by the field's current methods (mean expression, or binary percent-positive)."**

This mirrors the EV project's own Aim 9 headline discipline exactly: state one new, checkable claim, validate the supporting machinery, then point it at the claim. Framed the same honest way Aim 9 was — a real reordering/sensitivity finding on real data, not a claim to have "solved" antigen escape.

---

## 2. What gets reused directly (no rewrite needed)

| Existing file | Reused as | Change needed |
|---|---|---|
| `binding_model.py` — `q_bound`, `capture_probability`, `alpha_exact`, `selectivity_alpha` | The mathematical core: nonlinear response to receptor/antigen count | None to the math itself. New wrapper functions with CAR-T-domain names (see Section 3) that call the same underlying equations. |
| `binding_model.py` — `accessible_receptors` / corona-multiplier pattern | Template for a "not all antigen is functionally engageable" adjustment (e.g., antigen internalization/shedding) | Reused as a pattern, likely a new analogous multiplier for CAR-T-specific effects (see Section 5, stretch item). |
| `population_model.py` | Population-level averaging (mean capture/kill probability across a distribution of receptor counts) | Same, renamed conceptually. |
| `aim3_aim4_range_sweep.py`, `k_intra_z_plausibility_bounds.py` | Pattern for sweeping an uncalibrated parameter across a plausibility range rather than pretending to know it exactly | Directly reusable pattern; new plausibility bounds will need real anchoring data (see Section 4 — this is actually easier than the EV case, because CAR-T dose-response curves are directly, publicly measured, unlike K_INTRA). |
| `aim5_and_gate_dual_marker_selectivity.py` | AND-gate / dual-marker logic | Directly reusable concept — maps onto logic-gated CAR-T (a real, active research area: "AND-gate" CARs that only activate given two antigens, used to spare healthy cells). |
| Statistical tooling from `aim9_marker_ranking_order_check.py` (`auc_from_two_samples`, ranking-order comparison logic, negative-binomial reconstruction-from-summary-statistics convention) | Directly reusable for comparing "escape-fraction ranking" vs. "mean/positivity ranking" across candidate antigens | None — same statistical building blocks, new inputs. |
| The project's whole testing discipline (one test file per module, run every `__main__` after any change, never silently clip an impossible result) | Same discipline, unchanged | None — this is a working practice, not code. |

**Estimated reuse: roughly 70% of existing code/math carries over unchanged in form**, consistent with the pivot document's own estimate. This is not a rewrite — it's a relabeling-and-extension of a working system.

---

## 3. What gets built new

### 3.1 The kill-probability model (`aim10_car_t_kill_probability_model.py`)
- A renamed/extended wrapper around `capture_probability()`: instead of "vesicle with n receptors gets captured by a ligand surface," it's "tumor cell with n antigen copies gets killed by a CAR-T cell."
- Critically, CAR-T's dose-response curve is **directly measured and published** (Majzner et al. 2020, Walker et al. 2017 — already cited in `binding_model.py`), which is a real methodological upgrade over the EV project: instead of an uncalibrated K_INTRA swept across a plausibility range, we can fit the model's nonlinear response curve directly to real titrated-antigen-density-vs-killing data. This closes exactly the kind of calibration gap that stayed permanently open for K_INTRA in the EV project. This will be flagged honestly as a real methodological improvement in the write-up, not overclaimed as "fully solved" (dose-response data still has noise and is drawn from a specific engineered cell system, not universal).
- Output: given a full population of tumor cells' antigen-count distribution (not just the mean), compute the **escape fraction** — the fraction of cells below the functional kill threshold — as the model's core deliverable.

### 3.2 The ranking/scoring module (`aim11_antigen_target_ranking.py`)
- Takes candidate antigens (CD19, CD22, BCMA, and any others with public per-cell density data), and for each:
  - Computes the "naive" ranking (by mean expression, and separately by percent-positive — replicating exactly what the pivot document found current tools like the Nature Biotech logic-gate atlas and LogiCAR actually do).
  - Computes the "distribution-aware" ranking (by predicted escape fraction, using the real per-cell distribution and the fitted kill curve from 3.1).
  - Reports where these rankings disagree — this is the direct analog of Aim 9's "9 of 104 grid points reorder the four CRC markers" finding, now for a decision with real clinical stakes.
- Extends into AND-gate pairs using the `aim5` pattern: for two antigens targeted simultaneously, what's the joint escape fraction, and does it beat either single-antigen ranking?

### 3.3 Real data acquisition and preprocessing (`ev_capture_data/car_t_data/`, mirroring existing structure)
- Candidate real datasets to acquire and vet (in priority order, to be confirmed once we're inside implementation — this project's own standing discipline says "verify a claimed dataset independently before building on it," so exact data acceptance happens during implementation, not assumed now):
  1. **GSE241783** — a genuine single-cell CITE-seq atlas of CD19 CAR-T patients (large B-cell lymphoma, axicabtagene ciloleucel), confirmed this session to be publicly available via GEO/CellXGene/EGA. This is the primary target: real per-cell CD19 protein measurements from real patients.
  2. Quantitative flow cytometry antigen-density summary statistics already partially surfaced in the pivot document (e.g., median CD19 ~5,810 molecules/cell, IQR 3,442–8,725; dSTORM range 13–5,000 molecules/cell across patients) — usable the same way Guo et al.'s summary AUCs were used in Aim 9: real numbers, used honestly as summary-statistic-level evidence, clearly distinguished from raw single-cell data.
  3. Published titrated dose-response curves (Majzner 2020, Walker 2017) for fitting the kill-probability function itself.
- Each dataset gets the same provenance discipline as the EV project: confirm what it actually measures (protein vs. RNA proxy — CITE-seq protein counts are strongly preferred over RNA-only scRNA-seq, since RNA is an imperfect proxy for surface antigen density, a caveat the pivot document itself flags), confirm real grouping variables (e.g., responder vs. relapser, pre- vs. post-treatment), and never assume a figure/summary stat means more than it actually reports.

### 3.4 Validation and honest reporting (`docs/CAR_T_Pivot_Validation_and_Findings.md`)
- Same standard as `Prediction_Registry.md`: state the claim before checking it, report the real result whichever way it comes out (including if the ranking-reordering finding is small, absent, or inconclusive with the actual data found — this project has a strong, credible track record of reporting negative/partial results honestly, e.g. Aim 9's Check 2 and the two Prediction Registry checks, and that discipline continues here).

---

## 4. What is genuinely easier than the EV project (a real upgrade, not just a lateral move)

- **The kill-response curve is directly measurable from public data** (titrated antigen density vs. killing efficiency experiments exist), unlike K_INTRA, which stayed a permanently uncalibrated placeholder for the entire EV project. This is a real, structural improvement in rigor.
- **Real single-cell (not just summary-statistic) antigen-density data is more readily public** (GSE241783 and similar CAR-T atlases) than single-vesicle EV data was — the EV project had to fight through SVAtlas's hidden raw-data page, controlled-access approval walls (OMIX010139), and CORS/network blocks just to get two raw datasets. CAR-T single-cell atlases are a much more mature, better-supported public data ecosystem (GEO, CellxGene).

## 5. What is a real, honestly-flagged limitation of the pivot (say this in the write-up, don't bury it)

- This does not address CAR-T's other major failure modes: T-cell exhaustion, trafficking/infiltration into solid tumors, immunosuppressive tumor microenvironment. The claim is narrowly "better target prioritization by density-distribution," not "solves CAR-T."
- Antigen density for solid-tumor targets is often only measurable at the RNA level in public atlases, while the real clinical decision depends on protein-level density — CITE-seq (protein+RNA together) is the mitigation, and should be prioritized as the primary validation dataset over RNA-only atlases for exactly this reason.
- The superselectivity-to-immunology crossover is genuinely novel/early (currently opinion-piece level in the literature) — this is a strength (real novelty) framed honestly as also meaning "no existing sub-literature to lean on," so extra rigor is needed in explaining why the physics transfers.
- (Stretch, not required for the core claim) A CAR-T-specific analog of the EV project's "accessible_receptors" corona-occlusion multiplier could represent antigen shedding/masking/internalization — flagged as a real, optional extension, not a blocking requirement.

---

## 6. Concrete build order

1. **Data acquisition and vetting** — retrieve GSE241783 (or the best-vetted alternative found during implementation) and the quantitative antigen-density summary statistics; confirm exactly what each measures and its real grouping variables, the same way every EV dataset was vetted before use. Document what's usable and what isn't, honestly, before writing model code against it.
2. **Fit the kill-response curve** — use the real titrated dose-response data (Majzner/Walker) to anchor the nonlinear response function's parameters directly, replacing the EV project's permanently-swept K_INTRA/z with the closest thing to a real calibration this whole line of work has ever had.
3. **Build `aim10_car_t_kill_probability_model.py`** — the renamed/extended core model, with its own test file, following the existing one-test-file-per-module convention exactly.
4. **Build `aim11_antigen_target_ranking.py`** — the naive-vs-distribution-aware ranking comparison, plus the AND-gate dual-antigen extension reusing `aim5`'s logic.
5. **Run the real validation** — compute both rankings on real data, report the honest result (reorders or doesn't, and by how much), exactly as Aim 9 did for CRC markers.
6. **Update the project's headline framing and write-up materials** — extend `EV-Capture-CURRENT.md`'s "project in one sentence" pattern with the new headline claim from Section 1 above, present EV work as the validated foundation, CAR-T as the higher-impact application, per the pivot document's Recommendation 2.
7. **Full test suite run** — every new module gets its own test file; the full existing 149-test suite plus new tests must all pass; every new module's own `__main__` entry point run end-to-end, per this project's own standing discipline.
8. **Update `docs/AI_Use_Disclosure_Log.md`** with this session's entry, per ISEF's AI-disclosure requirement and this project's own consistent practice of logging exactly what Claude was used for each session.

---

## 7. What does NOT change

- No existing EV code, data, or test is deleted or overwritten. `ev_capture_model/src/aim1` through `aim9` and their tests remain exactly as they are — this is an *addition* (new aims 10+), not a replacement.
- The "What NOT to revisit" list in `EV-Capture-CURRENT.md` remains binding for all EV-specific settled decisions (K_INTRA/z/k_geom permanently uncalibrated for the EV case specifically, OMIX010139 deprioritized, etc.) — none of that reopens because of this pivot.
- The project's core working-style rules (plain language, explain fully, stress-test before executing, acknowledge corrections, report honest negative results, never narrow scope without justification) apply identically to this new work.

---

## 8. Open questions I will resolve during implementation, not guess at now

- Exact final choice and access method for the primary validation dataset (GSE241783 vs. an alternative found during real vetting) — decided during Section 6, Step 1, with the same honesty discipline as every EV dataset decision.
- Whether CITE-seq protein data for CD19 in GSE241783 (or the chosen dataset) has sufficient per-cell resolution and a real responder/relapser or pre/post grouping variable to support the ranking-reordering check, or whether it turns out (like several EV leads) to be a real, relevant, but structurally limited dataset — reported honestly either way.
- Exact set of candidate antigens beyond CD19 (CD22, BCMA, others) depends on what real per-cell density data is actually publicly available and vetted, not assumed in advance.

---

**This plan implements the pivot document's Recommendation 1 (primary pivot) and Recommendation 2 (keep EV work as backbone) in full. Recommendations 3 (MRD fallback) and 4 (do not pivot to scRNA-seq marker selection) are explicitly NOT pursued, per the pivot document's own reasoning: MRD/flow-cytometry is lower novelty and harder to beat entrenched consortia solo, and scRNA-seq marker selection is a crowded subfield that already treats distribution-awareness as its founding premise.**
