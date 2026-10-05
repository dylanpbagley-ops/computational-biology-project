# Prediction Registry

**Date registered: August 30, 2026.**

**Why this document exists:** the project's own novelty diagnosis (`docs/Novelty_Competitiveness_Diagnosis_and_Restructuring_Options.md`, Restructuring Option C) identified a real structural weakness: every check the project has done so far follows the pattern "find data, then see if the model fits it." That order of operations is weaker than the alternative: state a specific, falsifiable prediction FIRST, before knowing whether any data exists to check it, then go looking. This document does that. Everything below was written before any new search was performed to check it. Predictions are stated plainly, in the order they occurred to me, with no attempt to only list the ones I suspected would be checkable.

This document does not replace Aim 9. Aim 9 already found and reported a real result (the model reorders CRC markers relative to naive linear ranking at 8.7% of the plausibility grid). The predictions below are NEW, broader claims the model implies, not yet checked against anything.

---

## Prediction 1: markers with lower absolute receptor/copy count should show LARGER nonlinear-vs-linear ranking divergence than high-copy-count markers, for the same fold-change between conditions

**The specific claim:** the model's capture probability, `capture_probability(n, sigma, K_INTRA, z) = 1 - exp(-sigma * K_INTRA * z * n)`-family nonlinearity (see `binding_model.py`), means the relative IMPACT of a fixed fold-change in receptor count depends on where that count sits on the curve. A marker with a LOW mean receptor count per vesicle (say, tens of copies) sits in the steep, high-curvature part of the capture-probability curve; a marker with a HIGH mean receptor count (hundreds to thousands of copies) sits closer to the saturating plateau. This predicts: for two markers with the SAME fold-change between healthy and cancer populations, the LOWER-copy-number marker should show a LARGER difference between its "naive fold-change ranking" and its "nonlinear model-predicted ranking" than the higher-copy-number marker.

**Why this is a real, new, falsifiable claim:** this is a testable, quantitative, directional prediction about WHICH markers should show the reordering effect Aim 9 already found exists in aggregate (8.7% of the grid) — a prediction about WHERE that 8.7% concentrates, not just that it exists.

**Can this currently be checked?** Partially. This project's own data (ITGB3/ITGAM/ITGAL/ITGB2 from the CRC raw file, the tetraspanins CD9/CD63/CD81 from MASEV) spans a real range of positivity fractions/copy numbers already in hand. A future session could test this directly using data already in this project's possession, without any new search — this is flagged as the single most immediately actionable item in this registry.

**CHECKED, same day this registry was written — an honest "cannot meaningfully test with this data" result, not a pass or fail:** using the project's standing negative-binomial reconstruction convention (the same one Aim 9 uses) on the four real CRC integrin markers' positivity fractions, the reconstructed mean receptor counts came out as ITGAL=3.18, ITGB3=3.89, ITGB2=4.02, ITGAM=4.26 -- all within a narrow ~35% range of each other, not spanning a genuine "low-copy vs. high-copy regime" contrast. This is NOT evidence against Prediction 1 -- it means these four specific markers, in this specific dataset, don't differ enough in reconstructed copy number to distinguish the prediction one way or the other. A real test of Prediction 1 needs markers spanning a much wider copy-number range (e.g., tetraspanins, which are known to be present at higher absolute copy numbers than these integrins, per Aim 6/8's own data) -- comparing an MASEV tetraspanin against a CRC integrin would cross datasets/methods and introduce new confounds, so this remains a real, open, currently-unresolved item, reported honestly rather than forced to a conclusion the data doesn't support.

---

## Prediction 2: an EV panel that combines a LOW-copy-number marker with a HIGH-copy-number marker (an AND-gate) should show a LARGER selectivity gain than combining two markers of similar copy number

**The specific claim:** building on Aim 5's existing AND-gate model, this predicts that the "combination gain" (the improvement in tumor-vs-normal selectivity from requiring both markers vs. either alone) should be asymmetric — larger when the two markers being combined sit at different points on their respective nonlinear capture curves (one steep, one flatter), because the combination captures information from two structurally different regimes rather than two redundant ones.

**Why this is new:** Aim 3b's existing co-occurrence work found the "combination gain" number highly sensitive to an unconstrained absolute copy-number scale (a known, documented limitation — see "Known risks" in the handoff doc) and explicitly avoided citing a single point value for this reason. This prediction reframes that sensitivity as itself the finding: gain should scale with how DIFFERENT the two markers' copy-number regimes are, not be a single fixed number. This turns a known modeling limitation into a specific, testable directional claim.

**Can this currently be checked?** Not with a real two-marker joint distribution currently in this project's data (the same limitation Aim 5's own docstring already documents — the correlation between markers is assumed/swept, not measured). This is an honest example of a prediction that currently has NO known way to check it — included anyway, per this document's own stated purpose, rather than quietly dropped for being inconvenient.

---

## Prediction 3: the model predicts that a marker's DIAGNOSTIC AUC and a marker's RAW POSITIVITY FRACTION should NOT be monotonically related across a large enough marker panel

**The specific claim:** naive intuition says "a marker present in more vesicles overall should be a better diagnostic marker" (higher positivity → higher AUC). The model predicts this should NOT generally hold — a marker with a modest positivity fraction but a large SEPARATION between healthy and cancer distributions can outperform (in AUC) a marker with a higher overall positivity fraction but weaker separation, because AUC depends on distribution OVERLAP, not raw magnitude.

**Why this is new and different from Aim 9:** Aim 9's Check 1 tests whether the NONLINEAR MODEL reorders markers relative to LINEAR positivity. This prediction is about REAL, PUBLISHED data directly — does published AUC data itself (independent of this project's own model) already show non-monotonicity with positivity/prevalence? This is checkable using ONLY external published numbers, no model computation required, making it the "cheapest" prediction in this registry to check.

**Can this currently be checked?** Yes, immediately, using data already gathered this session: Guo et al.'s four real AUCs (ITGAM=0.9158, ITGB3=0.8681, EPS8=0.8498, CD9=0.8242) combined with the CRC raw file's positivity fractions for the two markers that overlap (ITGAM, ITGB3). A fuller check would need positivity data for EPS8/CD9 too, which — as Aim 9's own docstring already documents — does not exist in this project. This prediction can be partially checked now, fully checked only if EPS8/CD9 positivity data is ever found.

**CHECKED, same day this registry was written — a real, honest, negative result for the two-marker case currently in hand:** ITGAM has BOTH the higher real AUC (0.9158 vs. ITGB3's 0.8681) AND the higher real raw positivity fraction (14.62% vs. ITGB3's 10.98%). For these two specific markers, AUC and positivity ARE monotonically related — the opposite of what this prediction hypothesized. This is reported plainly as a real negative finding, not hidden or reframed. It does NOT disprove the broader prediction (n=2 markers is not enough to establish or refute a general non-monotonicity claim, and EPS8/CD9 — the two markers where this project has no positivity data — could behave differently), but it means the two markers currently checkable do NOT exhibit the predicted non-monotonicity. Report this honestly in any write-up that cites this registry: the cheapest, most immediately checkable prediction in this document was checked and did not hold for the available data, while the broader question remains genuinely open pending EPS8/CD9 data.

---

## Prediction 4: a real, upcoming single-EV cancer dataset — if one becomes available before the fair — will show that the model's predicted marker ranking (nonlinear yield ranking) agrees with real diagnostic performance MORE OFTEN for markers with lower mean copy number than for markers with higher mean copy number

**The specific claim:** this is a direct extension of Prediction 1 to a genuinely prospective (not-yet-existing-data) setting. It restates Prediction 1 as a claim about a FUTURE dataset rather than data already in hand — this is the prediction most suited to Option E (a publicly time-stamped, prospective claim), since by construction no data currently exists to check it.

**Why this is registered even though nothing can check it yet:** per this document's own principle (state the prediction before knowing if data exists), this is registered now, dated today, specifically so that if a matching dataset appears before the fair, this project has a genuine, dated, prospective claim to check it against — not a claim invented after the fact to fit whatever the new data happens to show. See the companion note on Option E below.

**Can this currently be checked?** No — by design. This is the one true "wait and see" prediction in this registry.

---

## What this registry does NOT do

This document does not claim any of the four predictions above are validated. None of them have been checked as of the date of this document. Prediction 3 is the cheapest to check and should be attempted first, using data already in this project's possession. Prediction 1 is the next most actionable. Prediction 2 is honestly currently unable to be checked with any real data this project has access to. Prediction 4 is a genuinely prospective claim, registered specifically to support a possible Option E (public prospective registration) pursuit — see the implementation plan for the honest assessment of whether Option E's timeline is currently workable (short answer: the mechanism — OSF preregistration, free and immediate — is workable; what is NOT currently in hand is an identifiable upcoming dataset to check Prediction 4 against, since the most promising 2025-2026 EV marker-ranking papers found this session are already fully published, not upcoming).

## Recommended next action

Prediction 3 should be checked first — it requires no new data gathering, only combining numbers already in this project's possession. If a future session does this, it should be checked BEFORE reading this conclusion any further, to preserve the predict-first discipline this whole document is trying to establish.
