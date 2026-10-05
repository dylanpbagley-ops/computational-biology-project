# Scope Narrowing Plan — September 29, 2026

Goal (yours): one clear accomplishment; cut anything peripheral or weakening; show only the most important limitations.

## 1. The one accomplishment

> **Before any drug or patient exists, normal-tissue data alone predicts which antigen targets cause dose-limiting on-target toxicity: a target is safe when it is confined to few normal cell types and its home cell type is expendable. A theorem explains why (the ROC-hull ceiling), and the rule replicates in independent solid-tumour data where the standard tumour-vs-normal comparison fails.**

Everything in the project is now either (a) a step in proving that sentence, (b) a check that it isn't an artefact, or (c) cut.
Working title idea (yours to choose): *Predicting on-target toxicity of cancer immunotherapy targets from normal tissue alone.*

## 2. What stays in the main paper (the "spine")

| # | Piece | Headline numbers (live) | Results-doc § |
|---|---|---|---|
| 1 | Theorem: ceiling = ROC hull vs. weighted normal mixture (short; proof in appendix) | 0 violations in 16,400 checks | 29 |
| 2 | Central test: confinement separates 19 safe vs 11 toxic targets | 0.880, p = 0.00014 | 23.1 |
| 3 | Adding expendability (haematopoietic home) | 0.914, p = 0.00003; direct test 8/19 safe vs 0/11 toxic, Fisher p = 0.013 | 23.2, Strengthening §2 |
| 4 | Not one target / not blood-only | leave-one-out 30/30 significant; solid-only 0.812, p = 0.013 | 23.3, 23.6 |
| 5 | Independent replication, other lab, matched normals, beats field standard | confinement 0.778 (p = 0.025, CI 0.52–0.96) vs standard ceiling ≈ 0.56 | 27 |
| 6 | Not "popularity": registry control | predicts toxicity, not trial advancement (0.448, chance) | 27.3 |
| 7 | Second normal reference + protein level | 0.751 (HPA bulk); 0.785 (antibody protein) | 23.7, 26.1 |

Items 5–7 are one Results subsection ("It is not an artefact"), one figure.

## 2b. Framing the expendability claim accurately (important)
Confinement carries the ranking (0.880). Expendability is a **real association** (p ≈ 0.01, no toxic target has a haematopoietic home) but a **small AUC gain** (+0.033, driven by CD22/CD33/CD38). Say exactly that. Overclaiming it is the one thing that could weaken the spine; stating it precisely does not.

## 3. Move to a short appendix (kept, not featured)

- Efficacy-fail specificity control (§23.4): now 0.746, p = 0.011. It no longer cleanly passes; the registry control (item 6) is the cleaner, externally compiled specificity test and does the same job. Keep the efficacy-fail result in the appendix with its two-cluster explanation. Do not delete it from the record.
- Mass-spec and whole-organ failure (§26): reduced to **one scope sentence** in Discussion ("the rule needs cell-type resolution; whole-organ averaging removes the information").
- Germ/placental sensitivity (§23.5), bug-fix history (substring bug, hull bug): one appendix table. Keep the honest bug notes short.
- Weight-perturbation / hygiene work: only if the screen is kept (see §4).

## 4. Cut from the main narrative

| Item | Why it goes | Where it lives |
|---|---|---|
| Discovery screen + DCT dossier (Angles 1, 4) | Different claim (hypothesis generation), soft below rank ~10, needs internalisation/epitope/dependence data we lack; invites the "so what's the drug?" question | Optional single Discussion sentence: "applied genome-wide, it recovers 10 of 16 known pursued antigens (57-fold)". Drop the candidate list from the paper/poster; keep the files. |
| Clinical-units ceiling (Angle 2) | Measures gating purity, not biology | One methods caution; keep T-cell numbers only if a reviewer asks |
| Cas9 variants (Angle 5) | Precommitted test failed the premise | 1–2 sentences under "scope of theorem" (or nothing) |
| Two-antigen extension (§17, 21, 28) and per-patient windows (§5) | A second theorem-sized story; dilutes the spine | Appendix / future work |
| Myeloma replication (§15), healthy PBMC/BM references (§16, 33) | Support the framework, not the safety claim | One supporting sentence at most |
| EV-capture Aims 1–9, ranking flip, n* calibrations, saturation curve | Origin story, not evidence | One sentence in Introduction |
| Angle 5 / DepMap use | Nothing depends on it | DepMap file stored, unused |

Nothing is deleted from the repository; "cut" means not in the paper, poster or 5-minute pitch.

## 5. Limitations to show (four, all stated once, plainly)

1. **Small panel, one person's labels:** 30 labelled targets (19 safe, 11 toxic); labels from literature, frozen and cited; outside review pending. Leave-one-out shows no single label drives the result.
2. **Expression, not drug behaviour:** confinement is necessary evidence, not a guarantee (a bound, not a prediction of clinical success); RNA and protein agree only weakly at gene level, though both reach the same conclusion.
3. **Needs cell-type resolution:** it fails at whole-organ resolution, as the theory predicts.
4. **Expendability adds modestly** to ranking (above).

Cut from limitations (still disclosed in appendix): efficacy-fail control drift, saturating CLDN18.2 curve (mention only as the theorem's scope), gating purity, screen instability.

## 6. What's left to do (in order)

1. Fix the four stale numbers (§23.4, §23.7, §27.1–27.4; "214"→ cut or 128). The verifier goes to 0 failures.
2. Confirm the abstract/pitch uses the one-sentence claim above.
3. Rewrite the paper skeleton around the spine (Results = items 1–7 only). I can do this next.
4. Send the label sheet to Professor Wideman; submit the OSF pre-registration; update the AI-use log.
5. Write the paper and poster in your own words.

DepMap: moved to `ev_capture_data/depmap/CRISPRGeneEffect.csv` (SHA-256 recorded next to it). No current analysis uses it.

## Addendum (Sept 29, later): two decisions after review
1. **Binary expendability stays primary.** Continuous weights (hand-assigned per cell type) were tested: 0.809 vs binary 0.914; regeneration-only 0.842; did not beat randomly shuffled weights (p = 0.15). Cause: the computed "home" cell type is unreliable for broadly expressed targets (EGFR -> placental tissue, scored fully expendable). The theorem's generality is stated as a scope note. Externally sourced weights were tried afterwards (item 4).
2. **Framework-vs-field material restored as framing, not new evidence** (`Framework_vs_Field_Actions_Sept26.md`): CD47/magrolimab ($4.9B vs last-place confinement) as the hook; HER2 CAR-T; "successes break the field's own criterion" (8/19 vs 0/11); the head-to-head (now 0.556 vs 0.778). Tissue-naming (about 4 of 11) kept as an optional caveated bonus; the misattributed-failure idea (Category E) stays out as speculative. The Sept 26 doc still shows the old numbers (17 targets, 0.528/0.750/0.495); use 18 targets, 0.556/0.778/0.448.

3. **Update (Sept 29, later): second continuous-weight round also failed the adoption rule.** Four theorem-based schemes using the whole expression distribution (E1–E4, fixed in advance) scored 0.818, 0.699, 0.885, 0.880 vs binary 0.914; E3/E4 equal the confinement-only level and do not beat shuffled weights (p ≈ 0.3). Binary stays primary. Results doc section 36; pinned by `tests/test_continuous_weights_v2.py`; verifier now 50 checks, all passing.

4. **Update (Sept 29, later still): third round, weights from an outside dataset (measured cell lifespans, Sender & Milo 2021), also failed.** Six precommitted schemes scored 0.498-0.833 vs binary 0.914; turnover alone is at chance (0.54, 0.50), and adding it to confinement lowers the ranking below confinement alone (0.880). Reason: turnover measures replacement speed, not survivability of the gap (B cells, 63 d, look less expendable than colon lining, 3.4 d, which is backwards). Binary stays primary. Results doc section 36.2; pinned by `tests/test_external_turnover_weights.py`; verifier now 57 checks, all passing. Twelve alternative schemes in total have been tried against the binary flag; record that in the methods.

5. **Update (Sept 29, night): the blood flag was cleaned up, and a search for a more effective rule found nothing better.** (a) The flag is now defined at lineage level (blood is the largest of ten lineage families) instead of "the single top cell type is blood"; it flags exactly the same 30-panel targets (8 safe, 0 toxic), so every headline number is unchanged (0.880 / 0.914), and it no longer depends on one arbitrary cell type for broad genes. Use this definition in the paper. (b) 44 further alternatives (entropy, top-3, lineage-level confinement, absolute abundance, vital-mass, absolute-threshold breadth; tau turned out to be rank-identical to top share) were tried with a family-wise correction; none beat 0.914, none promoted. (c) The remaining errors (16 of 209 pairs; CEACAM5 alone accounts for 8) look like modality (CAR-T / T-cell engagers respond to very low antigen levels) and tissue criticality, neither testable on 30 targets. Results doc sections 37-37.3; pinned by `tests/test_rule_search.py`; verifier now 73 checks. Next real gain: new labelled targets (record the modality), and the locked-prediction / OSF step.

6. **Update (Sept 30, early): new labelled targets, noise simulation, independent naming test, and a correction to the theorem bridge.** (a) Two literature searches found only 7 clean new targets (5 safe, 2 toxic) plus 5 hedged toxic; the frozen rule scores 0.875 on the scorable 4 v 2 (0.900 with a post-hoc absent-gene convention); far too small to be evidence alone (`Validation_Targets_Sept29.md`). (b) A simulation shows a genuinely better rule (by 0.02-0.04) scores lower than a coarser one on a 19 v 11 panel 8-22% of the time and is almost never detected (3-6%), so a lower score on 30 targets does not show a rule is worse; detecting 0.04 needs roughly 150-300 labelled targets. (c) A tissue-location test (does the model put the harm where it fell?) favours the hand-weighted model (median rank 4 v 7; better on 8 of 13, worse on 2; p = 0.11, suggestive not significant). Nothing beat the frozen rule; 52 alternatives tried in all. (d) **Correction:** section 29's claim that top_share is the "leading term" of the theorem's harm formula is too strong (the share-weighted formula is constant when weights are equal); confinement follows from the threshold form and is an empirical finding the theorem makes plausible. Use that wording. Results doc sections 38-38.5 and the correction at section 29; verifier 93 checks. (e) Label consistency: NECTIN4 is labelled safe despite a fatal-skin-reaction boxed warning; disclose that labels are agent-specific.

7. **Update (Sept 29, second session): are the expendability weights the bottleneck, and what is the most realistic rule? (Results doc section 39.)** (a) **No.** The best possible lineage-level weights, fitted with hindsight, raise separation only from 0.914 to 0.923; with cross-validation the same idea drops to 0.76-0.80. (b) A weight measured from outside data (FDA adverse-event reports: share of organ-damage reports that end in death, 28 organ systems) is uncorrelated with the hand weights (Spearman -0.008) because it measures acute lethality, not irreversible loss; in every rule form it was indistinguishable from randomly shuffled weights. (c) The rule I had declared most realistic (worst organ, exposure relative to the tumour level, measured weights) scored 0.694, clearly below the frozen rule, so my prior belief in that form was wrong on this panel. (d) Adding up (not maximising) relative exposure with the simple binary weights scored 0.895, statistically tied with the frozen rule (0.914; predicted in advance): this is the theorem's threshold form and is the best current restatement of the frozen rule; present it as a concordance check (family-wise p = 0.209), not a discovery. (e) Graded weights, hand or measured, do move the documented damaged tissue from about rank 9 to rank 4-5 of 154 (not significant, 13 targets): they describe where damage lands better, not whether a target is dose-limiting. (f) **74 alternatives have now been tried against the frozen rule; state that number in the methods.** The frozen rule stays primary. Verifier 198 checks; pinned by `tests/test_faers_and_relative.py`.

8. **Update (Sept 29, third session): searched for the data a better rule might need (Round 7; `Round7_Data_Search_Findings_Sept29.md`).** (a) Found: drug type, dose and toxicity per agent for all 30 panel and 12 validation targets (55 FDA-label rows, 111 trial-paper rows; 77 genes in total, 35 new); TCGA tumor and tumor-adjacent-normal expression for solid tumors (all 42 genes; myeloma, lymphoma, leukaemia not obtainable); UniProt membrane topology for all 42 genes. (b) **Main finding:** several frozen "toxic" labels are weakly supported or format-dependent in the primary papers (EGFR, EPCAM, NCAM1, CD47, FAP, ERBB2, CEACAM5's CAR-T vs TCR-T); the same target can be tolerated in one format, route or dose and toxic in another. This is agent-collected and needs a human check before any conclusion; **no label was changed and nothing was scored.** (c) FDA labels almost never attribute organ toxicity to target expression, so on-target attribution must come from trial papers. (d) Proposed next steps (each to be precommitted, and each counts toward the 74 alternatives): strict-label sensitivity analysis, an agent-level table with modality, a tumor/normal fold-change form on about 20 solid-tumor genes, and scaling up the labelled set.

9. **Update (Sept 29, fourth session): two more checks, no new rule (count stays 74). Results doc section 40.** (a) **Label fragility:** every set of 1-6 swapped labels was tried. No 3 swapped labels remove significance; significance is first lost at 4 swaps (349 of 27,405 sets), and all 349 move an approved blood-cancer target from safe to toxic. If 4 of 30 labels were wrong at random, 98.7% of sets still give p < 0.05. Swapping all eight weak labels (EGFR, EPCAM, NCAM1, FAP, CD47, ERBB2, CEACAM5 to safe; NECTIN4 to toxic) gives separation 0.744, p = 0.048; 97.3% of the 256 subsets stay below 0.05. (b) **Proved lemma:** top_share <= 1/(1 + k r) when k other cell types exceed r times the peak, so top_share is the tightest guarantee on exposure breadth and bounds the threshold-form harm; one-directional; loose on real genes; Spearman with real exposed count -0.92 to -0.95. This partly repairs the section 29 correction. Verifier 221 checks; pinned by `tests/test_fragility_and_bound.py`.

10. **Update (Sept 29, fifth session): locked-predictions check (results doc section 41), no new rule (count stays 74).** Precommitted and hashed first. Of 32 called targets: 0 forward (class A) verdicts, 5 scoreable pre-lock (class B) verdicts, 2 hits and 3 misses (CD70, LRRC15 hit; ITGB6, ROR1, LY6E miss, all weak ADC-tolerability evidence); 13 hedged, 14 unresolved. Below the floor of 8, so no conclusion either way. Paper wording: the lock is a standing design, not a result.
