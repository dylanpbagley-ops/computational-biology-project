# ISEF Paper Skeleton (rewritten Sept 29, 2026 around the one core claim)

**Replaces the Sept 10 skeleton.** That version had stale numbers (44-target headline panel, the 214-candidate screen as a main result, "leading term" wording for the theorem, and no mention of the 74 alternatives or the label limitations). Everything here follows `Scope_Narrowing_Plan_Sept29.md` and the results doc (`Aim13_Results_Two_Comparator_Ceiling.md`, §-numbers below).

**This is structure, not submittable text.** Society for Science bars generative AI from writing the words in your abstract, research plan or poster, so every sentence in the paper has to be yours. This document tells you what goes where, which number backs which claim, and where to check it. Claude's job after this is feedback on *your* drafts (facts, logic, what a judge will push on), plus figures and structure.

---

## Part 1. The one claim, and what kind of claim each part is

> **Before any drug or patient exists, normal-tissue data alone predicts which antigen targets cause dose-limiting on-target toxicity: a target is safe when it is confined to few normal cell types and its main lineage is one the body can spare (blood). A theorem explains what any therapy can achieve; the rule is an empirical result that replicates in independent solid-tumour data where the standard tumour-vs-normal comparison fails.**

Say each piece at exactly its own strength. This table is the honesty backbone of the paper:

| Piece | What it is | Strength |
|---|---|---|
| ROC-hull ceiling (§29) | Theorem (Neyman–Pearson), checked on 16,400 random cases with 0 violations | Proof plus numerical check |
| Confinement predicts safe vs toxic (§23.1) | Empirical result, 30 labelled targets | Strong statistics, small panel |
| Blood flag adds to it (§23.2) | Real association (no toxic target is blood-dominated) but a small gain in ranking | Modest, say so |
| Theorem explains confinement | **Partly proved, partly motivation.** §29 correction: the share-weighted formula does not reward confinement. §40.2 proves that top_share is the tightest possible guarantee on how many cell types can be exposed (harm in the threshold form is at most 1 + (1/top_share − 1)/r). One-directional: high top_share guarantees limited exposure; low top_share does not guarantee wide exposure | Proved lemma plus empirical rule |
| Relative-exposure sum form scores 0.895 (§39.5) | Concordance check, predicted in advance, family-wise p = 0.209 | Not a discovery |
| Replication in breast-cancer autopsy data (§27) | Independent lab, matched normals | Real, but 18 targets, CI 0.52–0.96 |

---

## Part 2. The story in six beats (tell it out loud before writing)

1. **Problem.** Targeted therapies (CAR-T, ADCs, bispecifics) kill any cell carrying the antigen. Most antigens sit on some normal tissue too. Sometimes that is survivable, sometimes fatal.
2. **What the field does.** Picks targets by tumour-vs-normal difference (fold-change). Hook: the field's greatest successes (CD19, CD20, BCMA…) have essentially no tumour specificity, and magrolimab/CD47 had the field's biggest bet ($4.9B, Gilead–Forty Seven, 2020) on a target the rule ranks last. Source: `Framework_vs_Field_Actions_Sept26.md`.
3. **Theory.** For any therapy whose kill probability only rises with antigen level, the best achievable (tumour killed, normal harmed) trade-off is bounded by the ROC hull of tumour vs a weighted mixture of normal cell types. Engineering slides you along the curve; it cannot move the curve.
4. **Prediction and test.** If normal-tissue distribution sets the ceiling, then normal-tissue data alone should separate targets with safe records from those with dose-limiting toxicity. Test: 30 labelled targets, 154 normal cell types, rule frozen before comparison.
5. **Is it an artefact?** Leave-one-out, solid tumours only, second reference, protein level, independent autopsy cohort where the field-standard measure fails, registry control.
6. **What it does not do, stated plainly.** Small panel; agent-specific labels; no drug type or tumour data yet; 74 alternative rules tried and none beat it.

---

## Part 3. Section-by-section skeleton

Required by the fair: 250-word Abstract, then Introduction, Methods, Results, Discussion, Conclusions, References (one page). Write in this order: Methods, Results, Discussion, Conclusion, Introduction, Abstract.

### 3.1 Abstract (250 words, write last)
Slots only: question; one-sentence method; the two headline numbers (0.880 and 0.914 with 19 v 11); the independent replication (0.778 vs 0.556); the limitation in one clause. No numbers beyond those.

### 3.2 Introduction
- **Clinical stakes, concrete cases** from the label file: HER2 CAR-T fatal event (Morgan 2010); FOLR1/mirvetuximab boxed warning (ocular); tisotumab vedotin (F3) boxed ocular warning; CEACAM5 colitis. Caveat: use the ones your labels support best (see Part 4).
- **Field practice and gap:** fold-change selection; Framework-vs-Field Categories C and D (successes break the field's criterion; the head-to-head).
- **Prior work, named honestly:** SCAN-ACT, Dannenfelser et al. logic-gate screen, the 2025 Frontiers CAR-T modelling review; none gives a computable ceiling from normal tissue alone (novelty search in results doc §7).
- **Hypothesis, stated as falsifiable.**
- **One sentence** on origin (EV-capture model, Aims 1–9). Bridge only.

### 3.3 Methods
1. **Theory (short; full proof to appendix).** Kill rule f(x) in [0,1], monotone in antigen; harm H = Σ w_c q_c E[f]; weighted mixture; hull statement (§29.1–29.2). *Say:* the theorem gives the ceiling and shows why expendability matters; confinement is an empirical finding the threshold form makes plausible (§29 correction box). Do **not** say top_share is the theorem's leading term. Add the lemma (§40.2, proof is three lines): if k other cell types exceed r × the peak, then top_share ≤ 1/(1 + k·r), so top_share upper-bounds exposure breadth and hence threshold-form harm. Say it is one-directional and that the bound is loose on real genes (real exposed count is about 20–30% of the bound) though top_share tracks real breadth closely (Spearman −0.92 to −0.95).
2. **The rule.** `top_share` = largest cell type's share of a gene's expression across 154 HPA single-cell types. Blood flag = blood is the largest of ten lineage families (`lineage_map.py`). Score = top_share + 1 if flagged; higher = predicted safer. No fitted parameter. Flag assigned from cell-type identity alone.
3. **Labels.** 19 safe, 11 toxic (main test); efficacy-fail targets kept as a side check; excluded ones listed. Labels from literature only, frozen in one file, cited. State plainly: **labels are agent-specific verdicts** (a target can be tolerated as one drug format, route or dose and toxic as another), and made by one person. Sources: `Label_Review_Sheet.md`, `Validation_Targets_Sept29.md`.
4. **Statistics.** Separation = AUC (fraction of the 209 safe/toxic pairs where the safe target scores higher; ties half). Exact/permutation p; leave-one-out; shuffled-label null; paired bootstrap intervals for comparing rules.
5. **Search discipline (required).** The rule was frozen before any alternative. **74 alternatives were then tried** (52 + 13 measured-severity variants + 9 relative-exposure variants, §38.5 and §39.7); none met the promotion rule (separation ≥ 0.944 and family-wise p < 0.05). Say that a lower score on 30 targets does not show a rule is worse (§38: detecting a 0.04 gain needs roughly 150–300 labelled targets).
6. **Data.** HPA single-cell RNA (154 types, primary); HPA bulk RNA (51 tissues) and HPA antibody staining (protein) as checks; De Schepper et al. 2025 rapid-autopsy cohort (30 patients, matched normals, 72 ADC targets, 623-trial registry). All public, all collected by others.
7. **Pre-registration.** OSF predictions submitted by you before checking; state what was locked and when.

### 3.4 Results (spine only; numbers in Part 5)
- **R1 Theorem.** 0 violations in 16,400 checks; verifier bug (hull loop) caught and fixed (§29.4). One figure: the hull (`fig_hull.png`).
- **R2 Confinement separates safe from toxic.** 0.880 (p = 0.00014); the ladder CD20 → BCMA → PSMA → CLDN18.2 → CEACAM5 → F3 → EGFR → HER2 → CD47 (§23.1). Figure: ladder.
- **R3 Adding the blood flag.** 0.914 (p = 0.00003); flag alone 0.711 (8 safe flagged, 0 toxic, Fisher p ≈ 0.013); gain of +0.033 is modest. CEACAM5 is the honest exception that motivates it. Report the substring bug (0.916 inflated → 0.914 corrected).
- **R4 Not one target, not blood-only, not label-fragile.** Leave-one-out min 0.905 / max 0.947, 30 of 30 significant; label fragility (§40.1): no 3 swapped labels remove significance, fragility index 4, every significance-losing 4-swap set moves an approved blood target to toxic, and swapping all 7+1 weak labels still gives 0.744 (p = 0.048); solid tumours only 0.812 (exact p = 0.0133, 8 v 10); germ/placental removed 0.914 unchanged.
- **R5 It is not an artefact (one subsection, one figure).** HPA bulk 0.751; antibody protein 0.785; autopsy cohort: confinement 0.778 (p = 0.025, CI 0.52–0.96) vs tumour-vs-normal ceiling 0.556 (p = 0.365), same 18 targets; registry control: confinement does not predict reaching phase 3 (0.448) but does predict toxicity (0.778).
- **R6 What the data cannot settle (this is a result, not an apology).** Weights are not the bottleneck: hindsight-optimal lineage weights reach only 0.923, cross-validated 0.76–0.80; measured FAERS severity weights are indistinguishable from shuffled ones; relative-exposure sum with binary weights reaches 0.895 (predicted in advance, concordance only); graded weights help locate the damaged tissue (median rank 9 → 4–5 of 154, not significant, 13 targets) but not safe/toxic separation (§39). The 18 misordered pairs (of 209) concentrate in CEACAM5 and CA9 among toxic targets and MET/NECTIN4 among safe ones.

Appendix table candidates: efficacy-fail control (0.746, p = 0.011; two-cluster explanation), mass-spec/whole-organ failure (§26), bug history, Round 7 label review.

### 3.5 Discussion
1. What is proven (theorem), what is empirical (rule), what is conjecture (why confinement works).
2. **The reframe:** successes come from targets whose normal home is a lineage the body can spare, not from tumour specificity; why blood cancers got easy targets and solid tumours did not.
3. **Alternatives ruled out, one each:** fold-change (0.556 vs 0.778); blood-only (solid-only 0.812); popularity (registry control); label choice (frozen, cited, leave-one-out).
4. **Scope of the theorem:** needs cell-type resolution (organ-level averaging destroys the information; one sentence, §26); requires monotone kill probability.
5. **Limitations, exactly these (stated once each):**
   - Small panel, one person's labels (but the result survives any 3 wrong labels, §40.1); **several labels are weakly supported or format-dependent in the primary papers** (EGFR, EPCAM, NCAM1, FAP, CD47, ERBB2; CEACAM5's CAR-T vs TCR-T; NECTIN4 is labelled safe despite a boxed warning for fatal skin reactions). You found this yourself in `Round7_Data_Search_Findings_Sept29.md`. Labels were not changed and no strict-label re-scoring was run.
   - Expression, not drug behaviour: no drug type, dose, route, or tumour expression in the rule. **Tumour data and drug type are future work** (tumour data is unobtainable for myeloma, lymphoma, leukaemia; about 11 of 19 safe targets are blood cancers).
   - RNA and protein agree only weakly at gene level (ρ = 0.167) though both reach the same conclusion.
   - Blood flag adds modestly.
6. What a bigger study needs: 150–300 labelled agents with modality recorded; the OSF-locked predictions as the forward test.
7. Optional single sentence: applied genome-wide the rule recovers 10 of 16 known pursued antigens (57-fold, §30.2). Do not show the candidate list.

### 3.6 Conclusion
Short. Thesis at its true strength; the two headline findings; what would settle the open question (more labelled agents, tumour data, drug type); why it matters, modestly.

### 3.7 References
One page. Take bibliographic details from the results doc citation lines, not from memory: De Schepper 2025; Human Protein Atlas (name each file); Neyman–Pearson; prior-art papers; Mizenko 2021 if the origin sentence stays; label sources from the label sheet.

---

## Part 4. Honesty items that must appear (you agreed to all four)
1. **74 alternatives** tried against the frozen rule; none promoted (Methods 5, Results R6).
2. **Labels are agent-specific verdicts**; a few are weakly supported in the primary papers (list above). Presented as a limitation you found.
3. **Tumour data and drug type are future work.**
4. **Theorem-to-rule bridge: partly proved, partly motivation.** The share-weighted "leading term" claim was retracted (§29 correction); what is proved is the breadth bound (§40.2), not that confinement is derived from the theorem.
Also state: the frozen rule is primary and no further rule search was done after this point; you skipped the strict-label sensitivity analysis on purpose (be ready to say why: 30 targets cannot support it, and each variant tested would count toward the 74).

---

## Part 5. Numbers in one place (all current; results doc §)

| Finding | Value | § |
|---|---|---|
| Theorem check | 0 violations / 16,400 | 29.4 |
| Confinement (top_share) | 0.880, p = 0.00014 (n_types 0.876, p = 0.00018), 19 v 11 | 23.1 |
| Frozen rule (top_share + blood flag) | 0.914, p = 0.00003 (≈191 of 209 pairs) | 23.2 |
| Blood flag alone | 0.711; 8 safe, 0 toxic; Fisher p ≈ 0.013 | Strengthening §2 |
| Among 22 non-blood targets | 0.851 (all 18 misordered pairs sit here) | 39 |
| Leave-one-out | 0.905–0.947; 30/30 significant | 23.3 |
| Label fragility | any 3 swaps: p < 0.05 always; 4 swaps: 98.7% of 27,405 sets; fragility index 4; weak labels all swapped 0.744 (p = 0.048), 97.3% of 256 subsets p < 0.05 | 40.1 |
| Breadth lemma | 0 violations (132 gene checks, 20,000 random); Spearman(top_share, exposed types) −0.92 to −0.95 | 40.2 |
| Solid tumours only | 0.812, exact p = 0.0133 (8 v 10) | 23.6 |
| Germ/placental removed | one-factor 0.856; frozen rule 0.914 | 23.5 |
| HPA bulk reference | 0.751, p = 0.0118; ρ with primary 0.592 | 23.7 |
| Antibody protein reference | 0.785, p = 0.0104 | 26.1 |
| Autopsy cohort, same 18 targets | confinement 0.778 (p = 0.025, CI 0.52–0.96) vs tumour-vs-normal 0.556 (p = 0.365) | 27.2, 27.4 |
| Registry control | phase 3: 0.448 (p = 0.675); toxicity 0.778 | 27.3 |
| Efficacy-fail control (appendix) | 0.746, p = 0.011 | 23.4 |
| Oracle lineage weights | 0.923 in-sample; 0.756 cross-validated | 39.2 |
| Relative-exposure sum, binary weights | 0.895 (paired diff −0.019, CI −0.086 to +0.053; family-wise p = 0.209) | 39.5 |
| Best graded-weight forms | 0.876 hand, 0.861 severe, 0.852 death | 39.5 |
| Naming test (13 toxic targets) | median rank 9.0 unweighted; 4.0–5.0 graded; not significant | 39.6 |
| Alternatives tried | 74 | 38.5, 39.7 |
| New validation targets | frozen rule 0.875 on 4 v 2; too small to be evidence | Validation_Targets |
| Panel composition | 19 safe, 11 toxic in the main test; the results doc counts 44 targets overall (12 efficacy-fail, 2 excluded). Check this total against the frozen label file before you write it | 23 |

---

## Part 6. What is cut from the paper (kept in the repo)
Discovery screen and candidate list (one sentence only); clinical-units ceiling; Cas9 variants; two-antigen extension and per-patient windows; myeloma, PBMC and bone-marrow references (at most one supporting sentence); EV-capture Aims 1–9, n* calibrations, saturation curve (one intro sentence). Category E (misattributed efficacy failures) stays out as speculation. Tissue-naming (about 4 of 11) is an optional caveated bonus, never a headline.

---

## Part 7. Compact glossary (learn these; explain them in your own words)
- **Antigen / on-target toxicity:** the protein the drug grabs / harm from the drug working as designed on normal cells that carry it.
- **Separation (AUC):** the fraction of the 209 safe-toxic pairs where the safe target scores higher; 0.5 = chance, 1.0 = perfect.
- **Permutation p-value:** shuffle labels many times; p = fraction of shuffles scoring at least as well.
- **`top_share`:** the top cell type's share of a gene's total expression across 154 types; near 1 = confined.
- **Blood flag:** blood is the gene's largest lineage family (of ten).
- **ROC hull / ceiling:** the boundary no monotone therapy can beat on the tumour-killed vs normal-harmed trade-off.
- **Monotone:** more antigen never lowers kill probability.
- **Frozen / pre-registered:** rule and predictions fixed and time-stamped before checking.
- **Family-wise p:** p-value corrected for having tried many alternatives (74 here).
- **Concordance check:** a result predicted in advance that agrees with the primary rule; it supports the story but is not a new discovery.

---

## Part 8. Questions a judge will ask (have a plain answer ready)
1. "Couldn't you have chosen labels to fit?" (frozen, cited, leave-one-out; and you disclose the weak ones)
2. "Isn't this just blood cancers?" (solid-only 0.812)
3. "You tried 74 rules; isn't 0.914 cherry-picked?" (frozen first; none beat it; family-wise correction)
4. "Why not a better rule with more data?" (the 30-target panel cannot detect a 0.04 gain; Round 7 shows label noise and format dependence)
5. "So what should a drug developer do with this?" (screen out broadly expressed targets early; it is a bound, not a guarantee)
6. "Is the theorem what predicts safety?" (the theorem gives the ceiling; confinement is an empirical finding it makes plausible)

## Part 9. Your to-do (not Claude's)
- Update the AI-use log.
- Send the label sheet and `Validation_Targets_Sept29.md` to Professor Wideman when you can reach him.
- Draft Methods first, in your own words; send it for feedback.
