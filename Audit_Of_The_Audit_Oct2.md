# Checking the independent assessment — October 2, 2026

Advisory material for Dylan, not submission text. I re-derived every checkable claim in `Independent_Assessment_Oct1.md` from the project's own code and data, and had the clinical and prior-art claims checked against primary sources. Commands and numbers below are reproducible.

## Verdict in one paragraph

The assessment is substantially correct, and on several points it is more accurate than the assessments I gave you in September. Its two mathematical counterexamples are exact. Its three label criticisms are real. Its statistical criticism of the head-to-head is right. Its record-keeping findings are right, including one error in a file I wrote for you yesterday. It also contains errors of its own, including a wrong citation title in the very paragraph where it criticises your citations. And on one point it is too *generous*: the prior-art situation is worse than it says. But it leaves out the most important number, which I computed: **the label errors it found do not move your headline result.** The central empirical finding survives every correction either of us can justify. What needs repair is the framing, the theorem wording, the sourcing, and the novelty claim.

---

## 1. Confirmed: the mathematics

I recomputed both counterexamples in exact rational arithmetic.

**Counterexample 1.** With tumour antigen distribution (0.45, 0.10, 0.45) and normal (0.10, 0.80, 0.10), a rule that kills the low and high antigen bins but spares the middle kills 90% of tumour and 20% of normal cells. Restricted to rules where killing never decreases as antigen rises, the best achievable at 20% normal kill is **23/45 = 51.11%**. Its figure is exact.

**Counterexample 2.** With uniform normal mass over four bins, tumour distributions (0.25, 0.25, 0.25, 0.25) and (0.45, 0.05, 0.05, 0.45) both have antigen-ranking AUC exactly 0.5, yet at 25% normal kill the first allows 25% tumour kill and the second 45%. Exact.

**Does your project actually make the claim being refuted?** Yes, in two places:
- `weighted_ceiling_theorem.py` line 51: "the best attainable trade-off is governed by a SINGLE number, AUC(p_T, p_W)."
- Results document line 1774: "the whole achievable region is governed by the single number AUC(p_T, p_W)."

And the conflation is real. Line 1743 introduces the ceiling as holding "for any kill rule monotone in antigen," but the theorem proved at line 1771 is the Neyman–Pearson statement about *all* rules, whose ROC is ordered by the likelihood ratio, not by antigen. Those coincide only when the likelihood ratio happens to increase with antigen. Your verifier enumerates all 2^n subsets of antigen bins (`weighted_ceiling_theorem.py` line 125), so its 16,400 passing checks confirm the all-rules statement, not the monotone one. The assessment is right about this.

**Two refinements it did not make, both in your favour.**

First, its counterexample 2 refutes the *antigen-ranking* reading of "single number AUC." Under the likelihood-ratio reading those two distributions have different AUCs (0.5 and 0.70), so that example does not touch it. The "single number" claim is still false under either reading — AUC is one scalar summarising an entire curve, and two curves with equal area can differ at any given point — but the specific demonstration offered is narrower than presented.

Second, and more useful: **your project already contains the correct version of this point.** Results document lines 152–155 state that where the empirical likelihood ratio is non-monotone a likelihood-ratio rule beats *every* density threshold, with measured effects (CD27 by +0.79 Youden's J, CD48 +0.38, CD38 +0.10, and exactly 0.000 for CD19, CD20, CD22, CD79B). Line 148 names your measured CLDN18.2 non-monotonicity as the boundary of the monotone assumption. So §29 contradicts your own earlier section. This is a wording repair, not new mathematics, and it is arguably a *strengthening*: you have measured the gap between the two ceilings on real antigens, which is more than the counterexample does.

## 2. Confirmed: the three label criticisms

Checked against primary sources and against your own `target_labels.py`.

| target | your record says | what the sources say | verdict |
|---|---|---|---|
| FAP | "FAP CAR-T… toxic… on-target toxicity to FAP+ bone-marrow stromal cells; cachexia and marrow failure" | Tran et al., *J Exp Med* 2013;210:1125 (PMID 23712432) is a **mouse** study. No human trial reports dose-limiting on-target FAP marrow toxicity. The two human FAP CAR-T trials were loco-regional and reported no CAR-related toxicity; both sibrotuzumab trials were well tolerated with no MTD reached. | **Confirmed. A mouse result is recorded as clinical history.** This fails your own labelling rule. |
| FOLR1 | "dose-limiting **on-target** ocular keratopathy (boxed warning)… corneal epithelium" | Matulonis et al., *Clin Cancer Res* 2019;25:1727 (PMID 30413525): "FRα expression was **absent** from human corneal tissues… phenotypically consistent with **off-target** effects on the cornea." A 2024 review attributes it to non-receptor-mediated macropinocytosis by corneal stem cells. | **Confirmed, and it is backwards.** Mirvetuximab is approved and its ocular toxicity is off-target, so by your own rule FOLR1 should be *safe*. |
| CEACAM5 | "CEA **CAR-T** (Parkhurst 2011)" | Parkhurst et al., *Mol Ther* 2011;19:620 (PMID 21157437) used **TCR-engineered** T cells, not a CAR. Dose-limiting grade 3 colitis and on-target attribution both hold ("probably due to lymphocyte recognition of the normal levels of CEA present in colonic mucosa"). | **Modality is wrong; the label itself survives.** |

Your own `Round7_Data_Search_Findings_Sept29.md` already flags ERBB2, EGFR, EPCAM, NCAM1 and CD47 for the same class of problem, and already concludes "part of the remaining error in the frozen rule is label noise, not a missing rule." You were ahead of this criticism; you just had not acted on it.

## 3. The number the assessment never computed — and it matters

The assessment calls labels "the largest scientific weakness" but never asks what the errors do to the result. I did, using unrounded top shares from `wholebody_v2.json`:

| scenario | panel | AUC | exact p |
|---|---|---|---|
| as published | 19 v 11 | 0.9139 | 0.00003 |
| drop FAP (mouse-only) | 19 v 10 | 0.9105 | 0.00006 |
| drop FOLR1 (off-target) | 19 v 10 | 0.9158 | 0.00004 |
| drop both | 19 v 9 | 0.9123 | 0.00009 |
| **most defensible: FAP excluded, FOLR1 → safe** | **20 v 9** | **0.9056** | **0.00011** |
| the above, plus NECTIN4 → toxic | 19 v 10 | 0.9158 | 0.00004 |
| the above, plus F3 dropped | 20 v 8 | 0.8938 | 0.00033 |

**The headline survives every correction I can justify, at p < 0.001 throughout.** The reason is instructive and worth saying out loud in your paper: FAP (0.167) and FOLR1 (0.205) score low, which is where the rule puts toxic targets, so they were never carrying the result. Three of the 18 misordered pairs involve them, so correcting them *removes* errors.

I also reproduced every headline figure exactly: AUC 0.9138755981, 191 concordant / 18 discordant / 0 ties / 209 pairs, exact one-sided p 0.0000278615; confinement 0.8803827751 and p 0.0001462456. My bootstrap 95% interval is 0.794–0.995 against its 0.789–0.995. Everything checks.

## 4. Confirmed: the head-to-head does not show superiority

This one I got wrong in September and you should change it. `Framework_vs_Field_Actions_Sept26.md` line 191 calls 0.556 vs 0.778 "the head-to-head." My own paired bootstrap over the same 9-versus-9 panel gives a difference of 0.2222 with a 95% interval of **−0.086 to +0.531**, and 7.7% of draws favour the field comparator. A significant association for one predictor and a non-significant one for the other is not a significant difference between them — that is a standard error and I made it.

Worse for that specific comparison: **six of its nine toxic targets are ones your own records flag** (FAP, FOLR1, EGFR, ERBB2, NCAM1, CEACAM5). It does hold up under correction (0.786 vs 0.529 after FAP is excluded and FOLR1 moved), but it cannot be presented as a demonstration that your framework beats field practice. Present it as directionally consistent and underpowered.

## 5. Confirmed: the record-keeping problems

- **The frozen TSV is not what the analysis uses.** `target_labels_frozen.tsv` has 40 rows; `target_labels.py` has 44. Missing: AXL, ERBB3, F3, GPNMB. **F3 is one of the 11 toxic targets in your scored panel.** Four precommitment hash files pin the 40-row TSV, so the hashes are valid but pin the wrong artifact. This is the single most damaging item on the list, because "the labels are frozen and hashed" is a claim a judge can check in thirty seconds.
- **"16 misordered pairs" should be 18.** It appeared twice in the skeleton. I fixed both. I also checked the surrounding claim still holds: all 18 lie among the 22 non-blood targets, and they do concentrate in CEACAM5 (8) and CA9 (6) among toxic and MET (5) and NECTIN4 (3) among safe, as written.
- **Your verifier tests substring presence** (`verify_documented_numbers.py` line 30: `ok = (s in text)`). So 221 passing checks mean the numbers appear in the document and match the JSONs. They do not audit the sentence around the number. "16 misordered pairs" is the proof: wrong, and never looked at.
- **Test count is 341, not 195.** I collected them; the docs are stale.
- **My error.** The handoff I wrote yesterday described the AUC as "the chance a random toxic target outscores a random safe one." It is the opposite: higher score means safer, mean safe score 0.791 against toxic 0.121, so 0.914 is the chance a *safe* target outranks a toxic one. I have fixed it. That is exactly the kind of slip a judge catches, and it was mine.

## 6. Where the assessment is itself wrong

- **Its FAP citation is wrong.** It gives the title "Immune targeting of fibroblast activation protein to treat murine thymoma." That is a different paper (Lee et al., *Cancer Res* 2005, a dendritic-cell vaccine study) which reported essentially *no* toxicity. The correct title is "…triggers recognition of multipotent bone marrow stromal cells and cachexia." Same journal, volume and pages; wrong title. It made a citation error inside the paragraph criticising your citation errors.
- **It mislabels PMID 30413525** as a "corneal mechanism study." It is primarily a phase I corticosteroid eye-drop prophylaxis study that contains mechanistic work. The substance holds; the description does not.
- **Its counterexample 2 is narrower than claimed**, as above.
- **Its probability criticism is fair but empty.** It is right that prize counts cannot calibrate an individual's chances, and right that I presented 85 / 72 / 22 with more precision than they deserved. But it offers no alternative and no method. Treat my September numbers as rough and too confident at the top end, not as refuted.

## 7. Where it is too generous: the prior art is worse than it says

It names two prior works. A dedicated search found four that matter more, and one of them is not merely prior art but a **contrary published result**:

- **Thomas et al., *Comput Biol Med* 2025 (PMID 40375426).** Built atlases from over 3 million cells across 35 healthy tissues, deliberately chose targets "for which clinical trial data are available," and compared target expression with clinical trial efficacy and toxicity. Finding, verbatim: **"Target expression did not significantly correlate with reported clinical toxicities in patients undergoing therapy."** Bonferroni-adjusted p > 0.05 throughout.
- **Madan et al., *Cancers* 2023 (PMID 37835579).** An explicit single-cell normal-tissue "safety score" from Tabula Sapiens, tested against clinical CAR targets: Wilcoxon **p = 0.95**, authors state "it is not discriminatory."
- **SCAN-ACT, *Genome Medicine* 2025 (10.1186/s13073-025-01514-9).** Already encodes **both of your features** as a coded filter: a count of HPA single-cell cell types with high expression, thresholded at the 75th percentile of clinically tested CAR-T targets, and lymphoid cells excluded "due to the potentially higher tolerability of their depletion." Published before your work.
- **BayesTS, *Cell Reports Methods* 2024 (PMID 39515334).** A normal-tissue-only specificity score *with* a discrimination statistic (AUPR 0.59), but benchmarked on expression-derived proxy labels, not toxicity outcomes.

So three things you may have thought were yours are not: the concept that normal-tissue expression predicts on-target risk, the two features themselves, and calibrating against clinical target lists.

**What is still genuinely unclaimed, after a hard look:** nobody has taken an expression-only score, applied it to a panel labelled by *observed clinical tolerability versus dose-limiting on-target toxicity*, and reported a discrimination statistic. BayesTS has the statistic with proxy labels. Thomas has the labels with no statistic. SCAN-ACT has the panel with neither. **That gap is real, and it is precisely what you have.** Your cell-type-level resolution and the two-feature parsimony are further distinctive points.

This reframing is narrower than "I found that normal tissue predicts toxicity," but it survives contact with the literature, and it gives you an answer to the obvious reviewer question — "didn't Kobold and Marr already test this and find nothing?" Your answers: their panel is almost entirely blood targets, so the lineage feature has almost no variance in it; they correlated expression *magnitude* against graded toxicity *rates*, where you rank a *concentration* metric against *binary* dose-limiting labels. Those are real differences. They must be stated explicitly, with Thomas et al. cited by name, or the work reads as rediscovering a published negative.

## 8. My corrected overall view

I over-rated this project in September, and the specific ways matter more than the overall grade:

- I called §29 a proved unification and the thing to aim at the Horvitz prize. It has a genuine gap that your own earlier section contradicts.
- I framed 0.556 vs 0.778 as your framework beating field practice. The difference is not significant.
- I wrote that novelty was "RESOLVED." It is not, and the honest claim is narrower.
- My probability estimates were too precise and too high at the top.

What I will not concede, because the data do not support it: that the central result is in doubt. Your empirical finding reproduces exactly, survives every label correction either assessment can justify, holds at p < 0.001 throughout, and has a bootstrap interval that never approaches 0.5. Across 341 passing tests and hash-locked precommitments, the *discipline* here is real and unusual for a high-school project. The problems are in what you claim the result *means* — the theorem's scope, the comparison with field practice, the novelty, and the sourcing of five or six labels. Those are all fixable in weeks, by you, without new data and without a new model.

The assessment's judgment — a worthwhile question with a strong retrospective signal, not a demonstrated general safety predictor, with several headline interpretations exceeding the evidence — is closer to right than mine was.

## 9. What to do, in order

1. **Fix the frozen-label provenance.** Regenerate the TSV from `target_labels.py` so all 44 rows are present, re-hash, and write a dated note saying the earlier TSV was missing four targets including F3 and that no label values changed. Do not quietly overwrite it.
2. **Correct the labels, as a documented deviation.** FAP: mouse-only evidence, exclude. FOLR1: approved drug, off-target ocular toxicity, move to safe. CEACAM5: TCR-T, not CAR-T. Then re-run and report the corrected headline (I get 0.9056, p = 0.00011) *alongside* the original, with the table from §3 above as a sensitivity analysis. This turns the criticism into a strength.
3. **Repair the theorem wording.** State two ceilings: the likelihood-ratio ceiling over all rules (what you proved and verified), and the monotone-in-antigen hull (what a CAR actually implements), with the first an upper bound on the second that can be loose. Drop "governed by a single number AUC" in both places. Cite your own lines 148–155 as the measured version. Ask Prof. Wideman to point you to someone mathematical for a read.
4. **Downgrade the head-to-head** to directionally consistent and underpowered, with the paired interval shown.
5. **Rewrite the novelty paragraph** around the real gap, and cite Thomas et al. and SCAN-ACT by name.
6. **Then write.** Methods first, in your own words. Do not start new analysis.

Nothing here requires new data, a new rule, or widening the project.

---

## 10. Follow-up: does Thomas et al. refute the project? No — but the defence is narrower than it first looks

Added after the first pass, in response to the obvious question.

**The overlap test.** Of the 17 targets in the Thomas panel, 10 appear in this project's labelled panel. Running the frozen rule on only those 10:

| panel | pairs | concordant | AUC | exact p |
|---|---|---|---|---|
| overlap as labelled (8 safe v 2 toxic) | 16 | 16 | 1.000 | 0.0222 |
| drop NCAM1 | 8 | 8 | 1.000 | 0.1111 |
| drop CD44 | 8 | 8 | 1.000 | 0.1111 |
| non-blood members only (2 v 2) | 4 | 4 | 1.000 | 0.1667 |

Every safe target scores above 0.257; both toxic ones below 0.074; the gap is 0.185.

**What this does establish.** The frozen rule ranks the shared targets correctly, so the Thomas null cannot be read as evidence that this predictor fails on those targets. The two studies do not conflict on shared data.

**What it does not establish, and must not be claimed.** The result rests on two toxic targets, **CD44 and NCAM1 — both of which this project's own `Round7_Data_Search_Findings_Sept29.md` flags as questionable attributions.** Removing either leaves p = 0.111. All six blood-lineage members of the overlap are labelled safe, so much of the separation is the blood flag restating field common knowledge. Restricted to non-blood members the direction holds with no power (2 v 2). **This is a consistency check, not a discovery, and presenting it as succeeding where a published study failed would not survive scrutiny.**

**The defensible claim is methodological.** Thomas et al. had clinical outcome data and normal-tissue single-cell data and asked about expression *magnitude*, using Pearson correlations against graded toxicity rates. This project asks about *concentration across cell types* and discriminates a binary dose-limiting outcome by rank. Magnitude and concentration are different quantities, and the evidence that concentration is the informative one comes from the full 30-target panel (AUC 0.906–0.914, p < 0.001), not from the 10-target overlap. Their panel is also almost entirely tolerated haematologic targets — of the 10 shared, 8 are safe and 6 are blood — so the blood-lineage feature has almost no variance in it.

**Open question — NOW ANSWERED, see §12.** (Left in place to show the reasoning order.) The characterisation of the Thomas method above came, at the time of writing this section, from the abstract and publisher page, **not the full text.** If that paper included any breadth, specificity or cell-type-count measure, this defence weakens substantially and the novelty claim narrows further. Obtaining and reading the Thomas full text is the single highest-priority verification task remaining.

---

## 11. Open questions for whoever reads this next

These are the things I could not settle, listed so they are not mistaken for settled:

1. ~~**Thomas et al. full text**~~ — **RESOLVED October 2, see §12.** Read in full; no breadth or concentration measure is used, and the endpoints are CRS/neurotoxicity/hematotoxicity rather than on-target tissue damage.
2. **A complete label audit.** I verified three labels against primary sources. The project has 30 scored labels and its own review flags at least six more (ERBB2, EGFR, EPCAM, NCAM1, CD47, CD44). Someone with clinical-toxicity expertise should review all 30 blind to the expression scores.
3. **The monotone-vs-likelihood-ratio repair** should be read by someone mathematically qualified. The repair itself is wording, but the two-ceiling statement should be checked by a human who knows ROC theory.
4. **Whether the project's own measured non-monotonicity** (results lines 148–155, CLDN18.2 and the likelihood-ratio gains on CD27/CD48/CD38) can be promoted into a positive result rather than a caveat. It may be the most original thing in the document and it is currently buried.
5. **OSF registration** — the handoff states the predictions were submitted; no public URL or timestamp was independently verified in either assessment.
6. **ISEF/AzSEF paperwork**, including whether Student Support Disclosure Form 2A applies. Not verified by me.

---

## 12. RESOLVED: the Thomas et al. full text (read October 2)

The highest-priority open question from §10 and §11 is now answered. I read the full text (Thomas M, Brabenec R, Gregor L, Andreu-Sanz D, Carlini E, Müller PJ, Gottschlich A, Simnica D, Kobold S, Marr C. *Computers in Biology and Medicine* 192 (2025) 110332). The defence in §10 holds, and a stronger one is available.

**What they actually measured.**

| | Thomas et al. 2025 | this project |
|---|---|---|
| predictor | mean normalised gene expression per cell type ("global CAR target expression profiles"); targets ranked by **Euclidean distance in gene expression space to the approved targets** CD19/BCMA, min-max normalised, via `scipy.spatial.distance.cdist` | top cell type's **share** of expression across 154 cell types, plus a blood-lineage flag |
| outcome | **hematotoxicity, neurotoxicity, and CRS**, as counts/rates of patients per trial. Hematotoxicity "extended beyond the conventional ICAHT grading by incorporating additional hematological parameters such as anemia, lymphopenia, thrombocytopenia, **infections**, and leukopenia." CRS graded 1–2 vs 3–4. | binary: dose-limiting on-target toxicity attributed by the authors to the target on a normal tissue |
| test | Pearson correlation, Bonferroni-adjusted | rank discrimination (AUC) with exact Mann–Whitney p |
| panel | 17 targets from FL, MM and B-ALL trials — essentially all blood-lineage | 30 targets spanning blood, epithelial, neural and stromal lineages |

**No concentration, breadth or specificity metric appears anywhere in the paper.** I searched the full text for specificity indices, entropy, tau, cell-type counts and equivalents. The words "confined" and "restricted" occur only as narrative prose describing individual targets ("CD22 expression was mostly confined to B cells"), never as a computed quantity. Their only cross-cell-type construct is similarity-to-CD19/BCMA by Euclidean distance, which is a nearest-known-safe-target measure, not a concentration measure. **The §10 defence is confirmed: they did not test the quantity this project tests.**

**The stronger point, which I had missed.** Their three toxicity endpoints — CRS, neurotoxicity and hematotoxicity-including-infections — are largely **not** on-target normal-tissue destruction. CRS is immune activation; ICANS is not predominantly on-target; infections and cytopenias are heavily driven by lymphodepleting conditioning. Their own discussion says these are "multifactorial in nature, influenced by various factors such as immune system activation, cytokine dysregulation, and interactions with the tumor microenvironment."

So on a panel of almost entirely blood-lineage targets — the ones whose on-target toxicity is survivable aplasia — they correlated expression against the toxicities that are *not* on-target tissue damage. A null there is consistent with this project's framework rather than contrary to it. This project's claim makes no prediction about CRS or ICANS rates at all.

**And their two near-misses point the right way.** They note that "some association between target expression on HSPCs and hematotoxicity (Fig. S2d) as well as target expression on immune cells in the lymph nodes and CRS (Fig. S2h) seems plausible" — the two places where the expression-to-toxicity mechanism is most directly on-target are where they saw the most signal, just not past Bonferroni correction.

**How to write this, without overclaiming.** Thomas et al. is a resource paper, not a refutation of anything, and it should be cited respectfully as the most closely related prior work with clinical outcome data. The honest sentences are: they tested whether expression *magnitude* predicts CRS, neurotoxicity and hematotoxicity, and found no significant correlation; this project tests whether expression *concentration* predicts dose-limiting on-target normal-tissue toxicity, a different predictor and a different endpoint; their panel is almost entirely blood-lineage, so a lineage feature has little variance in it; and the frozen rule ranks the ten shared targets correctly (§10), so the two results do not conflict. Do **not** write that this project succeeded where they failed — they did not attempt this test.

**What this does not change.** The §7 novelty finding stands untouched. SCAN-ACT still encodes both features as an uncalibrated filter, and the honest novelty claim is still the narrow one: those features have never been validated against observed toxicity outcomes with a discrimination statistic. Thomas et al. is no longer a threat to the project's correctness; it was never the threat to its novelty.

---

## 13. Note on provenance of this document

Written by Claude (Anthropic) on October 2, 2026, as advisory material for the student, in response to `Independent_Assessment_Oct1.md` (written by ChatGPT on October 1). It is not student-authored submission text and must not be submitted as such. Every numerical claim above was recomputed from the project's own `wholebody_v2.json`, `flag_cleanup.json`, `target_labels.py`, `adc_autopsy_ceiling.json` and `weighted_ceiling_theorem.py`; the mathematical counterexamples were recomputed in exact rational arithmetic; the clinical and prior-art claims were checked against primary sources. Where I could not verify something, I have said so rather than inferring it. Two factual errors found in Claude-authored working files were corrected in place: the AUC direction in `HANDOFF_START_HERE_Oct1.md`, and the misordered-pair count (16 → 18, two occurrences) in `ISEF_Paper_Skeleton.md`. No label, prediction, result file or frozen artifact was modified.

Update, October 2: §12 was added after the student supplied the Thomas et al. full text as a PDF. That paper was read directly and in full; every characterisation of its methods and endpoints in §12 is quoted or paraphrased from the text itself, not from the abstract.
