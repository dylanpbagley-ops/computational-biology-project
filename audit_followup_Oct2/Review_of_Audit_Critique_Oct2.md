# Review of the audit critique — October 2, 2026

Advisory material for Dylan, not student-authored submission text. Compared `Audit_of_the_Audit_Oct2.md` with `independent_audit_Oct1/Independent_Assessment_Oct1.md`, the operational labels and unrounded scores, and accessible primary publications. Historical labels, rules, predictions, and audit files were preserved. This is a post hoc review, not a new preregistered validation.

**My revised view:** the critique supplies a valuable missing sensitivity analysis and additional relevant prior work. I accept those additions. It does not overturn the original assessment: the project contains a strong retrospective association, with unresolved endpoint, novelty, and independent-validation limitations. Describing the project as having no worthwhile result would be wrong. Describing it as an established general safety predictor would also be wrong. Neither assessment can establish that it will or will not qualify for ISEF.

## The important addition I accept

I independently reproduced all seven scenarios using unrounded `wholebody_v2.json` values and the existing blood flag. Exact one-sided rank-test probabilities were calculated by counting possible rank allocations, without a statistical-library approximation.

| Scenario | Safe / toxic | AUC | Exact one-sided p |
|---|---:|---:|---:|
| Original | 19 / 11 | 0.913876 | 0.000027862 |
| Exclude FAP | 19 / 10 | 0.910526 | 0.000056765 |
| Exclude FOLR1 | 19 / 10 | 0.915789 | 0.000043435 |
| Exclude both | 19 / 9 | 0.912281 | 0.000092516 |
| Exclude FAP; classify FOLR1 as safe | 20 / 9 | 0.905556 | 0.000109036 |
| Previous scenario, also classify NECTIN4 as toxic | 19 / 10 | 0.915789 | 0.000043435 |
| Exclude FAP; FOLR1 safe; exclude F3 | 20 / 8 | 0.893750 | 0.000329783 |

**I should have included this targeted sensitivity table in the first audit.** It shows that these disputed entries do not explain away the observed ranking. The original audit already verified general label-flip/deletion robustness, but that was not a substitute for showing the specific source-driven scenarios.

These are hypothetical label scenarios, not my completed adjudication of all clinical records. Approval plus rejection of one on-target attribution does not by itself establish every part of the project's safe-label definition. A blinded, consistent review of species, agent, endpoint, severity, and attribution remains necessary.

The critique's explanation that correcting these labels simply removes errors is incomplete. Excluding both leaves 15 discordant pairs out of 171; excluding FAP and moving FOLR1 to safe leaves 17 out of 180. The latter AUC decreases slightly. The robust conclusion is that the association remains strong—not that correction necessarily improves it.

Small conditional p-values do not settle target selection, modality confounding, endpoint validity, or the effect of using the discovery panel to develop the rule. Nor are these sensitivity scenarios an independent replication or global correction for the project's exploration history. The original audit's independent-validation concerns therefore remain.

Reproduction: `recheck.py`, `label_sensitivity.csv`, and `recheck_results.json` in this folder. Input hashes are recorded in the JSON. Higher scores mean safe outranking toxic; AUC is not classification accuracy.

## Mathematics: agreement, with a useful extension

Both documents agree on the actual repair: distinguish arbitrary likelihood-ratio rules from rules constrained to increase with antigen, and remove the claim that one AUC determines the entire attainable frontier.

The critique correctly notes that my second example concerns antigen-ranking AUC. The original audit explicitly called it antigen-ranking AUC, so this is a useful scope clarification rather than a false counterexample. To remove any ambiguity about the broader claim, I checked a second, exact example:

- Normal masses: (1/3, 1/3, 1/3).
- Two alternative positive distributions: (0.60, 0.30, 0.10) and (0.55, 0.40, 0.05).
- Both likelihood-ratio ROC areas are exactly 2/3.
- At false-positive fraction 1/3, their optimal true-positive fractions are 0.60 and 0.55.

Thus equal likelihood-ratio AUC also does not imply equal frontiers. Exact fractions are in the reproduction output.

The project's earlier recognition of non-monotonicity is useful existing work and should be retained. It means some of the correct interpretation is already available. It does not make the contradictory later theorem wording correct, establish novelty of the classical likelihood-ratio principle, or turn annotation-based separation into independent functional validation. This repair need not require a new model, but it does change the scope of what the existing verification establishes.

## Source and novelty claims need careful distinctions

**SCAN-ACT is close prior art.** Its methods use normal-cell expression breadth and exclude lymphoid cells partly because depletion may be better tolerated. That supports the critique's warning about conceptual overlap. However, its cell-type-count filter and exclusions are not literally the same quantities as this project's largest-expression-share score and blood-home bonus. The concepts are established; whether this exact empirical implementation adds useful knowledge requires a fair comparison. [Primary methods](https://link.springer.com/article/10.1186/s13073-025-01514-9).

**Madan et al. is relevant, but the comparator matters.** The published Figure 2 caption compares safety-score distributions for clinical CAR targets with all surfaceome genes. That is not a tolerated-versus-toxic clinical-outcome discrimination test. Therefore the critique's cited p = 0.95 cannot be presented as another method failing the same clinical endpoint. I verified the comparison from the PubMed figure caption; the full-text endpoint for independently checking that exact p-value was inaccessible in this pass. [Publication and figure captions](https://pubmed.ncbi.nlm.nih.gov/37835579/).

**BayesTS already supplies a normal-reference score and a discrimination statistic.** The paper reports AUPR 0.59 for prioritizing a curated tumor-specific/safe proxy set, assembled from several sources including expression-based categories. That is not the same endpoint as this project's binary historical toxicity labels, but it substantially narrows broad novelty claims about normal-only scoring and quantitative benchmarking. [Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11705768/).

**Thomas et al. needs respectful treatment as related work.** The accessible abstract reports no significant association between target expression and reported clinical toxicities. The critique's final section distinguishes its predictors and endpoints from this project's; that is a better position than calling it a direct refutation. However, I did not independently obtain the full paper supplied to the other assistant, so I am not certifying its detailed claims about every metric, endpoint, or target count. Those full-text characterizations remain attributed to the critique. [Primary publication record](https://pubmed.ncbi.nlm.nih.gov/40375426/).

I do not accept “nobody has done this” as established by this literature search. A plausible narrow contribution is a transparent retrospective evaluation of a simple normal-expression score against carefully curated outcomes. Whether it is novel enough or useful enough for a particular journal is still a research judgment. Merely adding an AUC to a related biological idea does not establish substantial methodological novelty. Likewise, the shared-target check described in the critique is reuse of existing labels on an overlapping subset, not independent validation or evidence of beating Thomas et al.

## Two criticisms of my citations

The alleged FAP title error is not present in the supplied `Independent_Assessment_Oct1.md`. Its link label is “FAP CAR-T marrow toxicity and cachexia, Tran et al.” and points to the relevant JEM article. The disputed “to treat murine thymoma” title does not appear. I cannot substantiate that allegation in the file being reviewed; a different earlier version would need separate evidence.

“Corneal mechanism study” is an imprecise description of the Matulonis paper. It is primarily a prophylaxis study that includes mechanistic findings. I accept using its full title or a more precise description. That does not reverse the distinction between a serious ocular adverse event and proof of on-target ocular injury. [Publication](https://pubmed.ncbi.nlm.nih.gov/30413525/).

## What remains unchanged, and what this means for the project

The two assessments largely agree on the theorem mismatch, head-to-head uncertainty, stale frozen-table provenance, modality/source errors, and the need to narrow claims. Passing software tests and matching hashes are valuable evidence of computational discipline; they do not independently establish the scientific assertions around them. Corrections should preserve the original dated artifacts and use an explicit amendment, rather than making a newly generated file appear historically frozen.

The critique is right that substantial improvements to the current submission can be made without a new model or new dataset. But improved sourcing and writing cannot create independent outcome validation or establish superiority that the current comparison does not show. “All fixable in weeks” is reasonable for many document repairs, not a demonstrated timetable for resolving every scientific limitation.

The old numerical award probabilities lack a defensible calibration method. They should be withdrawn, not retained as approximately valid because nobody supplied replacement percentages. The useful assessment is conditional: a corrected exploratory project can be worthwhile; high-level competitiveness depends on novelty, evidence, execution, and the student's ownership and understanding.

**Practical consequence:** keep the current work as a defensible fallback and acknowledge the stronger sensitivity evidence. Continue a bounded search for an application with cleaner, independently testable outcomes because that is Dylan's stated goal—not because the original association has disappeared. The option-1/option-2 feasibility findings are reported separately in `Feasibility_Findings_Oct2.md`. So far they justify continued testing of imaging measurement reliability, not abandoning the old work or declaring a superior replacement already achieved.
