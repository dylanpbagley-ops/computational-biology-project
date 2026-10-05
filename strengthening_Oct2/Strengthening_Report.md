# Strengthening the existing biology project

October 2, 2026. Advisory work record. This pass preserves the biological question and historical results while repairing claims, testing specific objections, and preparing the next evidence step. It does not claim that every scientific weakness is resolved or that awards/publication are assured.

**The useful result remains. The best opportunity to strengthen it is now clearer: improve the clinical endpoint and evidence unit, then evaluate the unchanged rule on independently reviewed agent/study outcomes.** This can remain a computational-biology project and does not depend on obtaining an ASU wet lab.

## What was completed

1. Preserved copies and hashes before changes. Reconciled the 40-row frozen TSV with the 44-row operational table: AXL, ERBB3, F3, and GPNMB were omitted from the former; shared outcome values did not differ. Created a separately dated 44-row export. It retains historical evidence statements and is not a newly adjudicated truth table.
2. Corrected the mathematical source documentation and section 29. Added separate implementations and tests for monotone and arbitrary response bounds, explicit zero-weight handling, and counterexamples for frontier/AUC interpretation and simultaneous tissue constraints. The empirical score remains empirical.
3. Ran fixed-rule diagnostics specified before execution, without fitting or selecting a new rule. Added paired comparisons to simple baselines, a descriptive modality check, and broader label-uncertainty diagnostics.
4. Prepared a randomized, score-free review form for all 30 discovery-panel records, with candidate sources and blank reviewer judgments. The old review sheet exposed the original labels and is now marked superseded. A separate key preserves the historical labels.
5. Reviewed primary-source concerns and found an additional exclusion error involving CD45. Logged source corrections separately from uncertain clinical adjudication.
6. Inventoried all 166 prior source rows in `agent_evidence_index.tsv`, including original/validation/locked-list overlap and identity issues. Created `evidence_amendments_Oct2.tsv` with 44 historical records and separately documented source corrections; every final adjudicated outcome remains explicitly pending. Replaced the active outline and added a single current entry point, `PROJECT_STATUS_Oct2.md`.

## New numerical checks

These are post-discovery diagnostics conditional on the existing selected labels. A protocol and hash were recorded before their execution; prior results were already known. No new independent validation is claimed.

| Fixed predictor | AUC on original 19/11 panel |
|---|---:|
| Original specificity-plus-blood rule | 0.913876 |
| Specificity alone | 0.880383 |
| Blood flag alone | 0.710526 |
| Standard tau specificity | 0.880383 |
| Tau plus the same blood flag | 0.913876 |

The paired AUC gain over blood alone is **0.20335**, with a target-resampled 95% bootstrap interval of **0.07895–0.33014** (20,000 draws). This supports the limited statement that the association in this selected panel contains ranking information beyond the binary blood flag. It does not rule out disease-lineage selection or validate clinical prediction.

The gain over specificity alone is **0.03349**, interval **0.00000–0.10048**. The added flag's benefit is modest and not established independently. Tau's identical ranking is an algebraic identity, already acknowledged in the old section 37, not a new empirical discovery. The relevant contribution cannot be invention of a new specificity index.

Within the prespecified coarse modality groups, there are 42 eligible safe/toxic pairs and concordance is **0.92857**. The underlying groups are very small: CAR-T 3/3, antibodies 5/1, ADCs 7/4. Other categories have only one outcome class; three mixed-agent historical rows were excluded from this calculation. Pairs share targets and are not independent observations. This is descriptive evidence that the pattern is not entirely produced by cross-modality comparisons; it is not proof that modality, dose, or selection confounding is controlled. It still includes disputed historical labels.

The broader source-uncertainty test is less reassuring than the two-label correction alone:

| Diagnostic | Panel | AUC | Interpretation |
|---|---:|---:|---|
| Exclude FAP and FOLR1 | 19/9 | 0.91228 | Confirms the prior targeted sensitivity |
| Exclude those plus EGFR, EPCAM, NCAM1, CD47, F3 | 19/4 | 0.81579 | Residual positive association, but very few toxic examples |
| All binary assignments of those seven rows, fixed membership | 128 assignments | 0.70192–0.91388 | An uncertainty envelope; assignments are not equally plausible |

The 19/4 deletion's nominal exact one-sided p is 0.02722. This is not a new strict-label validation: remaining labels have not all been independently adjudicated, the scenario was informed by prior source concerns, and the number is not corrected for all research-history analyses. Deletion and relabeling answer different questions. We should not select the most favorable panel and replace the original headline with it.

## Clinical evidence: the main unresolved issue

The historical definitions do not cleanly separate categories. “Severe,” “dose-limiting,” “approved,” and “manageable with supportive care” can overlap. A clinical study can have dose-limiting events during escalation and still establish a tolerable regimen. A gene-only score gives the same value to agents with different mechanisms and clinical contexts. The outcome unit therefore needs to be explicit.

The score-free review packet recommends separate agent–study/regimen–event fields for severity, dose limitation, mechanism attribution, confidence, and adequacy of negative evidence. It leaves the final endpoint/category decisions blank. A knowledgeable reviewer must settle the definition before a revised analysis is scored. The assistant's prior exposure to scores prevents its triage from counting as independent blinded review.

Verified source corrections include FAP's species, FOLR1's corneal attribution, and Parkhurst's TCR-T identity. Additional findings matter:

- The CD45 exclusion says no cytotoxic agent reached patients, but SIERRA is a randomized phase III anti-CD45 radioconjugate trial. This stated exclusion reason is false. Transplant conditioning may warrant a distinct scope decision; the target should not simply be added as safe. [Primary trial](https://pubmed.ncbi.nlm.nih.gov/39298738/).
- A later CEA CAR-T colitis report exists. Correcting the earlier TCR-T citation cannot justify a blanket assertion that CAR-T colitis was never reported. This is a separate study and must remain separately sourced. [Primary report](https://pubmed.ncbi.nlm.nih.gov/39734708/).

The triage record distinguishes current source checks from inherited, not fully re-reviewed extraction. All 30 main records remain in the independent-review queue, including tolerated examples. There is no completed “corrected clinical truth” table yet.

## Mathematical repairs

The exact weighted-mixture identity requires a specified shared response and weights; it is not itself a clinical harm model. Nondecreasing responses use the antigen-threshold ROC hull. Arbitrary responses use the larger all-subsets feasible set, with a likelihood-ratio upper boundary. Neither entire frontier is specified by its AUC.

A further gap concerned multiple tissues. The minimum of separately optimized tissue ceilings can exceed what is jointly achievable, because each separate optimum can require a different response. The added two-bin example gives individual ceilings of **0.75**, but joint optimum **0.50**. The legacy one-hot/AUC check cannot establish joint attainability. These clarifications strengthen correctness and expose useful assumptions; they do not establish new clinical efficacy or a novel mathematical principle.

The corrected proof and tests are in `mathematical_scope.md`, `frontiers.py`, and `test_frontiers.py`. The historical source/output is preserved; the original numerical algorithm was not changed.

## What expansion is actually available?

The existing Round 7 files contain **166 source rows across 79 target strings**; two strings represent dual targets, so there are 77 single-target strings. There are 26 strings outside the original and validation/exclusion tables, but these include dual targets and a biosynthesis-proxy issue, and 12 already occur on the locked-prediction list. Across the full source pool, 23 strings occur on that list. None of these counts is a count of newly eligible independent outcomes.

This is nevertheless a concrete resource already collected within your project. The most useful extension is a transparent agent/study evidence benchmark with explicit eligibility and uncertainty, followed by evaluation of the unchanged score and relevant existing baselines. It can expose where gene-level specificity is informative and where agent context prevents a universal ranking. It need not involve building a more complicated predictor.

Before scoring an amended dataset: settle the endpoint, review both included and excluded records without scores, preserve source dates and prior exposure, and group any repeated agents by target for evaluation. Cases previously examined cannot be described as pristine validation. Existing locked predictions stay unchanged and must be evaluated using their own recorded rule version and genuine post-lock outcomes.

There is no defensible shortcut from 166 extracted rows to a large validated clinical dataset. The immediate external dependency is expert source/endpoint review, not advanced lab equipment. Relevant expansion is authorized, but no contact with a reviewer has been made.

## What remains unresolved

`weakness_ledger.csv` tracks 23 issues with their status. Major open items are independent endpoint adjudication, new eligible outcomes, selection/multiple-analysis uncertainty, fair comparison with established methods, RNA-to-protein and isoform limitations, reference partition/coverage, and biological calibration of harm weights. The small original panel cannot answer all these questions. More statistical resampling does not supply missing biology or independent subjects.

The external-data interpretation is now corrected: the autopsy comparison has independent comparator expression but reused outcome labels, and its paired difference interval crosses zero. It does not establish superiority. The old prospective results remain inconclusive, not successful prospective validation. Historical result documents still contain prior interpretations and are explicitly identified as historical; the current status, outline, and this report govern interpretation.

Student ownership also remains essential. The work log and understanding prompts identify what was AI-assisted; Dylan should verify, interpret, and write his own submission and disclosure with the Adult Sponsor. Current ISEF guidance requires independent student work and acknowledgment of AI assistance. [Official rules](https://www.societyforscience.org/isef/international-rules/rules-for-all-projects/).

## Verification and files

Fresh verification is recorded in `verification_results.json` and logs: new numerical tests, the affected legacy theorem entry point run in an isolated directory, the legacy documented-number checker, protocol hashes, and preserved historical artifacts. The old number checker checks formatted-value presence; its passing is not a contextual claim audit. The new explicit scope checks supplement it.

Use the score-free form and instructions for independent review. Do not send `review_key_do_not_send.tsv` or model outputs first. The original frozen labels, predictions, and numerical result files remain unchanged. `document_changes.json` records amended working documents; backups are under `before/`.
