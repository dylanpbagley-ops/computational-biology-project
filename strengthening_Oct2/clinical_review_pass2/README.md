# Clinical evidence repair, pass 2

October 2, 2026. Advisory research record; not submission prose.

This pass verifies selected clinical reports and separates study observations from target-level judgments. It does not replace the original labels, adjudicate a new endpoint, or calculate a new performance estimate. The assistant already knows the historical scores; this is not independent or blinded review.

## Scope and selection

Priority sampling covered disputed original examples, contrasting studies of the same target, and expansion candidates in the existing literature collection. This is a purposive source-repair pass, not a systematic literature review or an outcome-independent validation sample. Sources were consulted on October 2, 2026. Abstract-only extraction is explicitly marked. Searching a previously collected paper again does not make it new validation evidence.

`clinical_evidence.tsv` contains 15 evidence units from 12 publications across 9 target genes. Several units summarize event families or cohorts; they are not a completed individual-event database. No pooled event rate should be calculated from these rows. `study_source_form.tsv` is a separate, neutral 12-publication source inventory with blank reviewer fields. It contains no assistant outcome interpretations, original labels, or model scores. Supply the existing reviewer instructions with it; withhold the assistant evidence table until an independent first assessment is returned.

The five original-table genes are CA9, CD44, EGFR, EPCAM, and ERBB2. The other four are CDH6, CLDN6, PTK7, and ROR1; all four already appear on the locked prediction list. These older publications are not newly observed forward outcomes. Eleven rows of the 166-row prior source index have an exact target/PMID match to this partial extraction; that is coverage linkage, not eleven completed adjudications. Some additional publications concern targets already in the source pool.

The row-level observations and source locations are in `clinical_evidence.tsv`; the extraction input is `curated_evidence.json`. `coverage_and_checks.json` records linkage and integrity checks. Those checks establish file consistency, not clinical correctness.

## What changed scientifically

- **Severity, dose limitation, and target-mediated normal-tissue injury are different outcomes.** The EGFR study distinguishes mild target-mediated effects from severe conditioning-related blood effects. The CLDN6 study reports dose-limiting events while explicitly reporting no observed on-target/off-tumor injury. These are not contradictory statements. Sources: [EGFR primary report](https://aacrjournals.org/clincancerres/article/doi/10.1158/1078-0432.CCR-17-0432/274703/p/Phase-I-Study-of-Chimeric-Antigen-Receptor), [CLDN6 primary report](https://pmc.ncbi.nlm.nih.gov/articles/PMC10667102/).
- **The same gene does not identify a single clinical outcome.** The ERBB2 case report describes a fatal event with a proposed normal-lung mechanism; another ERBB2 study reports no dose-limiting toxicity in 19 patients. CA9 also has different observations across cohorts of the same study. These observations constrain interpretation of a gene-only score; they do not invalidate a possible target-level risk association. Sources: [Morgan](https://pubmed.ncbi.nlm.nih.gov/20179677/), [Ahmed](https://doi.org/10.1200/JCO.2014.58.0225), [Lamers](https://pubmed.ncbi.nlm.nih.gov/23423337/).
- **A halted study need not establish on-target injury.** The CDH6 study ended following neurologic adverse events, but the investigators did not identify their mechanism. ROR1's fatal pulmonary event also leaves normal-tissue recognition unresolved. Both must retain uncertain attribution. Sources: [CDH6](https://pubmed.ncbi.nlm.nih.gov/34515215/), [ROR1](https://pmc.ncbi.nlm.nih.gov/articles/PMC11788652/).
- **Even counts need source reconciliation.** The ROR1 results and discussion differ in their count of reversible grade 4 pulmonary events (two versus three). The extraction flags this discrepancy and does not turn it into a pooled incidence estimate. The fatal respiratory event is recorded separately from that count. [Primary report](https://aacrjournals.org/clincancerres/article/31/3/503/751210/Phase-I-Study-of-ROR1-Specific-CAR-T-Cells-in).
- **A fatal event cited by another study is not another death.** The breast-cancer bivatuzumab report references the fatal skin event in a parallel study. The table preserves the origin and does not attribute that death to the breast-cancer cohort. CD44v6 is also an isoform-specific antigen, not interchangeable with total CD44 expression. Sources: [Tijink](https://pubmed.ncbi.nlm.nih.gov/17062682/), [Rupp](https://pubmed.ncbi.nlm.nih.gov/17351401/).

## Recommended continuation within this project

The original result remains a retrospective target-level association requiring label repair. A useful extension is to evaluate how consistently that association survives explicit event definitions and treatment contexts. This preserves the expression analysis, baseline comparisons, sensitivity work, and relevant mathematical results. It is not yet evidence of a new successful predictor.

Use the coding guide in this folder to finish an event-level extraction across the full original outcome panel and a defined expansion frame. Resolve full-text gaps, adverse-event denominators, follow-up, antigen identity, and duplicate cohorts. Retain all screened records and reasons for exclusion. Select the endpoint before joining amended judgments to expression scores. Independent review should check both coding and whether a defensible negative category exists.

The first amended analysis should use the existing frozen score, confinement/tau, and blood-lineage baseline. Report endpoint sensitivity and within-modality performance; keep all records of a target together in any resampling or train/test division. Do not treat extra publications about the same gene as extra independent target observations. Distinguish previously known retrospective sources from genuinely subsequent outcomes under the original prediction lock.

**Next actionable work:** complete source extraction for the remaining original-panel targets using this schema, then resolve the reviewer questions in `coding_guide.md`. There is no need to wait for lab access to do this. A finished source-repair pass is not the same as finished independent outcome validation.
