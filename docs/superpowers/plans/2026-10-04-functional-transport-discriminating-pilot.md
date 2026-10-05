# Discriminating functional-transport feasibility pilot

> Execute with existing `/root/luna_functional_outcome_search`, gpt-6-luna Max. No additional agents. Parent writes plans and independently verifies; executor implements. This narrowly extends cycle4 from inspected schemas to one diagnostic calculation. It does not select a substantive extension or certify novelty.

**Goal:** Determine whether an exact, potentially useful biological counterexample merits further work, instead of another schema-only or generic agreement-map result.

**Architecture:** Fixed eligibility and descriptive gene-perturbation direction checks on two already downloaded public studies. No model fitting, performance optimization, AUC, clinical inference or outcome-driven parameter changes. Existing prior-art review remains applicable.

**Tech stack:** Existing NumPy/openpyxl/Python only, compressed TSV counts and one XLSX. No installs; scipy/xlrd unavailable. Main plan/addendum constraints and preservation apply.

**Spec:** `2026-10-04-functional-immune-outcomes-data-first.md` and `2026-10-04-functional-search-parallel-review.md` in this directory. Actual data paths below were opened by executor and independently previewed by parent. This is an exploratory feasibility test, not an adopted project analysis.

## Exact scientific question and limits

Within the separately selected Ramos validation gene panel, do perturbations with a consistently sensitizing CTL selection effect across at least three Lawson tumor models instead consistently promote resistance in BOTH Ramos in-vitro CAR-T conditions? Are any supported contrasts distinct from already analyzed comparisons, rather than known IFN/antigen/death/autophagy results, general context dependence or low-count artifacts?

This is a precisely scoped cross-study transport question. Tumor backgrounds, immune modality, antigen, timing and selection protocols differ: an observed contrast cannot identify which difference causes it. The known pathways are calibration/prior-art references, not forbidden research subjects; any worthwhile increment must be an exact untested contrast with biological value. One model's five contrasts are not five independent cancers. The in-vivo 10m-vs15m dose mismatch excludes those columns from the primary test. No cross-study human or clinical inference is allowed.

The schema/access gate passed for this calculation. Scientific originality and whether the extension would be substantial remain unresolved. Pilot is authorized solely to decide whether a meaningful exact contrast exists; a generic comparison or one unexplained gene hit is insufficient.

## Deadline

Main round dispatch04:37:35UTC, total75min INCLUDING parent review, ending05:52:35UTC. Executor start was reported04:38:08UTC (stage creation04:39:29). Finish this pilot AND stage closeout by05:38:08UTC, earlier if resolved. Stop all Dufva/other access checks and searches now. No remaining new resources. Target20min to implement/run plus closeout, no silent extensions.

## Files and preservation

ROOT `C:/Users/dylan.bagley/Desktop/outputs`; STAGE `opportunity_search_Oct4/cycle4_functional`.
Inputs, unchanged:
- `data/GSE149933_mTKO_readCounts_4May20.txt.gz`:94,528guide rows,7annotation fields plus115MUS count columns.
- `data/GSE149933_mTKO_library.txt.gz`:94,528annotation rows,5fields; existing count-file Gene.symbol can be used for gene aggregation. Do not create missing labels from sequence matching.
- `data/GSE149933_family.soft.gz` and derived `data/sample_mapping.tsv`:97mappedcolumns; unmapped18 remain unused. Preserve source/title/antigen context. Six matched CTL models are EMT6,CT26,4T1,B16,MC38,Renca-HA. Antigen-free Renca cytokine arms are NOT pooled with Renca-HA controls.
- `data/Ramos_SupplementaryData1.xlsx`:1,021selectedvalidationgenes; Gene plus5MAGeCKLFC/FDRpairs. Primary only `cd19_1to2_v_egfr_1to2` and `cd19_1to10_v_egfr_1to10`. Read exact workbook headers, do not guess partial-column positions.

Write new files only under `STAGE/pilot/`, and finish main stage records/validator/hashes/live status as already authorized. Add this plan to protected plan references before computation. Suggested new files: `protocol.json`, `inspect_transport.py`, `test_transport.py`, `gene_effects.tsv`, `transport_flags.tsv`, `summary.json`, `decision.md`. Protocol records exact thresholds/claim limits and its SHA before output inspection. Previously read source outcomes are development knowledge; never call this blinded or independent prospective validation.

## Fixed calculation

1. Read source mappings and validate six distinct CTL/control model contexts with two measured timepoints and three culture replicates per arm. Use the numerically latest matched timepoint per model as primary. Keep earlier timepoints for the specified consistency check; different calendar durations are a limitation. Verify exact source contexts, not collapsed cell-line names. If required matched arms are missing, report unavailable; don't invent them.
2. Compute whole-library column totals from all guide rows. For each guide/sample, use `log2(1e6*(count+0.5)/(column_total+0.5*n_guides))`. This is an exploratory library-normalized log frequency, NOT the published drugZ statistic, a new predictor, direct percent killing, or MAGeCK reproduction. Do not compare absolute effect magnitude across methods/studies.
3. Guide measurement filter fixed before results: count at original model-specific T0 >=20; raw count in each of the three primary control samples >=20. A gene requires >=3 retained unique guide IDs. Excluded genes/guides remain in records with reason. Do not reclassify excluded/unmapped/missing as zero effect. This filter addresses noisy depletion/floor effects and restricts interpretation to evaluable guides; not chosen for outcome agreement. Do not alter thresholds after seeing patterns.
4. For each eligible gene/model, calculate all9 treated-culture vs control-culture differences by taking the median guide-level normalized log-frequency difference across retained guides. Report median of9, minimum and maximum. These are DESCRIPTIVE observed-culture envelopes, NOT confidence intervals, independent donor replication or universal bounds. Negative means sensitization/depletion; positive means resistance/enrichment. Label sign-stable sensitizing only if the maximum of9 is strictly<0; sign-stable resistant only if minimum>0; otherwise mixed. Very small sign-stable values still need effect-size/noise caution.
5. Compute descriptive control-vs-single-T0 guide-median differences as a baseline-fitness diagnostic, with the single-T0 limitation. Do not subtract fitness twice. Primary CTL contrast already compares treated/control. Immune-specific effect is not necessarily free of every proliferation/timing confound.
6. Join by exact gene symbols to the1,021Ramosrows (same species); record missing/nonunique symbols rather than silent case-folding or alias substitution. Retain all panel genes, exclusions and exact input values. For Ramos primary columns, define supported resistance as LFC>0 AND authors' reported FDR<=0.05 in BOTH E:T conditions. Check reported FDR meaning/direction in workbook/source if ambiguous; otherwise label inference unresolved. Do not call FDR a probability of biological truth or treat nonsignificance as no effect.
7. Primary exploratory flag: sign-stable CTL sensitization in>=3of6primarymodelcontexts AND supported Ramos resistance in both in-vitro contrasts. Predefine a symmetric descriptive resistant-CTL/sensitizing-CAR column but do not switch primary to it because it looks better. Do not search a menu of metrics or choose thresholds to get flags. Repeat the same guide/filter/sign criteria at earlier matched CTL timepoints solely as a prespecified consistency check; report instability, don't discard it to improve a claim.
8. Known comparison anchors for checking sign/meaning (not forced outcomes): B2m,Ifngr1,Ifngr2,Jak1,Jak2,Stat1,Irf1,Ptpn2,Socs1,Fitm2,Fadd,Casp8,Atg3,Atg5,Atg7,Atg12,Becn1,Rb1cc1. Respect actual eligibility/panel membership. If basic signs contradict source figures, stop and inspect implementation/mapping before believing a novel reversal. Authors' scores and this exploratory frequency estimate can differ, so report this rather than silently changing normalization.

## Verification and decision gate

- Before relying on calculations, test with a tiny synthetic matrix: identical treatment/control frequencies =>0; a depleted guide =>negative; enriched=>positive; scaling every count in a column should leave normalized frequencies nearly unchanged (pseudocount caveat); zero/missing orduplicate keys don't become eligibility; Renca-HA controls never join antigen-free cytokine controls. Test median/envelope independently with hand-calculated examples. Test only these meaningful failure modes, not an elaborate suite.
- Run helper end-to-end on actual files. Outputs finite JSON and exact TSV widths, with guide/model/panel counts and all exclusion reasons. Record no tuned parameters. No plots required at this gate; parent reviews raw records first.
- If zero primary flags, known-pattern-only findings, insufficient eligibility, unstable early/late directions, questionable calibration or only generic agreement: explicitly say this did NOT establish a substantial extension. No further fitting/access/search here.
- If a robust flag/class emerges, report exact genes and observed contexts to parent. A flag alone is not a new mechanism or publication-ready finding; existing selection, many inspected genes, score reconstruction and absent independent experimental replication limit it. Parent will review closest primary evidence and decide whether there is a genuinely worthwhile follow-up. No hit inflation, clinical prescriptions or automatic mechanism claims.
- Close out main stage: finish schemas/data cards/gates/decision, hash every scientifically relevant output including frozen prior_art/pilot files, validate new stage and priorcycle3, and update compact liveSTATUS. Record actual start/end, title-screen allocation overrun by prior-art agent, exact main query/title counts, resource/byte totals, and stage outcome. Preserve all originals/earlier outputs. Return concise results and explicit path decision by executor deadline.
