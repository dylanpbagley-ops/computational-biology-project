# Feasibility findings — October 2, 2026

Advisory research-planning material, not student-authored submission text. These are exploratory pilots. Protocols were saved before numerical comparisons; results do not establish a publishable method or competition prospects.

**Decision: stop developing the tested option-1 candidate; retain option 2 for a bounded, more demanding benchmark. Neither pilot yet demonstrates a project significantly stronger than the existing work.** Option 2 passes the access and measurement-problem checks, which is enough to keep investigating it rather than move immediately to another direction.

## Option 1: small protein panels across donors

Tested whether selecting measurements by their weakest training-donor separation improves recovery in a new donor. Compared that simple rule with standardized mean differences and average donor-level AUC, using the same classifier and measurement budgets. This screens a concrete candidate; it does not test every possible uncertainty or measurement-selection method.

The source contains 228 protein measurements and 161,764 cells from eight donors. Donors, including all their time points, stayed together: six training donors, one calibration donor, and one test donor in each of eight rotations. Marker selection used training data only. The primary task was recovery of the annotated dendritic-cell group; a B-cell task supplied an easier reference. Annotation recovery is the endpoint, not clinical prediction or physical cell sorting.

The primary comparison used four measurements and a nominal training budget of 1,000 cells per donor, supplemented where needed to include up to 50 positive examples and reweighted for sampling. Threshold calibration targeted 1% false-positive rate, which is not equivalent to 1% contamination.

| Prespecified primary criterion | Required | Observed |
|---|---:|---:|
| Median paired recall improvement over mean-difference selection | At least 5 percentage points | 0.60 percentage points |
| Donors with improved recall | At least 6 of 8 | 4 of 8 |
| Test donors with false-positive rate at most 1.5% | At least 6 of 8 | 5 of 8 |
| Median paired improvement over average-AUC selection | Evidence beyond that simple baseline | 0.00 percentage points; wins in 2 of 8 |

The candidate failed the screening criteria. The easy B-cell task was already near saturation: median recall approximately 99.5% for both the conventional and proposed approach. I did not search alternative endpoints or tune the failed candidate after seeing these results. A full comparison with established panel-selection packages was not run because this first candidate failed its preliminary screen.

The labels were produced using joint RNA/protein analysis, so they are not independent ground truth for a protein predictor. Eight donors remain eight donor-level units; the large cell count does not fix that limitation. These facts make a small gain insufficient grounds to commit to this direction.

Data checks confirmed barcode and feature order and inspected 3,358 raw coordinates across the source matrix. One historical cache entry was missing because the loader skipped three lines despite a two-line header. The pilot restored that one raw value in memory under a recorded addendum; the historical cache and source were preserved. This is a reproducibility correction, not a reason to reinterpret the negative screen as positive.

Sources: [GEO source sample](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5008738), [reference annotation documentation](https://azimuth.hubmapconsortium.org/references/). Reproduction files: `measurement_feasibility_Oct2/precommitment.md`, `precommitment_addendum_1.md`, `data_audit.json`, `pilot.py`, `results.csv`, and `results.json`.

## Option 2: how much imaging is enough for a reliable measurement?

Question: can distribution-aware uncertainty tell us when another image field is needed, improving reliability at a fixed acquisition budget? The reusable principle is that a population average can conceal heterogeneity and sampling uncertainty. The proposed contribution must be a useful, validated measurement decision—not simply using distributions.

The public JUMP Cell Painting pilot provides individual-object measurements with field, well, and plate identities. I retrieved a fixed selection of six compound conditions and six solvent-control wells on each of two same-batch plates, nine fields per well. All 216 expected files were present. The extracted eight features covered 21,692 nuclear objects, with 45–188 per field; all selected values were finite. Object identifiers and extract hashes were checked. Raw-image segmentation quality has not been independently reviewed.

For each well, I compared its full nine-field distance from controls with estimates from three fields. An optimistic computational comparator sampled the same number of cells scattered throughout the full well. Two summaries were fixed in advance: feature means, and feature 10th/50th/90th percentiles. Thirty resampling draws per well describe sensitivity; they are not independent biological replicates.

| Summary | Median absolute error: three fields | Same-size scattered-cell sample | Error ratio |
|---|---:|---:|---:|
| Means | 0.02772 | 0.01588 | 1.746 |
| Quantiles | 0.02720 | 0.01643 | 1.655 |

Errors are in standardized distance-to-control units. Thus whole-field sampling had approximately 75% and 65% greater median error respectively. Both exceeded the prespecified 20% flag for investigating field structure. This is consistent with field-level heterogeneity or artifacts; it does not establish their cause. Scattered-cell sampling is a computational benchmark, not a deployable acquisition method.

**The important negative result:** full-profile rank agreement between the two plates was weak: Spearman correlation 0.168 for means and 0.042 for quantiles across the 12 matched wells. Those wells include six controls; these correlations are descriptive, not estimates of compound-prediction accuracy. The pilot gives no convincing advantage for quantile profiles, and the two profile types have different dimensions. Slightly smaller error on one summary is not evidence of a superior biological representation.

There are only six compound conditions and two same-batch plates. Full nine-field measurements are same-well reference values, not independent truth. Controls and scaling use full reference data, so this pilot is not a complete low-budget acquisition simulation. An adaptive stopping policy has not been implemented or validated. The experiment establishes that access and retrospective sampling analysis are workable, and that ignoring acquisition structure can be optimistic.

Sources: [Cell Painting Gallery](https://github.com/broadinstitute/cellpainting-gallery), [JUMP pilot repository](https://github.com/jump-cellpainting/2024_Chandrasekaran_NatureMethods_CPJUMP1), [data structure](https://broadinstitute.github.io/cellpainting-gallery/data_structure.html). Exact source URLs, retrieval metadata, and hashes are retained in `imaging_feasibility_Oct2/download_manifest.json` and per-file provenance records. Protocol, fixed selection, retained extracts, sampling implementation, detailed results, and summary are in that directory.

## What would justify committing to option 2?

The next study should answer one practical question: **Can a rule identify when a measurement needs another field, and save imaging while retaining reliable calls on independent replicate experiments?**

Before further outcomes are inspected, fix a larger development/test split by whole plates or batches, with enough distinct conditions and replicates to support the endpoint. First establish repeatable signal after documented image/object quality checks. If full-budget measurements themselves do not replicate, repair the measurement endpoint or stop; an adaptive policy cannot rescue an unreliable reference.

Then compare fixed field budgets, a simple variance-based rule, and the proposed distribution-aware rule at matched acquisition cost and error tolerance. Hold all policy tuning in development data. Evaluate the whole stopping policy on held-out plates, include the control-imaging budget, and retain uncertain/no-call results in the comparison. Any practical success threshold must be fixed before that evaluation. Public retrospective experiments can test measurement reliability; actual microscope-time savings would still need appropriate acquisition validation.

This is a recommendation for the next bounded test, not a claim that a novel algorithm has been found. The literature already includes distribution profiles and imaging-budget studies; see the sources in `direction_landscape_Oct1/Stronger_Project_Options_Oct1.md`. A successful contribution would need an advantage over those relevant methods or a well-validated failure boundary that changes measurement decisions.

For now, keep the current project as a usable fallback and retain its mathematics, data-processing experience, and distributional reasoning. Do not discard it or invest months in the imaging pivot on the strength of this pilot alone. Option 3 remains a possible mathematical component rather than a demonstrated independent project; option 4 still depends on suitable experimental access. Neither dependency needs reopening while option 2's next test remains viable.
