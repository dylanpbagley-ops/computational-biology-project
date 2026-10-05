# Stronger project directions: research-planning assessment

October 1, 2026. Advisory material for Dylan, not student-authored submission text. No historical scoring rules, labels, predictions, or project results were changed. This is a candidate-selection assessment, not evidence that any proposed method already works or is novel.

**Recommendation.** First investigate distribution-aware measurement design: determining when a small protein-marker panel can reliably distinguish a cell population in a new donor, when to decline to make a call, and which additional measurement would most improve the decision. Keep reliable drug-screen analysis as the alternative. Both have a clearer route to empirical evaluation than extending the current clinical-toxicity ranking. EV capture-bias correction remains an attractive experimental option if a suitable measurement system becomes available.

The objective is not to combine these into a larger project. A short comparison should select ONE question for sustained work. An award probability cannot be inferred from an appealing question; independent results and the student's understanding remain decisive.

**What transfers from the existing project.** Its reusable assets are distribution overlap, nonlinear response averaging, joint marker information, sensitivity to assumptions, attainable selectivity limits, and experience processing single-cell measurements. The clinical safety labels are not required for these new questions. The old ROC mathematics may inform an assay trade-off, but assay identification is not equivalent to therapy action. Arbitrary classifiers, monotone responses, and physically realizable gates have different constraints.

The October handoff supersedes the September headline. The current project is a retrospective normal-expression safety association with important endpoint and validation limitations; the earlier EV and CAR-T work supplies methods and motivation. These proposed directions do not depend on defending its clinical interpretation.

**Landscape screened.**

| Direction | Potential contribution | Present assessment |
|---|---|---|
| Small marker panels under donor variation | Determine when a panel is trustworthy and which extra measurement is useful | Best first feasibility pilot; direct reuse and local data |
| Drug-screen measurement reliability | Reproducible detection of heterogeneous cell responses at a limited imaging budget | Strong alternative; ample public infrastructure, substantial prior work |
| Missing joint-marker measurements | Bound what can be concluded and select informative co-measurements | Useful mathematical component; classical bounds alone are not enough |
| EV capture-bias correction | Recover or bound population composition distorted by selective capture | Strong experimental story, currently conditional on matched measurements |
| Cancer-drug response prediction | Predict an independent functional endpoint from pretreatment heterogeneity | Potentially valuable, but matched data and patient count are limiting |
| Sustained CAR-T prediction | Forecast later functional performance from early measurements | Prior local audit found essential missing data; do not restart model building |
| Larger normal-tissue safety score | More reliable agent-specific safety prediction | Retains the current project's hardest label and confounding problems |
| Other heterogeneous systems, such as adsorption or environmental exposure | Transfer nonlinear averaging to a measurable physical endpoint | Possible, but this search established no specific advantage over the closer options |

This is a bounded landscape search, not an exhaustive enumeration of applications or a proof of novelty.

**1. Measurement design for small marker panels — recommended first pilot.**

Question: Given a fixed number of measured protein markers, can distribution overlap and donor variation identify when a cell population can be recovered reliably, and whether another marker or more cells would improve the result?

For example, two panels might distinguish a cell population equally well in pooled training data. One may fail in a particular donor because its desired and undesired cell distributions overlap there. A useful method would detect that limitation before deployment, select a more stable panel, or explicitly decline a confident answer. This is a proposal, not an observed finding.

Verified local starting material: the GSE164378 protein matrix header has 228 rows and 161,764 cell columns. Its barcodes match all metadata cell identifiers in exact order. Metadata contains eight donors, sample/batch information, and three cell-label resolutions. These checks establish a usable structure, not correct marker naming, independent labels, or completed preprocessing. The existing loader hard-codes marker names; source feature-order verification is required before analysis.

Eight donors are eight donor-level units, not 161,764 independent validation subjects. Some annotated populations have very few total cells, making them unsuitable as the initial endpoint. The cohort has been used in earlier work, so it is development material; a new analysis split must not be described as a never-seen prospective cohort.

The proposed contribution is a validated reliability and measurement-selection rule, not merely selecting markers. COMET already proposes combinatorial marker panels; GateFinder already finds simple gates; scPanel already uses patient-level splitting for sparse gene panels; CellMarkerPipe already varies marker count, population imbalance, and cell count. ADTnorm addresses protein-measurement shifts. These are nearby methods or essential context, not all interchangeable direct competitors.

Sources: [COMET](https://pmc.ncbi.nlm.nih.gov/articles/PMC6811728/), [GateFinder](https://pmc.ncbi.nlm.nih.gov/articles/PMC6247943/), [scPanel implementation](https://github.com/carissaxie/scPanel), [CellMarkerPipe](https://www.nature.com/articles/s41598-024-63492-z), [ADTnorm](https://www.nature.com/articles/s41467-025-61023-6).

Pilot design: choose a prespecified, adequately represented population; audit label provenance; compare a conventional marker set, a simple classifier, an applicable established selection method, and the proposed overlap-based selection. Keep panel size and training information comparable. Hold out whole donors with all their time points together. Fit selection and normalization decisions inside training partitions. Evaluate recovery of desired cells at a fixed contamination constraint, separately by donor, plus uncertainty calibration and the fraction of samples for which a call is made.

Do not obtain apparently reliable performance simply by declining almost every sample. Compare methods at the same accepted-call rate or show the complete error-versus-call-rate curve. More cells reduce sampling uncertainty; they do not necessarily resolve biological overlap. A new marker can add distinguishing information, but redundant markers may add little. The scientific result would be a tested rule for distinguishing those situations.

Main validity risks: annotations may derive partly from the same proteins being evaluated; withholding donors alone does not remove that circularity. Seek an orthogonal labeling route and a separate cohort. CITE-seq antibody-derived counts are not calibrated antigen molecule counts and do not validate a physical flow-cytometry panel without a transfer experiment. Artificially hiding columns measures an information limitation, not reagent interference, staining changes, or real assay cost. Report a marker-count/cell-count trade-off unless actual cost information is obtained.

Go forward only if the pilot shows a specific, repeatable benefit over strong baselines or a reproducible failure boundary that improves measurement choices. Stop if the result reduces to rediscovering familiar lineage markers, if labels cannot support the claim, or if existing methods solve the selected task equally well. A mathematical impossibility statement from an estimated distribution must acknowledge estimation uncertainty; finite data cannot prove that no better classifier exists without assumptions.

**2. Reliable drug-screen decisions from heterogeneous cell images — strongest alternative.**

Question: Can a method decide whether a drug-induced change is reproducible, and whether more imaging is needed, without losing signals confined to a minority of cells?

This transfers the distribution principle from antigen abundance to cell morphology: size, shape, texture, and other measured features. The practical outcome would be fewer incorrect screen decisions at a fixed measurement budget, or comparable reliability using fewer images. Neither outcome would by itself establish clinical drug efficacy.

Cell Painting is not an empty field. A 2019 study combined average, dispersion and covariance information; a 2022 framework uses distribution distances; SPACe uses Earth Mover's Distance. Published imaging-system and assay-optimization studies already examine cell counts and numbers of imaged sites. Therefore, “use the full distribution,” “use Wasserstein distance,” and “measure more cells” are not sufficient novelty claims.

Sources: [heterogeneity-aware profiling](https://www.nature.com/articles/s41467-019-10154-8), [distribution framework](https://www.nature.com/articles/s42003-022-04343-3), [SPACe](https://www.nature.com/articles/s41467-024-54264-4), [imaging-system study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10841730/), [assay optimization](https://pubmed.ncbi.nlm.nih.gov/37344608/).

The narrower candidate is an uncertainty-aware decision or stopping rule: after initial fields have been observed, decide whether a response is sufficiently reproducible to call or whether another field is needed. Test on real replicate wells and held-out plates; use blocks of fields rather than pretending individual cells are independent acquisition units. Repeated looks require a valid sequential procedure or evaluation of the entire fixed stopping policy on held-out experiments.

Compare mean/median profiles, mean-plus-dispersion, an established distribution method, and fixed imaging budgets. Include cell count, plate location, segmentation errors, and plate/lab effects as possible explanations. Agreement with all images from the same well is only a fidelity check; independent wells and known perturbation relationships provide a stronger endpoint. A retrospective replay is an acquisition simulation, not a prospective microscope experiment.

The [Cell Painting Gallery](https://github.com/broadinstitute/cellpainting-gallery) and [JUMP pilot repository](https://github.com/jump-cellpainting/2024_Chandrasekaran_NatureMethods_CPJUMP1) provide documented public resources. Repository availability was checked; a complete site-level numerical subset and its joins were not downloaded and validated in this turn. That is the first access gate. Start with processed measurements from a small subset, not the entire image collection. If only well averages are available for a chosen subset, that subset cannot test the proposed distribution method.

Advance only if an advantage persists against existing distribution-based methods on independent wells/plates. A gain over means alone, or a gain created only by synthetic mixtures, is insufficient for the intended stronger claim.

**3. Measurement value when joint distributions are missing — a possible mathematical core.**

Question: When separate measurements of markers leave several possible joint distributions, can we bound the range of achievable identification performance and choose the co-measurement that most reduces decision uncertainty?

This directly extends the project's existing coexpression limitation into a measurement question. Start with real jointly measured protein data, conceal selected joint information, derive allowed outcome ranges, then reveal the withheld data to test the usefulness of those ranges. Compare the recommended added measurement with random selection and simple correlation-based selection on development donors, followed by external evaluation.

Classical probability bounds and linear programming are established tools; applying them does not constitute inventing a theorem. The new evidence would need to be that the bounds meaningfully flag unreliable conclusions and guide useful measurements. The earlier project already studied identifiability, so repeating its impossibility result is not a new direction. This is best considered within option 1 rather than added as a second project.

**4. EV capture bias as the object of study — conditional experimental option.**

Question: How much does capture distort the apparent composition of an EV mixture, and can independently calibrated capture measurements recover or bound its original composition?

This preserves the most of the original EV work. However, selective capture bias is already demonstrated in [Mizenko et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC8379740/). A stronger result would require a correction or bound tested against known input composition and independent preparations, not another demonstration that capture markers differ.

Necessary evidence includes independently characterized input and recovered material, usable marker measurements, and controls for measurement bias. A population that is never captured cannot generally be reconstructed from captured material alone without extra assumptions or measurements. The project's uncalibrated EV parameters must not be silently treated as measured capture probabilities.

Do not reopen the retired image archive or closed datasets as though they were new leads. No suitable lab arrangement is confirmed. A bead or other model-particle experiment could validate a physical selection model, but not automatically tumor-EV recovery or cancer diagnosis. This becomes competitive with the computational options only when an actual measurement system and a distinct experimental test are available.

**Why not prioritize cancer-response forecasting now?**

Predicting an independent measured response is more valuable than rearranging modeled target rankings, but public summaries are not automatically matched individual-cell/outcome data. The prior CAR-T feasibility audit documented missing inputs. For the AML pharmacoscopy lead, the paper provides processed single-cell releases, while raw sequencing access is restricted; the locally inspected release metadata lists roughly 5.2 GB of CyTOF and 307 MB of RNA data. Neither file size nor a public-data statement establishes the necessary patient/assay matching. No complete prediction dataset was validated here. [AML source and availability statement](https://www.nature.com/articles/s41467-024-53535-4).

**Decision sequence.** Allow roughly two weeks for feasibility, not a commitment to a large model. First establish trustworthy feature names, labels and donor grouping for option 1, reproduce a relevant baseline, and define the distinct measurement question. In the alternative track, inspect one small Cell Painting subset for individual-cell, field, well and replicate identifiers. Select one route based on actual access and baseline difficulty. Before outcome comparisons, write and hash a protocol with the endpoint, comparisons, held-out units and stopping rules. Discovery results remain discovery results. Seek external-cohort confirmation for the selected route and leave time for Dylan's own analysis interpretation and writing.

What would make the final project substantially stronger is a coherent chain: a concrete measurement problem, a justified method, a fair comparison, a result on independent units, and a bounded practical claim. This investigation identifies plausible ways to obtain that chain. It has not yet demonstrated the performance improvement or novelty needed to claim it is achieved.
