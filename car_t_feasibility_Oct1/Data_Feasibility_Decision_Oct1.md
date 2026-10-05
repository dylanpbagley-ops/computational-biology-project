# CAR-T direction: data feasibility and decision

Advisory assessment for Dylan, October 1, 2026. This is a research-planning memo, not a student-written paper or submission. Original project files and locked predictions were preserved.

## Decision

**Do not commit to the proposed distribution-aware sustained CAR-T model yet.** The biological motivation survives scrutiny; the accessible data do not yet support the exact experiment needed to establish its added value. Proceed only with a short, targeted attempt to obtain the missing measurements and metadata. There is also a potentially useful forecasting-methods question in the AML material, but it needs additional data before it can become a defensible project.

This is a no-go for building the large model now, not proof that the broader project cannot become competitive. It changes my earlier recommendation from a promising direction to a conditional candidate. The investigation has not established good award odds, and neither a polished reproduction nor added model complexity would do so.

## What was actually checked

- Downloaded the complete B7-H3 paper and its numerical source-data ZIP, extracted and inventoried 93 text tables, and inspected the antigen, killing, expansion, serial challenge and imaging tables.
- Downloaded the AML paper, repository metadata and code archive; inspected 55 code files. Checked the latest repository version: record 12823159, modified January 12, 2025, was still the latest returned by the API at inspection time.
- Inspected the central directory of the AML 2.4GB data archive using an HTTP byte range. It contains exactly three entries: two saved posterior-chain files and LT6.csv. Extracted LT6.csv and the small retained-chain file directly from the fetched bytes, with ZIP CRC validation. The complete large archive and its full-archive checksum were not downloaded or verified.
- Read the Majzner CD19 paper and its methods. The linked supplement download returned HTML instead of a DOCX, so that supplement was not inspected. Its response is explicitly labeled as an unsuccessful download. No complete numerical event-data release was verified for this candidate.
- Checked nearby primary literature for prior work on antigen-density thresholds, continuous antigen distributions and dynamic CAR-T expansion.
- Ran an original descriptive reproduction of three B7-H3 source tables, preserving donor labels and averaging the three labeled wells within each donor before averaging donors. No predictive model was trained; no new prediction accuracy or independent validation is claimed.

Downloads have URL, retrieval time and SHA256 provenance. The downloaded Julia code was read, not executed. Scripts and results are in this folder. There was no author contact.

## Why the exact question cannot yet be tested

The proposed question was: **Do measured distributions of antigen across individual tumor cells, combined with measured CAR-T killing and expansion responses, predict sustained control better than mean antigen density or a static killing threshold?**

That question requires more than a histogram in a paper. We need numerical individual-cell values or usable distribution counts, calibration information, and a clear link to the functional experiments. We also need independent experiments for evaluation, with simple competitors trained and evaluated under the same information constraints.

| Candidate | Verified usable material | Missing or limiting information | Current role |
|---|---|---|---|
| Barisa et al., B7-H3, 2025 | 93 figure-input tables; calibrated density values for target models; donor-labeled killing, T-cell counts and residual antigen summaries | No FCS/event files in the source ZIP; density table 3a gives one value per model, not full single-cell distributions; unclear replication mapping in key tables; few donors and constructs | Useful discovery and reproduction material; conditional primary candidate if raw flow data and metadata are obtained |
| Shah et al., AML, 2026 | Source code and LT6.csv: 120 rows, 24 experimental combinations, five times through day 16; tumor/T-cell abundances | LT15 training file required by code absent from inspected archives; LT6 has no replicate identifiers or measured antigen-density column; validation evaluation needs closer scrutiny | Conditional candidate for an honest forecasting benchmark; not a distribution-based validation dataset |
| Majzner et al., CD19, 2020 | Readable paper, calibrated antigen-density methodology, killing/proliferation and animal results | No verified matched numerical release of single-cell densities and sustained functional outcomes; linked supplement not retrieved successfully | Biological reference; not currently an executable primary dataset |

The absence of raw data in an inspected release does not prove the authors lack it. B7-H3 explicitly allows requests to the corresponding author. Conversely, a public-data statement is not evidence that all required files are actually present.

The AML file's CD33 column is tumor-cell abundance, not CD33 molecules per cell. The CD3 column is T-cell abundance. Confusing these with antigen density would create a fundamentally wrong analysis.

## Numerical reproduction: what it establishes

The following values come from 6d.txt and 6e.txt. Each is an equal-weight mean over the three donor labels after averaging three wells per label. They are exploratory summaries of released files, not a claim to reproduce the paper's complete biological sample.

| Construct | Day 1 specific cytotoxicity | Day 7 specific cytotoxicity | Day 5 T-cell count / initial 10,000 |
|---|---:|---:|---:|
| MGA.271 | 29.0% | 32.9% | 5.48-fold |
| TE9 | 61.0% | 71.1% | 13.49-fold |
| 376.96 | 64.0% | 81.5% | 18.71-fold |

These values establish that the files contain measured functional differences worth investigating. They do not establish that a new distribution-aware model improves predictions. The construct ordering at day 1 already matches the ordering at day 7 in these group means. There are individual-label differences, but choosing an interesting exception after inspecting it would be discovery, not independent validation.

The paper already relates avidity, antigen density, cytokines, proliferation and sustained activity. Recalculating those relationships, fitting another regression or showing that the constructs differ would mostly reproduce its contribution.

Important constraints:

- Figure 6 describes four donors across four experimental replicates; 6d/e/i expose three donor labels with three well labels each. Their mapping must be clarified before treating the tables as complete or estimating independent sample size.
- Nine wells per construct are not nine independent donors. Time points from a culture are also dependent.
- Specific cytotoxicity is a normalized assay quantity. Negative values and values above 100 elsewhere in the tables should not be silently clipped or interpreted as cell-death probabilities.
- Residual antigen in 6i is geometric mean fluorescence intensity, not a calibrated full density distribution. A geometric mean does not determine how many cells lie below a functional threshold.
- Weekly irradiated-target assays, daily challenges and imaging after sorting/resting/resetting the effector-to-target ratio should not be concatenated as one uninterrupted biological time course.

Machine-readable calculations: `descriptive_donor_means.csv` and `reproduction_results.json`. Antigen-density source: `barisa_source/3a.txt`. Replication and calibration limitations are more consequential than the number of numerical rows.

## Novelty: what would and would not be a contribution

Antigen-dependent CAR-T function is established. Majzner showed that constructs differ in antigen sensitivity and that killing alone need not imply the proliferation needed for sustained activity. Barisa already links functional avidity to antigen thresholds and serial-challenge performance. Santurio's published 2024 model already represents continuous antigen levels and resistance mechanisms. Shah models joint killing and expansion across constructs and conditions.

Therefore, the general statements “antigen heterogeneity matters,” “include expansion,” and “use a dynamic model” cannot be the novelty claim. A useful contribution would need to establish a specific advantage under a demanding test, such as:

1. A measured distribution contributes predictive information beyond its mean, using matched functional experiments and genuinely held-out runs or donors.
2. A model trained on early observations reliably forecasts later control or construct ordering under a specified challenge protocol, outperforming simple early-readout and growth baselines.
3. An experimental-design result identifies which additional measurement or perturbation resolves an ambiguity, and that result is tested with withheld measurements or new experiments.

None of these advantages has been demonstrated in this feasibility investigation.

## A distinct AML opportunity, with a serious qualification

The archived validation routine computes the agreement of each simulated trajectory with validation outcomes, then chooses the trajectory with the largest combined distance correlation. In `aml_code/plotting_functions.jl`, this is visible in `validation_ts`, especially lines 366–431: `best_pair = argmax(target_cor + effector_cor)`, followed by plotting and reporting the selected trajectory and correlations. `plot.jl` invokes that routine for Figure 3.

Parameters may have been trained on separate data, but selecting a particular prediction using later outcomes means that selected trajectory is not an outcome-independent forecast. This measures the best match among sampled trajectories rather than the performance of a prespecified central forecast. It can give an optimistic impression of predictive ability. The magnitude of any optimism has not been quantified here.

The published paper describes maximal-correlation trajectories for its best-fit examples and validation using posterior trajectories. This warrants checking the final figure generation, rather than assuming the older archive exactly reproduces every final reported result. There is also a mismatch between the paper's description of posterior predictive uncertainty including measurement error and the inspected routine's visible trajectory/standard-deviation calculation. Do not conclude a final-paper error or misconduct from this inspection alone.

A possible new question is: **How reliably do early CAR-T/tumor abundance measurements predict later control when forecasts are fixed before later outcomes are observed, and which experimental conditions cause predictive failure?**

Use all appropriate posterior draws for a prespecified forecast distribution; evaluate absolute/log-scale count errors and interval calibration, not correlation alone. A favorable correlation does not establish correct counts or correct decisions. Compare with simple baselines using exactly the same initial measurements. Treat outcome-selected trajectories only as an explicitly labeled oracle comparator.

This could be scientifically useful if it produces a generalizable measurement or modeling lesson. A single corrected plotting choice or one paper-specific criticism would be a much weaker project. We cannot perform the required full reproduction now because the training file is absent from the inspected releases. Saved posteriors might allow a partial audit, but they are not a substitute for a reproducible training dataset. The 120 validation rows also contain no replicate column, so uncertainty about biological replication remains.

## Next steps and stopping rules

**Immediate priority: secure one complete experimental dataset and metadata before substantial model development.** Keep the existing project available as a bounded fallback while testing this opportunity.

For B7-H3, the specific missing items are:

- Numerical raw antigen flow files for the relevant initial and residual tumor populations, with negative controls and calibration beads where applicable.
- Matching donor, run, target, construct and well identifiers; exact challenge timing, input cell numbers, sampling/removal and replenishment rules.
- Raw functional readouts and controls where possible; clarification of the Figure 6 donor/replicate discrepancy.
- Enough independent runs to reserve evaluation data. Data already read in this investigation are discovery data, even if they later receive a train/test label.

For AML, the specific missing items are:

- LT15 training data and its replicate identities, plus confirmation of how LT6 values were summarized.
- Final code corresponding to the published Figure 3 and uncertainty calculations.
- If available, a distinct donor/run or cell-line dataset for a meaningful transfer test.

A data request would need to be sent by Dylan or explicitly authorized separately. A mentor can help assess whether the request is realistic and scientifically useful; specialized lab access is not a prerequisite for that discussion.

Use a roughly one-week access check, then a short pilot if the required files arrive. This is a proposed planning limit, not a scheduled task or a promise that authors will reply within a week. The pilot must compare a prespecified question against simple baselines and report how much independent evidence is available. Stop or narrow if distribution data are unavailable, matching is ambiguous, or evaluation is only within the same few observed constructs. Do not keep adding parameters to compensate for missing measurements.

If neither dataset becomes complete, I recommend comparing genuinely different questions in an established mentor system or in another fully accessible computational dataset. School-lab experiments should be considered only after its capabilities and the relevant assay are confirmed. An unrelated demonstration would not validate CAR-T efficacy or EV capture.

## Implication for the ISEF goal

The proposed pivot has not earned commitment yet. With currently verified files, it supports careful reproduction and exploratory analysis, but lacks the evidence needed for the stronger predictive contribution we were aiming for. Obtaining suitable raw data and a meaningful independent test could change that assessment substantially.

I do not conclude that the existing project is incapable of reaching a high level under any changes. That sweeping conclusion would require ruling out many directions we have not tested. I do conclude that ambition, available time, additional equations and presentation quality cannot resolve the present data and novelty limits. The next valuable step is an evidence-based access decision, followed by a small decisive pilot or an early change of question.

## Primary sources

- Barisa et al. [Functional avidity of anti-B7H3 CAR-T constructs predicts antigen density thresholds for triggering effector function](https://www.nature.com/articles/s41467-025-61427-4). Article and linked numerical source-data ZIP downloaded.
- Majzner et al. [Tuning the Antigen Density Requirement for CAR T-cell Activity](https://pmc.ncbi.nlm.nih.gov/articles/PMC7939454/). Paper downloaded; linked supplement not obtained.
- Santurio et al. [Mechanisms of resistance to CAR-T cell immunotherapy: Insights from a mathematical model](https://www.sciencedirect.com/science/article/abs/pii/S0307904X23003773). Published 2024; DOI 10.1016/j.apm.2023.08.029. Used for prior-art scope, not data extraction.
- Shah et al. [Quantifying target antigen-dependent CAR T-cell performance against AML](https://pmc.ncbi.nlm.nih.gov/articles/PMC13420862/). Paper downloaded; DOI 10.1186/s12935-026-04419-8.
- [AML code/data repository, Zenodo 12823159](https://zenodo.org/records/12823159). Latest-version API checked; downloaded code and inspected data archive entries.
