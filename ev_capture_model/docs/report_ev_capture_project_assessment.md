# Independent Red-Team Assessment of the EV-Capture Project

**Assessment date:** 26 August 2026  
**Audience:** solo high-school ISEF project  
**Scope:** scientific novelty, mathematical correctness, physical applicability, data support, external validation, reproducibility, and next-step viability

## Executive verdict

The project is **not currently defensible as a calibrated predictor of real extracellular-vesicle (EV) capture yield or selectivity**. Its central parameters (`K_INTRA`, `z`, and `k_geom`) are uncalibrated and partly confounded; the EV–antibody contact regime is not derived from the physical system; most Lennon point clouds are reconstructed; pooled vesicles lack donor labels; and no untouched experimental dataset validates quantitative yield.

The project **is defensible after reframing** as:

> **A conditional, distribution-aware design framework and measurement-gap analysis for EV affinity capture.**

That version asks a valuable question: *how can per-EV marker heterogeneity, EV size, contact-patch assumptions, and reconstruction choices change predicted capture rankings?* The strongest surviving empirical example is the **CA19-9 plasma comparison**, whose ranking change remains under strict provenance filtering. The analogous **EGFR inversion does not survive** removal of density-estimated and split-cluster points and should be retired.

The broad novelty claim that “population integration is the core mathematical novelty” is contradicted by prior work on heterogeneous multivalent particles and ligand mixtures [20, 96, 97, 105]. A narrower, search-bounded contribution may remain: no direct paper was found in the 32-paper assessment set that combines measured joint per-EV marker-count and size distributions with an equilibrium EV-capture ranking sensitivity analysis. That is **not proof of priority** and should not be phrased as “first ever.”

## What was independently checked

- The uploaded project was preserved unchanged.
- The unmodified test suite was rerun: **39 tests passed in 3.39 seconds**.
- All 13 CSVs were profiled for dimensions, ranges, provenance, and anomalies.
- The selectivity derivative was derived analytically and checked against an independent stable finite difference.
- The current exponential partition factor was compared with a fixed-total-site alternative.
- Lennon rankings were stress-tested over marker, provenance, patch exponent, two physical formulations, eight `K_INTRA` values, source weighting, gate shifts, and coordinate jitter.
- The stem-cell crossover was tested across 63 reconstructed scenarios and eight `K_INTRA` values.
- **32 papers** were included in the novelty assessment; **16 pivotal papers were reviewed from locally retrieved full text**. Abstract-only records were not used for detailed quantitative claims.

## Six layers of evidence

| Layer | Verdict | Meaning |
| --- | --- | --- |
| 1. Software correctness | Strong | The implementation is reproducible under its own tests. |
| 2. Mathematical correctness | Mixed but repairable | The equilibrium equation is coherent; exact alpha replaces the coarse finite difference. |
| 3. Physical applicability | Weak/uncertain | The chosen ensemble and contact-patch scaling are not uniquely justified for mobile EV receptors on immobilized antibodies. |
| 4. Data support | Conditional | CA19-9 has one robust assumption-sensitive pattern; EGFR does not. Data are reconstructed and pooled. |
| 5. External validation | Insufficient | Current comparisons are internal, tuned, qualitative, or use a different physical system. |
| 6. Novelty | Narrow | The EV-specific integration may be useful; the general mathematics and heterogeneity concepts are established. |

## Reconstructed project state and data limitations

The project uses

```text
q_b = exp(n sigma K_INTRA) - 1
theta = z q_b / (1 + z q_b)
n_contact = k_geom n diameter^(p - 2)
```

where `p=0` represents a fixed-absolute-patch/density-driven endpoint and `p=2` a fixed-fraction-patch/size-blind endpoint. These endpoints are separately normalized and should be presented as **sensitivity assumptions**, not as calibrated versions of one experiment.

The four pooled Lennon Figure 3c reconstructions contain 2,050 EGFR PDAC, 815 EGFR healthy, 2,653 CA19-9 PDAC, and 644 CA19-9 healthy vesicle points. However, density-estimated points make up 77.0%, 88.8%, 82.0%, and 93.8% of those tables, respectively. The strict individually resolved subsets are much smaller: EGFR PDAC 148, EGFR healthy 37, CA19-9 PDAC 129, and CA19-9 healthy 12.

The source paper used plasma from **6 healthy and 5 PDAC donors**, normalized each sample to `7×10^8 EVs/150 µL` before qSMLM, removed EVs with radius below 7 nm or fewer than 2 molecules, and reported simulation-based underestimation for fewer than 15 molecules per EV and large sparse EVs; diameters below 40 nm could be overestimated [39]. The plotted reconstruction has no donor labels. Therefore, vesicles must not be treated as independent biological replicates.

## Mathematical audit

### Exact selectivity

For `x=n sigma K_INTRA`, the exact logarithmic selectivity is

```text
alpha(x,z) = x exp(x) / [(exp(x)-1) (1 + z(exp(x)-1))].
```

The independent stable finite-difference check had maximum absolute error **1.09e-08**. On the project’s 104-cell grid, the exact peak alpha is **10.465**. The project coarse finite-difference result differs from the exact alpha by as much as **0.702** at a grid point.

The often-repeated **59.6%** value is simply 62 of 104 chosen cells with exact peak `alpha>1.01`. It is not a probability that real EV capture is superselective. Figure 1 replaces that arbitrary percentage with the continuous analytical `alpha=1` boundary.

![Exact selectivity phase map](ev_capture_alpha_phase_map.png)

### Identifiability

The output depends on combinations such as

```text
n × sigma × K_INTRA × k_geom × diameter^(p-2).
```

Consequently, `sigma` and `K_INTRA` cannot be estimated separately from one ligand-density condition. At `p=2`, `k_geom` and `K_INTRA` are confounded. In dilute regimes, `z`, geometry, and binding strength become even more confounded. Sweeping these values does not estimate them.

### Physical regime

The project expression `exp(n a)-1` is the weak-per-partner approximation to the fixed-total-site form `(1+a)^n-1`. In a controlled system with mobile DNA linkers on both colloid and flat surface, Linne et al. still used a fixed total ligand number per colloid and treated mobility through recruitment and configurational free energy [30]. Smith and Seifert modeled a vesicle with a fixed total number of mobile ligands, immobilized substrate receptors, explicit contact-zone allocation, and bending energy [19]. Therefore, receptor mobility does **not** by itself prove that the project’s exponential count ensemble is correct.

For total exponent `x≤1`, the two simplified formulations differ in predicted theta by as much as **41.8%** if low valency is allowed, and **7.25%** for `n≥10`. They happen to give similar Lennon ranking directions, but that local agreement does not settle the EV physics.

## Lennon robustness audit

| Marker | Formulation | Strict p=0 PDAC/healthy | Strict p=2 PDAC/healthy |
| --- | --- | --- | --- |
| CA199 | exponential | 0.219–0.407 | 1.262–1.447 |
| CA199 | fixed-n | 0.219–0.407 | 1.262–1.447 |
| EGFR | exponential | 0.265–0.378 | 0.959–0.990 |
| EGFR | fixed-n | 0.265–0.379 | 0.959–0.990 |

**CA19-9:** Every strict-resolved p=0 result is below 1 and every strict-resolved p=2 result is above 1. The critical crossover is approximately `p=1.58–1.62`. This is the strongest surviving conditional result.

**EGFR:** The hybrid reconstruction crosses near `p≈1.76–1.84`, but strict-resolved and non-density-estimated p=2 results stay below 1. The EGFR inversion is driven by reconstruction choices and is not robust.

Coordinate jitter of approximately ±1 nm diameter and ±0.5 copy caused narrow intervals relative to the large p=0 versus p=2 separation in the hybrid reconstruction. That test addresses digitized coordinate precision only; it does not address donor sampling, missing low-copy EVs, affinity-enrichment bias, or reconstruction-model uncertainty.

![Lennon patch-exponent robustness](ev_capture_lennon_patch_robustness.png)

Gate sensitivity reinforces the provenance result. In strict-resolved points, EGFR raw-copy ratios range 0.879–0.909 while density ratios range 0.067–0.084. CA19-9 raw-copy ratios range 1.360–1.423 while density ratios range 0.052–0.064. Thus, the lower reconstructed PDAC marker density is robust in direction, while raw copy-number conclusions depend on marker.

## Stem-cell extension

Across 63 combinations of positivity endpoint, detection threshold (4, 7, or 10 molecules), count family, overdispersion, and a 500-copy cap, CD81 exceeded CD63 at the largest tested `K_INTRA` in **100%** of scenarios. Yet this occurred at only one of eight `K_INTRA` values, or **12.5%** of all scenario-by-parameter rows.

This is a useful **conditional mechanism demonstration**, not an independent biological discovery. It is reconstructed from summary statistics and pooled positivity ranges; the fractions are not tissue-specific, the result is not disease-related, and the affinity endpoint is not externally calibrated.

## External validation audit

- **Monte Carlo:** internal consistency only. It re-samples the same assumed Poisson process.
- **DNA nanostars:** related multivalent behavior, but tuned and physically different; not EV validation [33].
- **Mizenko:** strong qualitative evidence that capture-marker choice biases apparent EV phenotypes [1]. It does not provide verified 10×/5× capture-efficiency calibration; the paper’s approximately 10× statement refers to SP-IRIS detecting fewer EVs than fluorescence.
- **Chiodi:** closest EV-specific modeling prior art. It measured real-time capture for three antibodies at six surface concentrations and fit fast/slow binding components [112]. It does not validate the present joint size/count equilibrium model, but it shows what an informative calibration experiment looks like.

## Novelty verdict by component

| Component | Verdict | Defensible description |
| --- | --- | --- |
| Biological problem | Not novel but important | EV capture markers bias which heterogeneous vesicles are observed [1,5]. |
| Use of multivalency/superselectivity | Established | Core concepts and equations have extensive theory and controlled experiments [30,31,38,48]. |
| Analytical selectivity derivation | Useful correction, not broad novelty | Exact alpha is a clean project-specific derivation of an established model class. |
| Integrating heterogeneous particle populations | Established generally | Prior nanoparticle and ligand-mixture models already integrate valency distributions [20,96,97]. |
| Joint per-EV size/count integration for EV capture | Possibly narrow/search-bounded novelty | No direct match was found in the assessed set; avoid universal priority language. |
| Provenance-sensitive reversal analysis | Potentially distinctive | The CA19-9 versus EGFR contrast and reconstruction dependence are a credible methodological finding. |
| Measurement-gap analysis | Useful synthesis | Explicit identifiability and calibration requirements are actionable. |

## Viability scorecard

Scores use 1 = currently inadequate and 5 = strong.

| Criterion | Score / 5 | Reason |
| --- | --- | --- |
| Scientific premise | 4 | Capture-marker bias and per-EV heterogeneity are real and important. |
| Equation-to-system applicability | 2 | A useful simplification, but the ensemble/contact mechanics are not derived for EV-on-antibody capture. |
| Data adequacy | 2 | Rich vesicle counts, but mostly reconstructed, pooled, affinity-enriched, and donor labels are absent. |
| Parameter identifiability | 1 | K_INTRA, z, k_geom, and some geometry effects are uncalibrated and confounded. |
| External validation | 1 | No untouched EV capture-yield prediction has been tested. |
| Numerical reliability | 3 | Core code is stable and tested, but alpha should be analytical rather than coarse finite difference. |
| Reproducibility | 4 | The unmodified suite passes 39 tests and the audit is scripted. |
| Scope control | 2 | Current handoff overreaches; a narrower framework is defensible. |
| ISEF competitiveness now | 2 | Interesting mechanism and strong computational work, but current predictive/novelty claims are too broad. |
| ISEF competitiveness after repair | 4 | Could become a strong project with one external test and disciplined claim language. |

### Viable now

1. The software is internally reproducible.
2. Exact analytical alpha and the `alpha=1` phase boundary are correct for the stated equation.
3. Per-EV distributions can change conditional population-average rankings.
4. CA19-9 provides a provenance-robust example of sensitivity to patch scaling.
5. Stem-cell CD81/CD63 provides a conditional endpoint mechanism example.
6. The project clearly identifies which measurements are missing.

### Viable after repair

1. Replace finite-difference alpha with the exact stable expression.
2. Present exponential and fixed-total-site formulations as alternatives.
3. Treat `p` as a continuous sensitivity axis and remove calibrated-geometry language.
4. Center the project on CA19-9; retain EGFR as a negative control showing reconstruction sensitivity.
5. Rewrite novelty as EV-specific application/synthesis, not invention of population integration.
6. Add one untouched external test or collaborator experiment.

### Not viable on current evidence

1. Quantitative real-world capture-yield prediction.
2. A robust EGFR plasma ranking inversion.
3. Donor-level or clinical claims from pooled vesicles.
4. “First-ever population integration” language.
5. Nanostar or Monte Carlo analyses as external EV validation.
6. A marker-specific claim that EV receptor mobility has been measured.

## Claim audit summary

| Claim | Status | Allowed wording |
| --- | --- | --- |
| C01 | software-only | The archived implementation is internally reproducible under its current tests. |
| C02 | supported | For the project’s stated equilibrium equation, alpha has the exact closed form reported in the audit. |
| C04 | contradicted | On the project’s chosen 104-cell grid, 59.6% of cells exceeded alpha=1.01; this is not a prevalence estimate. |
| C06 | unsupported | The exponential form is one simplifying regime; fixed-total-site and contact-zone models are credible alternatives that should be compared. |
| C09 | contradicted | Monte Carlo reproduces the assumed stochastic model and is an internal consistency check. |
| C11 | conditionally supported | In this reconstructed dataset, the conditional CA19-9 ranking changes with the assumed size scaling and survives strict provenance filtering. |
| C12 | contradicted | The EGFR inversion is reconstruction-sensitive and should not be claimed as robust. |
| C15 | conditionally supported | A CD81-over-CD63 crossover is a robust conditional endpoint demonstration across the tested reconstructed families. |
| C18 | unsupported | The model provides conditional rankings and sensitivity maps, not calibrated capture yield. |
| C19 | contradicted | The mathematical ingredients are established; the possible novelty is the EV-specific combination and diagnostic sensitivity analysis. |
| C20 | conditionally supported | In the literature searched through 26 August 2026, no direct EV study was found that combines measured joint per-EV size/count distributions with equilibrium capture-ranking sensitivity. |
| C23 | conditionally supported | A conditional framework showing how heterogeneity, size scaling, and reconstruction provenance affect capture rankings. |

The full 24-row audit, including failure modes and required evidence, is in `ev_capture_claim_audit.csv`.

## Six-month ISEF roadmap

### Weeks 1–2: repair the scientific story

- Freeze the current archive and retain the passing test baseline.
- Implement exact alpha with stable `expm1`, `log1p`, and overflow guards.
- Replace “59.6% plausible space” with the analytical phase boundary.
- Rename the two binding formulations descriptively rather than calling either one “correct.”
- Rewrite all conclusions using the allowed wording in the claim-audit CSV.

**Go/no-go:** if the project narrative still requires calibrated yield or a universal inversion, stop and reframe before adding more computation.

### Weeks 3–5: make uncertainty the result

- Use CA19-9 as the primary positive demonstration.
- Use EGFR as a negative control: the inversion disappears under strict provenance filtering.
- Report hybrid, non-density-estimated, strict-resolved, and density-only results side by side.
- Avoid vesicle-level p-values because donor labels are missing.
- Keep the stem-cell analysis in a short mechanism section.

**Go/no-go:** if CA19-9 loses its direction under a newly justified physical model, pivot fully to “why current data cannot determine rankings.” That remains a legitimate methodological result.

### Weeks 5–8: obtain one real external test

Priority order:

1. Ask Lennon study authors for donor-labeled diameter/copy tables and reconstruction details.
2. Ask Mizenko/ExoView authors for single-EV event tables and exact capture-surface conditions.
3. Ask Chiodi/IRIS authors for probe-density–resolved binding curves or processed tables.
4. Seek an untouched single-EV dataset with size and marker counts [43,44].

Pre-register the prediction before seeing the answer: marker, direction, p range, formulation set, and pass/fail criterion.

### Optional minimal collaborator experiment

Use one EV source and one marker first. A practical minimal design is:

- **EV marker abundance:** low and high source/preparation, independently quantified if possible.
- **Antibody surface density:** at least 3 levels spanning weak to strong capture.
- **EV concentration:** at least 2 levels.
- **Time:** kinetic readout or at least 3 incubation times.
- **Controls:** isotype/non-binding surface, monovalent or Fab-like control if feasible, and a capture-independent particle count.
- **Replication:** 3 independently prepared surfaces per condition; biological preparation replication if available.
- **Readouts:** bound EV count, size distribution, post-capture marker distribution, and input distribution.

This design begins to separate surface density, effective affinity, transport/time, and population heterogeneity. It is much more informative than adding another unconstrained parameter sweep.

### Final months

- Reproduce from a clean environment.
- Lock claim language and figures.
- Document AI use, data provenance, digitization, and all exclusions.
- Build the poster around one positive result (CA19-9), one negative control (EGFR), and one measurement-gap diagram.

## Recommended one-sentence project claim

> **Using reconstructed single-EV size and marker-count distributions, I developed a conditional multivalent-capture framework and showed that predicted target rankings can depend on contact-patch scaling and data reconstruction provenance; the analysis identifies the measurements required before quantitative EV-capture prediction is possible.**

## Bottom line

This is a **promising computational research project after disciplined narrowing**, not a failed project. Its strongest contribution is not a universal target recommendation. It is the demonstration that EV heterogeneity, size scaling, provenance, and non-identifiability can change what a capture model appears to conclude—and a concrete plan for collecting the missing evidence.

## Assessed literature set

The paper-level evidence depth, overlap, and novelty implication are provided in `ev_capture_novelty_evidence.csv`.

- 1. Mizenko et al. (2021). Tetraspanins are unevenly distributed across single extracellular vesicles and bias sensitivity to multiplexed cancer biomarkers. Journal of Nanobiotechnology. DOI: 10.1186/s12951-021-00987-1. Evidence: full text.
- 3. Jiang et al. (2024). Multiparametric profiling of HER2-enriched extracellular vesicles in breast cancer using Single Extracellular VEsicle Nanoscopy. Journal of Nanobiotechnology. DOI: 10.1186/s12951-024-02858-x. Evidence: full text.
- 5. Zarovni et al. (2025). Stoichiometric constraints for detection of EV-borne biomarkers in blood. Journal of Extracellular Vesicles. DOI: 10.1002/jev2.70034. Evidence: full text.
- 10. Liu et al. (2021). Rapid Capture of Cancer Extracellular Vesicles by Lipid Patch Microarrays. Advanced Materials. DOI: 10.1002/adma.202008493. Evidence: abstract/metadata.
- 11. Jahnke et al. (2024). Membranes on the move: The functional role of the extracellular vesicle membrane for contact-dependent cellular signalling. Journal of Extracellular Vesicles. DOI: 10.1002/jev2.12436. Evidence: full text.
- 12. Zaborowski et al. (2019). Methods for Systematic Identification of Membrane Proteins for Specific Capture of Cancer-Derived Extracellular Vesicles. Cell Reports. DOI: 10.1016/j.celrep.2019.03.003. Evidence: abstract/metadata.
- 13. Kucuk et al. (2026). Thermally Modulated Metasurface Sensor for Dynamic and Time-Resolved Isolation of Extracellular Vesicles. Advanced Materials. DOI: 10.1002/adma.202522964. Evidence: abstract/metadata.
- 19. Smith et al. (2005). Effective adhesion strength of specifically bound vesicles. Physical Review E. DOI: 10.1103/physreve.71.061902. Evidence: full text.
- 20. Li et al. (2014). Evaluating Binding Avidities of Populations of Heterogeneous Multivalent Ligand-Functionalized Nanoparticles. ACS Nano. DOI: 10.1021/nn406455s. Evidence: abstract/metadata.
- 24. Merminod et al. (2021). Avidity and surface mobility in multivalent ligand–receptor binding. Nanoscale. DOI: 10.1039/d1nr02083h. Evidence: full text.
- 30. Linne et al. (2021). Direct visualization of superselective colloid-surface binding mediated by multivalent interactions. Proceedings of the National Academy of Sciences. DOI: 10.1073/pnas.2106036118. Evidence: full text.
- 31. Martinez-Veracoechea et al. (2011). Designing super selectivity in multivalent nano-particle binding. Proceedings of the National Academy of Sciences. DOI: 10.1073/pnas.1105351108. Evidence: abstract/metadata.
- 32. Xia et al. (2023). The Role of Receptor Uniformity in Multivalent Binding. JACS Au. DOI: 10.1021/jacsau.3c00052. Evidence: abstract/metadata.
- 33. Linne et al. (2024). Optimality and cooperativity in superselective surface binding by multivalent DNA nanostars. Soft Matter. DOI: 10.1039/d4sm00704b. Evidence: full text.
- 38. Dubacheva et al. (2023). Determinants of Superselectivity─Practical Concepts for Application in Biology and Medicine. Accounts of Chemical Research. DOI: 10.1021/acs.accounts.2c00672. Evidence: full text.
- 39. Lennon et al. (2019). Single molecule characterization of individual extracellular vesicles from pancreatic cancer. Journal of Extracellular Vesicles. DOI: 10.1080/20013078.2019.1685634. Evidence: full text.
- 43. Ko et al. (2021). Sequencing-Based Protein Analysis of Single Extracellular Vesicles. ACS Nano. DOI: 10.1021/acsnano.1c00782. Evidence: full text.
- 44. Saftics et al. (2023). Single Extracellular VEsicle Nanoscopy. Journal of Extracellular Vesicles. DOI: 10.1002/jev2.12346. Evidence: full text.
- 48. Dubacheva et al. (2019). Multivalent Recognition at Fluid Surfaces: The Interplay of Receptor Clustering and Superselectivity. Journal of the American Chemical Society. DOI: 10.1021/jacs.8b12553. Evidence: full text.
- 49. Curk et al. (2017). Optimal multivalent targeting of membranes with many distinct receptors. Proceedings of the National Academy of Sciences. DOI: 10.1073/pnas.1704226114. Evidence: abstract/metadata.
- 53. Morzy et al. (2022). Significance of Receptor Mobility in Multivalent Binding on Lipid Membranes. Angewandte Chemie International Edition. DOI: 10.1002/anie.202114167. Evidence: full text.
- 94. Liu et al. (2020). Combinatorial entropy behaviour leads to range selective binding in ligand-receptor interactions. Nature Communications. DOI: 10.1038/s41467-020-18603-5. Evidence: abstract/metadata.
- 95. Tito et al. (2016). Communication: Simple approach for calculating the binding free energy of a multivalent particle. The Journal of Chemical Physics. DOI: 10.1063/1.4948257. Evidence: abstract/metadata.
- 96. Waddell et al. (2010). Origin of broad polydispersion in functionalized dendrimers and its effects on cancer-cell binding affinity. Physical Review E. DOI: 10.1103/physreve.82.036108. Evidence: full text.
- 97. Tan et al. (2021). A general model of multivalent binding with ligands of heterotypic subunits and multiple surface receptors. bioRxiv (Cold Spring Harbor Laboratory). DOI: 10.1016/j.mbs.2021.108714. Evidence: abstract/metadata.
- 98. Angioletti-Uberti (2017). Exploiting Receptor Competition to Enhance Nanoparticle Binding Selectivity. Physical Review Letters. DOI: 10.1103/physrevlett.118.068001. Evidence: abstract/metadata.
- 105. Woythe et al. (2021). A quantitative view on multivalent nanomedicine targeting. Advanced Drug Delivery Reviews. DOI: 10.1016/j.addr.2020.11.010. Evidence: abstract/metadata.
- 106. Feng et al. (2024). High-throughput capture and in situ protein analysis of extracellular vesicles by chemical probe-based array. Nature Protocols. DOI: 10.1038/s41596-024-01082-z. Evidence: abstract/metadata.
- 109. Hendrix et al. (2023). Extracellular vesicle analysis. Nature Reviews Methods Primers. DOI: 10.1038/s43586-023-00240-z. Evidence: abstract/metadata.
- 112. Chiodi et al. (2021). Multiplexed Affinity Measurements of Extracellular Vesicles Binding Kinetics. Sensors. DOI: 10.3390/s21082634. Evidence: full text.
- 173. Morales et al. (2022). Future of Digital Assays to Resolve Clinical Heterogeneity of Single Extracellular Vesicles. ACS Nano. DOI: 10.1021/acsnano.2c04337. Evidence: abstract/metadata.
- 175. Carney et al. (2024). Harnessing extracellular vesicle heterogeneity for diagnostic and therapeutic applications. Nature Nanotechnology. DOI: 10.1038/s41565-024-01774-3. Evidence: abstract/metadata.
