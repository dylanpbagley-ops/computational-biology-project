# Outcome uncertainty: what survives the source corrections?

October 2, 2026. Advisory analysis, not student submission prose.

**The specific new corrections do not erase the original ranking signal, but broader unresolved labels prevent a firm clinical-performance claim.** This pass uses the unchanged score in hypothetical scenarios. No new outcomes are independently adjudicated.

## Numerical results

AUC measures ranking: 0.5 is chance-level pairwise ranking and 1 is perfect ranking for the labels used. It is not patient-level accuracy or risk.

| Scenario | Frozen rule | Confinement alone | Blood flag alone |
|---|---:|---:|---:|
| Original 30 targets | 0.914 | 0.880 | 0.711 |
| Exclude FAP's mouse-only injury example | 0.911 | 0.879 | 0.711 |
| Exclude FAP; make GPRC5D unfavorable | 0.884 | 0.843 | 0.722 |
| Exclude FAP; make NECTIN4 unfavorable | 0.924 | 0.889 | 0.722 |
| Exclude FAP; make both unfavorable | 0.902 | 0.858 | 0.735 |

All other labels remain unchanged. These changes are hypothetical, not amended ground truth. GPRC5D's individual rash mechanism is not proven. Selecting the best scenario would be inappropriate. See `scenario_results.tsv`.

The ten uncertainties are FOLR1, EGFR, EPCAM, NCAM1, F3, ERBB2, CD47, GPRC5D, NECTIN4 and CLDN18. Excluding FAP, all 1,024 assignments were evaluated while holding 16 favorable and 3 unfavorable labels fixed. Those fixed labels remain assumptions.

| Conditional sensitivity analysis | Rule AUC range | Rule minus blood-only AUC |
|---|---:|---:|
| Ten disputed labels variable | 0.517–0.924 | −0.150 to +0.206 |
| CD47 and NECTIN4 unfavorable; eight others variable | 0.643–0.924 | −0.039 to +0.205 |

**These are attainable mathematical ranges, not confidence intervals.** Assignments have no probability weights. Unknown attribution is not evidence of safety; some assignments are clinically implausible or unsupported. The broad minimum assigns eight uncertain targets favorable and GPRC5D/CLDN18 unfavorable, producing 24 favorable and 5 unfavorable targets. This witness explains the bound; it is not endorsed.

The 256-assignment constrained analysis was added after viewing the broad results and is explicitly post hoc. Its source-based rationale is in `constrained_sensitivity_note.md`. Neither range includes all uncertainty in the fixed labels, study selection, exposure, antigen identity or transportability.

The rule exceeds confinement alone throughout these assignments, by 0.028–0.092 in the broad set and 0.032–0.078 in the constrained set. This does not independently validate the blood flag: most favorable blood-lineage labels remain fixed. Its advantage over the blood-only baseline is not guaranteed.

## Attribution progress

`attribution_decisions.tsv` records eight working evidence decisions. Under the proposed severe-injury criterion, CD47 anemia and NECTIN4 skin injury support unfavorable observations in their reported contexts. The former has author-attributed on-target anemia despite a DLT exemption; the latter has severe skin events with author-probable attribution. Sources: [CD47](https://pmc.ncbi.nlm.nih.gov/articles/PMC7186585/), [NECTIN4](https://pmc.ncbi.nlm.nih.gov/articles/PMC7106979/).

GPRC5D's rash DLT and skin-expression observations support plausibility, but not proof of that patient's event mechanism. The tisotumab ocular analysis describes a target-mediated mechanism as potential. Lorvotuzumab's neuropathy DLT is verified but its target-versus-payload mechanism remains unresolved in the inspected evidence. Sources: [GPRC5D](https://www.nejm.org/doi/full/10.1056/NEJMoa2204591), [F3](https://www.sciencedirect.com/science/article/abs/pii/S0090825822001202), [NCAM1](https://pmc.ncbi.nlm.nih.gov/articles/PMC4859861/).

## Negative observations and mathematical robustness

For zero events in n independent observations with perfect detection and a common binomial risk, the one-sided 95% upper limit solves (1−p)^n=0.05, giving p=1−0.05^(1/n). Illustrative upper limits are 52.7% for n=4, 14.6% for 19, 10.5% for 27, 7.2% for 40 and 3.0% for 100. Clinical studies need not meet these assumptions. These are not risk estimates for the named treatments or retrospectively chosen eligibility cutoffs. They explain why an explicit no-event report can be useful without establishing near-zero risk.

On the 29-target panel, any one arbitrary label change can reduce AUC to 0.823, and two to 0.750. These adversarial checks have no claim of clinical plausibility. Witnesses for exactly 0 through 5 changes are in `arbitrary_flip_bounds.tsv`.

For each fixed number of changes in either direction, the worst rank sum removes the highest favorable ranks and adds the lowest unfavorable ranks. Enumerating those counts gives exact bounds. Five automated tests check direction/ties, direct pairwise agreement, exhaustive small-case bounds and witnesses, invalid inputs, and the zero-event formula. The original AUC is reproduced; integrity checks retain 19 historical files and 17 files recorded in preceding source-review manifests.

## Consequence for the project

Keep testing the existing project and frozen score. These findings support continued investigation, not a declaration of clinical validation or guaranteed competitiveness. The next priority is coherent endpoint review and adequate negative observations, particularly managed blood-lineage depletion. Fitting a more complicated score to uncertain outcomes would not resolve this problem.

An amended evaluation needs its own endpoint definition, source version and uncertainty decisions. If adequate negative evidence cannot be established, retain event/observational/indeterminate categories and narrow the claim. The original retrospective association and scoped mathematical work remain usable components; no replacement binary truth has been manufactured here.
