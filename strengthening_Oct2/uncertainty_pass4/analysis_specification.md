# Outcome uncertainty analysis specification

October 2, 2026. Written before this pass's numerical execution, with full knowledge of earlier scores, labels and diagnostics. This is a post-discovery sensitivity analysis, not preregistration or independent validation.

Question: how sensitive is the existing target-ranking result to label assignments that its clinical source review has left unsettled?

Keep the original score (top_share + M2_flag), confinement (top_share), and blood-flag baseline unchanged. Positive labels in the AUC calculation mean the historical favorable/safe category; higher scores predict that category. No new training or selection of scoring parameters.

1. Reproduce the original 30-target AUC, then exclude FAP alone because its cited injury evidence is preclinical. Treat this exclusion as a human-evidence scope diagnostic, not a complete amended panel.
2. On the 29-target panel, change GPRC5D alone, NECTIN4 alone, and both from favorable to unfavorable. These are source-motivated hypothetical scenarios. They do not assert that protocol DLT, severe toxicity and causally established normal-tissue injury are interchangeable. Hold all other historical labels unchanged.
3. Enumerate every binary assignment for ten identified uncertainties: FOLR1, EGFR, EPCAM, NCAM1, F3, ERBB2, CD47, GPRC5D, NECTIN4, CLDN18. Keep FAP excluded and the other 19 labels fixed. The ten cases span uncertain attribution, mismatched endpoints, context, and antigen identity, not a probability model. Report the full attainable AUC range, witnesses, and the paired differences from both baselines. The 16 favorable and 3 unfavorable fixed labels remain assumptions, not verified ground truth; this envelope is conditional and does not include all clinical uncertainty.
4. Compute worst/best AUC after exactly k arbitrary label changes for k=0,...,5 on that same 29-target panel. These are adversarial mathematical robustness checks, not medically justified relabelings. No significance testing or probabilities of assignments.
5. Explain why zero observed events in a small study does not establish zero risk. For an ideal independent binomial sample with perfect detection and a common risk, derive the one-sided 95% upper limit after zero events as 1-0.05^(1/n). Use illustrative n=4,19,27,40,100. This is an illustration, not calibrated clinical risk estimation or a retrospective eligibility cutoff.

Do not convert any output into a final amended clinical label or overwrite frozen artifacts. A broad range means existing evidence cannot identify the metric tightly under the stated assumptions; it does not mean all assignments are equally plausible. A favorable sensitivity result would not replace independent endpoint review.
