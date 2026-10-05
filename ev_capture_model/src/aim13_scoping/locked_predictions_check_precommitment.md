# Locked-predictions check: precommitment (written Sept 29, 2026, before any search)

Object: `prospective_predictions.json` (generated 2026-09-08). 52 targets: 14 TOLERATED, 18 ON-TARGET TOXICITY, 20 no call.
The predictions are NOT re-scored, re-derived or edited. The rule that produced them never saw outcomes.

Labelling standard (identical to `Validation_Targets_Sept29.md`): safe = approval or documented clinical tolerability, no
dose-limiting on-target normal-tissue toxicity. toxic = documented severe, dose-limiting or fatal toxicity attributed by the
authors (or, for an approved drug, a boxed warning stating the tissue) to the target on normal tissue. Excluded: non-protein
antigens, early-only evidence, same gene with an approval and a toxicity record, blockers rather than cell-killers.
Each label needs a primary source (paper, FDA label or company statement) with a date. Hedged attribution (payload,
lymphodepletion, off-target) = reported as HEDGED, never scored.

Evidence date classes (fixed now)
 A. FORWARD: the verdict evidence first became public after 2026-09-08 (the lock date). Genuinely prospective.
 B. PRE-LOCK: the evidence was public before 2026-09-08 but was not used in the prediction (the rule uses only normal-tissue
    data). Valid as an out-of-sample check of the rule, NOT as a time-stamped forward test. Reported separately.
 C. UNRESOLVED: no verdict meeting the standard. Expected for most targets.
Only the 32 called targets (14 + 18) are scored. No calls are not scored (a separate informational note may list any
no-call target that has a clear verdict).

Scoring: hit = TOLERATED and safe, or ON-TARGET TOXICITY and toxic; miss = the opposite; hedged/unresolved = not scored.
Report: counts by class A/B/C, every hit and miss with source and date, both classes shown separately and combined,
and the same table restricted to targets NOT already in the main panel or the validation set (F3, AXL, ERBB3 are in the main
panel as toxic / efficacy_fail / efficacy_fail; they are reported but flagged as not independent).
Expectation stated now: only three weeks separate the lock from today, so class A will be empty or nearly empty; class B is
the realistic yield; most targets will be C. No conclusion is drawn from fewer than 8 scored targets beyond listing them.
Misses are reported exactly as prominently as hits. No label enters the frozen panel from this check.
