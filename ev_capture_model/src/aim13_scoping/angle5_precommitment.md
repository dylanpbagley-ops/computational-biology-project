# Angle 5 (Cas9 variants) -- analysis plan, written BEFORE any outcome was computed

Written 2026-09-28, after the data audit and before any test below was run.
Data (SHA-256):
  25ca0328d4b543aa33c18afb9b8b6918feb07e4ae0ce5ac74f280b9e5a57dbd2  crisprSQL_100720.zip
  5291917cc72eea15031dde55a8a48759ff9dff8e7b87e2cc55b30e815745c105  gkad702_Supplemental_Files.zip

## What the audit changed (facts, not outcomes)
* NAR 2023 Table S10 has 306 guides x up to 7 DESIGNED mismatch sites (three mismatch positions
  per guide, all 1-, 2-, 3-mismatch combinations): 1,924 pairs. It is not a genome-wide natural
  off-target set and gives at most 7 points per guide. Within-guide frontiers (the original
  plan) are therefore NOT feasible. Every test below pools across guides.
* crisprSQL: 25,632 sites, 87 guides, 17 studies; only 27.6% of sites carry a gene ID, so
  essentiality weighting would rest on a minority of sites. Essentiality weighting is dropped
  from the primary analysis for that reason and because DepMap could not be downloaded.

## The question
The ceiling theorem's premise is that a therapy's kill probability is a MONOTONE function of one
latent binding score. If that holds for Cas9 variants, then a variant changes stringency (a
threshold shift) but cannot re-rank sites, so every site's variant rate is one monotone function
of its wild-type rate, and on-target and off-target sites lie on the SAME curve. If it fails,
the variant changes discrimination itself (it moves the curve), which the theorem does not cover.

## Tests (pooled; guide-grouped 5-fold cross-validation so no guide is in train and test)
T1  On-target, one map. Isotonic (monotone) regression of HiFi on-rate and LZ3 on-rate on WT
    on-rate across the 306 guides. Report cross-validated R^2 and Spearman rho.
T2  Off-target, one map. Same regression across all 1,924 off-target pairs, overall and within
    mismatch class (M1 / M2 / M3).
T3  Same curve for on- and off-target? Fit the monotone map on OFF-target pairs only; predict
    each guide's ON-target variant rate from its WT on-rate; restrict to the range of WT rates
    where both site types occur (5th-95th percentile of the overlap). Report the median signed
    residual (observed minus predicted, in rate units) with a guide-level bootstrap 95% CI.
    Positive = the on-target site retains more activity than off-target sites did at the same
    WT rate, i.e. the variant DISCRIMINATES beyond a pure stringency shift.

## Predictions (my honest priors, stated so they can be wrong)
P1  Cross-validated R^2 for T1 will be moderate, between 0.5 and 0.85 (the source paper's
    headline is guide-specific efficiency loss, which is scatter beyond a single map).
P2  T2 will show a tighter relationship than T1 (off-target rates are dominated by mismatch
    count, which every variant respects).
P3  T3 median residual will be POSITIVE with a CI excluding 0: high-fidelity variants improve
    discrimination and do not merely shift a threshold.

## Interpretation rules, fixed now
* R^2 >= 0.90 in T1: premise strongly supported and the paper's 'guide-specific' claim is
  mostly WT-level effects. R^2 < 0.40: premise fails for on-target rates.
* T3 CI includes 0: variants behave as a pure stringency shift; the ceiling applies cleanly.
  T3 CI above 0: variants move the curve; the write-up must say the ceiling bounds a FIXED
  score and that engineering a new score is outside it. This is a limitation of the transfer,
  not a confirmation, and must never be written as validation.
* Whatever happens is reported, in this section's order, with these predictions alongside.

## Decision rule for the angle's fate
* Angle 5 stays ONE bounded section (the plan's own limit). It is written up ONLY as 'a test of
  the theorem's premise in a second domain' -- never as a second clinical validation.
* If T1 and T3 both fail the premise, the section becomes a documented limitation and the angle
  is otherwise dropped rather than stretched.
