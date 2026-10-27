# Pre-registration: Target-ranking validation against clinical ground truth

---

## Purpose and scope

This document registers the labeling rules, analysis plan, and decision rules for the
target-ranking validation study.

---

## 1. Question

Does ranking candidate immunotherapy target antigens by the overlap of their tumor and normal
per-cell surface-density distributions predict which antigens have produced approved drugs,
better than ranking them by mean expression on tumor cells?

## 2. Hypotheses

- **H1.** Antigens with an approved therapy will rank higher by the overlap measure than
 antigens known to be broadly expressed on normal cells.
- **H2.** The mean-expression ranking will separate those two groups less well than the
 overlap ranking.
- **H3 (cross-disease).** The same pattern will hold in a second disease with different
 antigens and a different measurement technology.
- **H4 (boundary condition).** Adding further normal-cell compartments will lower ceilings
 and never raise them, so antigens excluded on a narrow compartment set stay excluded.

## 3. Labeling rules

An antigen is labeled **VALIDATED** if, for the disease in question, a therapy directed at
that antigen has full regulatory approval (any modality: CAR-T, bispecific engager, monoclonal
antibody, antibody-drug conjugate).

An antigen is labeled **BROAD / EXPECTED-RISK** if it is documented as expressed across
multiple normal leukocyte lineages or on all nucleated cells, independent of this project's
data.

An antigen is labeled **TRIAL** if a directed therapy has entered human trials without
approval, and **UNLABELED** otherwise. Only VALIDATED and BROAD enter the primary test.

Labels are assigned from published clinical status only, never from any number this project
computed. The full label list is fixed in section 6 and will not be revised after seeing
results; any antigen whose status changes will be reported as a deviation.

## 4. Primary analysis

For each antigen, compute the probability that a randomly chosen tumor cell carries more of it
than a randomly chosen normal cell, separately for each normal compartment, and take the
minimum across compartments as that antigen's ceiling. Rank antigens by ceiling. Rank the
same antigens by median tumor mean expression.

**Primary outcome:** the fraction of (VALIDATED, BROAD) pairs in which the VALIDATED antigen
ranks higher, computed for both rankings. H1/H2 are supported if the overlap ranking's value
is higher than the mean ranking's.

**Significance:** for perfect separation of *v* validated above *b* broad antigens, the
probability under random ranking is 1 / C(v+b, v), reported exactly.

## 5. Falsification criteria

- The overlap ranking separating VALIDATED from BROAD no better than the mean ranking.
- Approved targets clustering at the bottom of the overlap ranking.
- The second disease showing the opposite ordering.
- Any ceiling rising when a normal compartment is added (this would break the theorem).

These will be reported whichever way they come out.

## 6. Fixed label lists

**Disease 1 — diffuse large B-cell lymphoma** (mass cytometry, 37 antigens)

VALIDATED: CD19 (axicabtagene/tisagenlecleucel/lisocabtagene, blinatumomab); CD20 (rituximab,
mosunetuzumab, glofitamab, epcoritamab); CD22 (inotuzumab ozogamicin); CD79B (polatuzumab
vedotin).

BROAD / EXPECTED-RISK: CD44, HLA-ABC, CD48, CD27, CD43, CD38, CD84, HLA-E, CD25, FasR.

*Deviation noted in advance:* CD45RA was included in the BROAD list in the exploratory
analysis and is removed here, because naive B cells are CD45RA-high, so it is not cleanly a
pan-leukocyte marker. Both versions will be reported.

**Disease 2 — multiple myeloma** (CITE-seq, 63 antigens)

VALIDATED: CD269/BCMA (idecabtagene, ciltacabtagene, teclistamab, elranatamab); CD38
(daratumumab, isatuximab); CD319/SLAMF7 (elotuzumab).

TRIAL: CD307/FcRH5 (cevostamab). Reported separately, not in the primary test.

BROAD / EXPECTED-RISK: CD45, CD44, CD3E, CD4, CD8A, CD14, CD11B, CD11C, CD16, CD33, CD56.

NEGATIVE CONTROLS: FITC, PE (panel fluorophores with no biological target; must sit at
chance).

**Prior-failure case, registered as a stated expectation:** CD138 (indatuximab ravtansine,
discontinued) is predicted to score well on bone-marrow data and to be a false positive,
because its dose-limiting expression is on epithelium, which bone marrow does not contain.
This is registered in advance as a demonstration of the framework's stated boundary.

## 7. Analyses completed prior to registration

The DLBCL analysis in section 6's first list was completed before this registration and gave
perfect separation (1.000) for the overlap ranking versus 0.705 for the mean ranking. This
result is not pre-registered. The registration governs: the healthy-PBMC comparator analysis,
the myeloma analysis under these exact labels, and any further disease.

## 8. Deviations

Any departure from this plan will be listed in the final write-up with the reason, including
analyses that were run and did not work.

---

# Addendum 1 (September 7, 2026) — the whole-body extension

**Registration status.** The analyses described in this addendum were completed before the
date of registration; their results are not pre-registered and must not be presented as such.
This addendum freezes the label set, the measures, and the decision rules for all subsequent
analyses, and provides a dated record of the label definitions.

## A1. Frozen definitions

**Label file.** `target_labels.py`, exported verbatim to
`target_labels_frozen.tsv` and to `docs/Label_Review_Sheet.md`. 44 targets:
19 `safe`, 11 `toxic`, 12 `efficacy_fail`, 2 stated exclusions. Every label is clinical safety
history, carries a one-line citation, and no analysis script may write to that file.

**Measures**, both parameter-free:
`top_share` = share of a gene's total normal expression in its single highest cell type;
`n_types` = number of cell types above 10% of that maximum.

**Expendability factor.** A home cell type counts as expendable if and only if it is a
haematopoietic lineage, matched by the word-boundary pattern list in `target_labels.py`. This
is assigned from the cell type's identity alone and never from a drug outcome.

**Two-factor rule.** `score = top_share + 1.0 if the home is haematopoietic`. No fitted
parameter, and the `+1` dominates by construction, so the rule is exactly "expendable first,
confinement second".

**Test.** Rank separation (equivalently AUC) between the `safe` and `toxic` groups, with a
one-sided permutation p-value; exact enumeration where C(n,k) ≤ 250,000 and 200,000 random
permutations otherwise, with the result stating which was used.

## A2. Decision rules for all subsequent runs

- **H5.** On a protein-level normal reference, `top_share` will separate `safe` from `toxic`
 with separation > 0.70 and p < 0.05. Below that, the RNA-level result is reported as
 reference-dependent.
- **H6.** The two-factor rule will beat the one-factor rule on any new reference. If it does
 not, the expendability factor is reported as not replicating.
- **H7.** On shuffled outcome labels, separation will fall to 0.5 ± 0.1. If it does not, the
 test is broken and every result depending on it is withdrawn.
- Any target added to the panel after this date is recorded with the date it was added, and
 results are reported both with and without the additions.
- A label may be changed only by an outside reviewer, using the review sheet, and any change
 is reported as a deviation with the before-and-after numbers.

## A3. Deviations recorded prior to registration

One errors was identified and is recorded in the results document:

1. **Haematopoietic substring bug** (results §23.2). A plain substring test classified
 "conjunctival goblet cells" as haematopoietic because the string contains `"t cell"`. It
 inflated the two-factor separation from the correct 0.905 to 0.916. Matching is now
 word-boundary anchored with an assertion-based guard table.

## A4. Falsification criteria, whole-body extension

- Confinement failing to separate `safe` from `toxic` on a protein-level reference.
- The result disappearing when blood-cancer targets are removed. (Already tested: it does
 not — 0.792, exact p = 0.023 on solid tumours alone.)
- An outside reviewer's corrected labels dropping the separation below significance.
- The efficacy-failure control separating as strongly as the toxicity comparison, which would
 mean the measure tracks "programme failed" rather than "on-target toxicity".

## A5. Outcomes of the registered protein-level hypotheses (recorded September 7, 2026)

**H5 — split result.** On the immunohistochemistry reference (76 tissue × cell-type groups):
separation 0.785, p = 0.0104 — met (0.775, p = 0.0158 with uncertain antibodies dropped).
On the mass-spectrometry reference (20 bulk tissues): separation 0.650, p = 0.1148 — not met.

Two post-hoc explanations for the mass-spectrometry failure were proposed and both were
tested and rejected: sparse measurement (separation falls to 0.393 as sparsely-measured
targets are removed) and insufficient granularity (randomly coarsening the 154-cell-type RNA
reference to 20 groups still yields 0.837, with 0.650 outside the range of all 60 random
partitions).

Per the registered rule, the section 23 result is reported as reference-dependent, and the
claim is scoped to restriction measured across normal cell types rather than stated
universally. Results document §26.

**H6 — not testable.** Both protein references are indexed by tissue or tissue × cell type,
so no target's home cell type can be a haematopoietic lineage and the two-factor rule
collapses onto the one-factor rule by construction. Recorded as untestable on these
references, not as a failure of the hypothesis.

**H7 — met, on the two-antigen machinery rather than the whole-body test.** With outcome
labels shuffled, the single-antigen separation collapses from 0.692 to 0.002. Results
document §25.

---

# Addendum 2 (October 2026) — locked-prediction resolution protocol

**Registration status.** The 52 prospective predictions governed by this protocol were locked
on September 9, 2026 (`prospective_predictions.json`, filed in this OSF project). No
prediction has resolved as of the date of this addendum. This addendum specifies the
resolution rules in advance of any resolution. It registers no results.

## B1. Scope

The 52 locked predictions in `prospective_predictions.json` (September 9, 2026), each pairing
a cell-surface therapeutic target with a predicted clinical tolerability outcome (TOXIC or
TOLERATED), under the frozen label set and two-factor rule defined in Addendum 1.

## B2. Resolution criteria

A locked prediction resolves when a pre-specified evidentiary bar is met:

- A TOXIC prediction is **confirmed** by a human trial report of dose-limiting or severe
  on-target toxicity attributed to the target (not the payload or modality class).
- A TOLERATED prediction is **confirmed** by regulatory approval, or by a phase-complete
  trial with no dose-limiting on-target toxicity.
- A TOXIC prediction is **refuted** when a phase-complete trial with pre-specified safety
  reporting shows no dose-limiting on-target toxicity attributed to the target. The same
  evidence, symmetrically, refutes a TOLERATED prediction when a trial report documents
  dose-limiting on-target toxicity attributed to the target.
- **Censoring.** A program discontinued for efficacy before adequate safety exposure leaves
  the prediction unresolved, with the reason recorded as "censored — program ended before
  toxicity could be assessed."

Attribution follows the endpoint contract from pass 5 (managed depletion / severe
consequences / supportive care / recovery recorded as separate fields).

## B3. Authoritative sources

Trial registries, primary publications, and FDA/EMA labels and review documents.
Preclinical-only evidence does not resolve a prediction.

## B4. Adjudication

Resolutions are extracted by someone blind to the prediction where possible. Ambiguous cases
are recorded as unresolved, never forced.

## B5. No-call handling

No-calls are excluded from scoring and reported as a coverage limitation.

## B6. Minimum-resolutions rule

No formal conclusion is drawn until at least 8 predictions have resolved. Below that,
results are reported as descriptive only.

## B7. Primary metric

Rank separation (AUC) between resolved TOLERATED and resolved TOXIC predictions, with the
exact permutation p-value. Secondary: hit rate with a binomial interval.

## B8. Schedule

Formal resolution checks occur quarterly; each check is logged with its date and sources.
Interim program monitoring (monthly, plus targeted checks after major conference
proceedings) is logged but confers no resolution authority; only a quarterly check under
§B2–§B4 resolves a prediction.

## B9. Monitoring log

| Date | Event | Effect on predictions |
|---|---|---|
| 2026-09-25 | CD276 / I-DXd BLA withdrawn by the sponsor (efficacy/regulatory grounds; no on-target toxicity finding) | CD276 prediction unaffected; recorded as a censored-type monitoring event — program ended before toxicity could be assessed |
| 2026-10-25 | Scheduled re-check of 4 monitored programs (CD276, CD70/cusatuzumab, CDH6/R-DXd, VTCN1/B7-H4): no resolutions; all other programs active | All 52 predictions remain open |

---
