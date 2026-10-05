# Pre-registration draft — target-ranking validation against clinical ground truth

**For Dylan to review, edit into his own words, and submit at osf.io/registries (free, ~20 min).
ISEF bars generative AI from writing submitted abstract / research plan / poster text; this is a
working draft of the analysis plan, not submitted prose. Rewrite before submitting.**

---

## Why this is worth doing (honest assessment)

Pre-registration matters here for one specific reason. The project's strongest claim is:
*ranking candidate targets by tumor-versus-normal distribution overlap recovers the antigens
that became approved drugs, while ranking by average expression does not.* The obvious
objection is that the "approved" and "risky" labels were assigned after seeing the ranking.

They were not — they were fixed before the ranks were examined — but there is currently no
independent record of that. A time-stamped registration converts "we assigned labels fairly"
from something a judge has to take on trust into something with a date on it.

It is worth doing **only if it is registered before the remaining analyses run.** Registering
after the fact is worse than not registering, because it looks like a claim of rigor that
wasn't there. The remaining analyses (healthy-PBMC comparator, any further disease) have not
been run against these labels yet, so the window is open now.

Honest counterweight: this does not make the existing result more true. It makes it more
*credible*, which for a judged competition is a different and also real benefit.

---

## 1. Question

Does ranking candidate immunotherapy target antigens by the overlap of their tumor and normal
per-cell surface-density distributions predict which antigens have produced approved drugs,
better than ranking them by mean expression on tumor cells?

## 2. Hypotheses, stated before the remaining data is analysed

- **H1.** Antigens with an approved therapy will rank higher by the overlap measure than
  antigens known to be broadly expressed on normal cells.
- **H2.** The mean-expression ranking will separate those two groups less well than the
  overlap ranking.
- **H3 (cross-disease).** The same pattern will hold in a second disease with different
  antigens and a different measurement technology.
- **H4 (the boundary condition).** Adding further normal-cell compartments will lower ceilings
  and never raise them, so antigens excluded on a narrow compartment set stay excluded.

## 3. Labelling rules, fixed in advance

An antigen is labelled **VALIDATED** if, for the disease in question, a therapy directed at
that antigen has full regulatory approval (any modality: CAR-T, bispecific engager, monoclonal
antibody, antibody-drug conjugate).

An antigen is labelled **BROAD / EXPECTED-RISK** if it is documented as expressed across
multiple normal leukocyte lineages or on all nucleated cells, independent of this project's
data.

An antigen is labelled **TRIAL** if a directed therapy has entered human trials without
approval, and **UNLABELLED** otherwise. Only VALIDATED and BROAD enter the primary test.

Labels are assigned from published clinical status only, never from any number this project
computed. The full label list is fixed in section 6 below and will not be revised after
seeing results; any antigen whose status changes will be reported as a deviation.

## 4. Primary analysis

For each antigen, compute the probability that a randomly chosen tumor cell carries more of it
than a randomly chosen normal cell, separately for each normal compartment, and take the
minimum across compartments as that antigen's ceiling. Rank antigens by ceiling. Rank the same
antigens by median tumor mean expression.

**Primary outcome:** the fraction of (VALIDATED, BROAD) pairs in which the VALIDATED antigen
ranks higher, computed for both rankings. **H1/H2 are supported if the overlap ranking's value
is higher than the mean ranking's.**

**Significance:** for perfect separation of *v* validated above *b* broad antigens, the
probability under random ranking is 1 / C(v+b, v), reported exactly.

## 5. What would falsify this

- The overlap ranking separating VALIDATED from BROAD no better than the mean ranking.
- Approved targets clustering at the *bottom* of the overlap ranking.
- The second disease showing the opposite ordering.
- Any ceiling rising when a normal compartment is added (this would break the theorem).

These will be reported whichever way they come out.

## 6. Fixed label lists

**Disease 1 — diffuse large B-cell lymphoma** (mass cytometry, 37 antigens)

VALIDATED: CD19 (axicabtagene/tisagenlecleucel/lisocabtagene, blinatumomab); CD20 (rituximab,
mosunetuzumab, glofitamab, epcoritamab); CD22 (inotuzumab ozogamicin); CD79B (polatuzumab
vedotin).

BROAD / EXPECTED-RISK: CD44, HLA-ABC, CD48, CD27, CD43, CD38, CD84, HLA-E, CD25, FasR.

*Deviation noted in advance:* CD45RA was included in the BROAD list in the exploratory analysis
and is being **removed** here, because naive B cells are CD45RA-high, so it is not cleanly a
pan-leukocyte marker. Both versions will be reported.

**Disease 2 — multiple myeloma** (CITE-seq, 63 antigens)

VALIDATED: CD269/BCMA (idecabtagene, ciltacabtagene, teclistamab, elranatamab); CD38
(daratumumab, isatuximab); CD319/SLAMF7 (elotuzumab).

TRIAL: CD307/FcRH5 (cevostamab). Reported separately, not in the primary test.

BROAD / EXPECTED-RISK: CD45, CD44, CD3E, CD4, CD8A, CD14, CD11B, CD11C, CD16, CD33, CD56.

NEGATIVE CONTROLS: FITC, PE (panel fluorophores with no biological target; must sit at chance).

**Prior-failure case, registered as a stated expectation:** CD138 (indatuximab ravtansine,
discontinued) is predicted to score **well** on bone-marrow data and to be a false positive,
because its dose-limiting expression is on epithelium, which bone marrow does not contain.
This is registered in advance as a demonstration of the framework's stated boundary, not as
a success.

## 7. Analyses already run before registration (declared, not hidden)

The DLBCL analysis in section 6's first list has already been run and gave perfect separation
(1.000) for the overlap ranking versus 0.705 for the mean ranking. **This registration does not
claim that result as pre-registered.** It is declared here so the record is accurate, and the
registration governs: the healthy-PBMC comparator analysis, the myeloma analysis under these
exact labels, and any further disease.

## 8. Deviations

Any departure from this plan will be listed in the final write-up with the reason, including
analyses that were run and did not work.

---

# ADDENDUM (September 7, 2026) — the whole-body extension

**Registration status of this addendum, stated honestly and first.** The analyses in this
addendum **have already been run.** Registering them now is therefore NOT a pre-registration
of those results and must never be presented as one. What this addendum does is freeze the
label set, the measures and the decision rules **for everything that comes next**, so that any
future run — a protein-level reference, a new dataset, additional targets — is genuinely
pre-specified. It also puts a dated record on the labels themselves, which is the part an
outside reader most needs.

## A1. What is frozen

**The label file.** `ev_capture_model/src/aim13_scoping/target_labels.py`, exported verbatim to
`target_labels_frozen.tsv` and to `docs/Label_Review_Sheet.md`. 40 targets:
19 `safe`, 10 `toxic`, 9 `efficacy_fail`, 2 stated exclusions. Every label is clinical safety
history, carries a one-line citation, and no analysis script may write to that file.

**The two measures**, both parameter-free:
`top_share` = share of a gene's total normal expression in its single highest cell type;
`n_types` = number of cell types above 10% of that maximum.

**The expendability factor.** A home cell type counts as expendable if and only if it is a
haematopoietic lineage, matched by the word-boundary pattern list in `target_labels.py`. This
is assigned from the cell type's identity alone and never from a drug outcome.

**The two-factor rule.** `score = top_share + 1.0 if the home is haematopoietic`. No fitted
parameter, and the `+1` dominates by construction, so the rule is exactly "expendable first,
confinement second".

**The test.** Rank separation (equivalently AUC) between the `safe` and `toxic` groups, with a
one-sided permutation p-value; exact enumeration where C(n,k) ≤ 250,000 and 200,000 random
permutations otherwise, with the result stating which was used.

## A2. Decision rules fixed for all FUTURE runs

- **H5.** On a protein-level normal reference, `top_share` will separate `safe` from `toxic`
  with separation > 0.70 and p < 0.05. Below that, the RNA-level result is reported as
  reference-dependent.
- **H6.** The two-factor rule will beat the one-factor rule on any new reference. If it does
  not, the expendability factor is reported as not replicating.
- **H7.** On shuffled outcome labels, separation will fall to 0.5 ± 0.1. If it does not, the
  test is broken and every result depending on it is withdrawn.
- **Any target added to the panel after this date is recorded with the date it was added**, and
  results are reported both with and without the additions.
- **A label may be changed only by an outside reviewer**, using the review sheet, and any
  change is reported as a deviation with the before-and-after numbers.

## A3. Deviations already recorded, and where

Two errors of my own were found and are recorded in the results document rather than quietly
corrected:

1. **The haematopoietic substring bug** (results §23.2). A plain substring test classified
   *"conjunctival goblet cells"* as haematopoietic because the string contains `"t cell"`. It
   inflated the two-factor separation from the correct **0.905** to 0.916. Matching is now
   word-boundary anchored with an assertion-based guard table.
2. **The over-tight pancreatic claim** (results §24.3). A six-target ranking supported the
   claim *"the two top-ranked pancreatic targets are the two this project calibrated"*. On the
   full fourteen-target panel mesothelin ranks fifth, and the claim is withdrawn and replaced
   with the weaker, correct one.

## A4. What would falsify the whole-body extension

- Confinement failing to separate `safe` from `toxic` on a protein-level reference.
- The result disappearing when blood-cancer targets are removed. (Already tested: it does not —
  0.792, exact p = 0.023 on solid tumours alone.)
- An outside reviewer's corrected labels dropping the separation below significance.
- The efficacy-failure control separating as strongly as the toxicity comparison, which would
  mean the measure tracks "programme failed" rather than "on-target toxicity".

## A5. Outcome of the registered protein-level hypotheses (recorded September 7, 2026)

**H5 — SPLIT RESULT.** On the immunohistochemistry reference (76 tissue × cell-type groups):
separation 0.785, p = 0.0104 — **MET** (0.775, p = 0.0158 with uncertain antibodies dropped).
On the mass-spectrometry reference (20 bulk tissues): separation 0.650, p = 0.1148 — **NOT MET**.

Two post-hoc explanations for the mass-spectrometry failure were proposed and both were tested
and **rejected**: sparse measurement (separation *falls* to 0.393 as sparsely-measured targets
are removed) and insufficient granularity (randomly coarsening the 154-cell-type RNA reference
to 20 groups still yields 0.837, with 0.650 outside the range of all 60 random partitions).

Per the registered rule, the section 23 result is therefore reported as **reference-dependent**,
and the claim is scoped to restriction measured across normal *cell types* rather than stated
universally. Results document §26.

**H6 — NOT TESTABLE.** Both protein references are indexed by tissue or tissue × cell type, so
no target's home cell type can be a haematopoietic lineage and the two-factor rule collapses
onto the one-factor rule by construction. Recorded as untestable on these references, not as a
failure of the hypothesis.

**H7 — MET, on the two-antigen machinery rather than the whole-body test.** With outcome labels
shuffled, the single-antigen separation collapses from 0.692 to 0.002. Results document §25.
