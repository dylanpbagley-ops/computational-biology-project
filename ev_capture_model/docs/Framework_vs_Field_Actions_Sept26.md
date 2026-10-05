# Where the framework contradicts what the field actually DID

**Written September 26, 2026.** Not what scientists believed — what they *did*: trials run, money
spent, programs launched. Cases where the framework says "no window exists," the field acted
anyway, and reality adjudicated.

**The load-bearing logical point, stated once because everything below rests on it:** *running a
Phase 1 trial is itself documented evidence of prior belief.* Nobody doses human patients with an
agent they expect to poison them. An IND filing, an acquisition, a launched program — each is a
revealed preference, far stronger evidence than a quote in a review. So every toxic target in the
30-target safe/toxic panel is, by construction, a case where the field committed resources and reality
disagreed. What varies is how *loud* each case is, and whether the specific failure mechanism
matches what the framework names.

Five distinct categories below. They are not equally strong and the weak ones are labeled as such.

---

**Label caveat added Sept 29, 2026 (read before using any case below).** The clinical verdicts behind these cases are agent-specific: the same target can be tolerated in one drug format, route or dose and toxic in another. A later data search (`Round7_Data_Search_Findings_Sept29.md`) found that several "toxic" labels are weakly supported or format-dependent in the primary papers (EGFR, EPCAM, NCAM1, FAP, CD47, ERBB2, and CEACAM5's CAR-T vs TCR-T). No label was changed and nothing was re-scored, but for the paper and poster prefer the best-supported cases (CD47 as a *bet*, which the acquisition documents, is fine; HER2 CAR-T fatal event is a documented single-patient case) and describe the rest as "labelled toxic in the frozen panel". The label limitation is disclosed once, in the paper's Discussion.

---

## Category A — Expensive commitments the framework ranks at the bottom

Targets where the field spent heavily, the framework assigns a near-zero ceiling, and the program
failed on on-target toxicity.

### A1. CD47 / magrolimab — the strongest case in the entire project

- **What the field did:** Gilead acquired Forty Seven, Inc. for **$4.9 billion** in March 2020,
  specifically to obtain magrolimab (anti-CD47). Publicly announced, widely covered, heavily
  premium-priced. This is about as unambiguous as revealed preference gets.
- **What the framework says:** `top_share = 0.030` across **120 cell types** — the *lowest
  confinement of any target in the panel*. No expendable home. The framework ranks it dead last.
  The ceiling is not low; it is effectively absent.
- **What happened:** program discontinued 2024. On-target anemia is a documented, mechanistically
  predicted consequence — CD47 is the "don't-eat-me" signal carried on every red blood cell, so an
  agent that blocks it attacks erythrocytes by design, not by accident.
- **Honesty note:** the discontinuation involved futility and an observed increase in death risk,
  not anemia alone. State it as "the program failed and on-target hematologic toxicity was a
  known, mechanistically predicted feature," not "anemia killed it." Overclaiming here would be
  easy and unnecessary — the $4.9B against a last-place ranking carries the point by itself.

### A2. HER2 CAR-T (Morgan et al. 2010)

- **What the field did:** dosed a human patient. HER2 was selected on the field's stated
  criterion — high tumor overexpression, excellent fold-change.
- **What the framework says:** `top_share = 0.041` across **72 cell types**, non-expendable home.
  Near-zero ceiling.
- **What happened:** fatal pulmonary event within hours of infusion. Trastuzumab separately causes
  on-target cardiotoxicity.
- **Why it matters:** this is the cleanest illustration that fold-change and safety are different
  axes. HER2 has an outstanding fold-change and no achievable window.

### A3. CD44v6 / bivatuzumab mertansine
Fatal on-target skin toxicity. Framework: `top_share = 0.061`, 45 cell types, non-expendable home.

### A4. CAIX CAR-T (Lamers et al. 2006)
Grade 2–4 on-target liver enzyme toxicity; trial halted. Framework: 8 cell types but a
non-expendable home — the two-factor rule rejects it on the expendability axis, not confinement.

### A5. EpCAM / catumaxomab
On-target hepatotoxicity; catumaxomab withdrawn. Framework: 35 cell types, non-expendable home.

### A6. NCAM1 (CD56) / lorvotuzumab mertansine
Dose-limiting peripheral neuropathy; program halted. Framework: 31 cell types, neural home.

### A7. EGFR CAR-T / cetuximab
Dose-limiting on-target skin and pulmonary toxicity. Framework: 60 cell types.

### A8. FAP CAR-T
On-target toxicity to FAP+ marrow stromal cells; cachexia and marrow failure. Framework: 11 cell
types, home = fibro-adipogenic progenitors (a marrow mesenchymal population).

**What still needs gathering for A2–A8:** the clinical outcome is already documented and cited in
`target_labels.py`. What is NOT yet assembled is the *prior-optimism* half — the paper, press
release, or review where each target was championed before the trial. This is straightforward
literature work and it is what converts these from "the framework retrodicts failures" into "the
framework contradicts documented commitments." A1 already has this (the acquisition). Doing it for
A2 and A3 would probably suffice; all eight is unnecessary.

---

## Category B — Cases where the framework named the tissue that actually broke

Stronger than predicting *that* something fails: predicting *where*. Computed from RNA alone, with
no clinical input.

**Exact matches (2 of 11):**

| Target | Computed home cell type | Tissue actually damaged |
| --- | --- | --- |
| **CEACAM5** | colonocytes (`top_share` 0.283, 7 types) | severe dose-limiting colitis |
| **F3 / tissue factor** | ocular epithelial cells | boxed-warning ocular toxicity, up to 60% of patients |

**Tissue-level matches (2 more):**

| Target | Computed home | Actual |
| --- | --- | --- |
| **FAP** | fibro-adipogenic progenitors (marrow mesenchymal) | FAP+ marrow stromal cells, marrow failure |
| **NCAM1** | oligodendrocytes, OPCs, Bergmann glia (all neural) | peripheral neuropathy — neural tissue, but CNS glia vs. peripheral nerve is not an exact cell-type hit |

**Honest accounting — the other 7 do NOT match:** ERBB2 (computed esophageal, actual lung/heart),
EGFR (computed syncytiotrophoblasts, actual skin/lung), CA9 (computed early spermatids, actual bile
duct), CD47 (computed alveolar type 1, actual erythrocytes), CD44 (computed pancreatic acinar,
actual keratinocytes), FOLR1 (computed breast lactating cells, actual corneal epithelium), EPCAM
(computed colonocytes, actual hepatocytes/biliary).

**So the specific-tissue prediction lands in roughly 4 of 11. Do not claim more.**

**But there is a real structural reason for the misses, and it should be stated rather than
hidden.** Look at the `top_share` values among the mismatches: 0.041, 0.051, 0.030, 0.061, 0.092.
For a target spread across 30–120 cell types with no cell type holding more than ~9% of total
expression, *"the home cell type" is not a meaningful prediction* — nothing dominates. The
framework's claim for those targets was never "this specific tissue will break." It was "no dose
separates tumor from normal, so something will break." **And something did, in every one of
them.** The home-tissue prediction is only a claim the framework makes for *confined* targets, and
among confined toxic targets it hits CEACAM5 and FAP and misses CA9 and FOLR1 — 2 of 4.

Written that way this is an honest, interesting, defensible result. Written as "the framework
predicts which tissue gets damaged" it is an overclaim a judge will dismantle.

---

## Category C — The inverse: successes that violate the field's own stated criterion

This direction is as strong as Category A and is easy to overlook.

The canonical target-selection review (*J Hematol Oncol* 2019, DOI 10.1186/s13045-019-0758-x)
lists **specificity** as a primary criterion and does not list expendability at all.

- **CD19 is expressed on 100% of normal B cells.** It has essentially *no* tumor specificity. It is
  also the most successful CAR-T target in history.
- The same holds for **CD20, CD22, CD79B, BCMA, SLAMF7, CD38, CD52** — all approved, all sitting on
  a normal lineage the drug wipes out.

By the field's stated criterion these should have been deprioritized. They became the field's
greatest successes. The framework explains why — the home lineage is expendable — and the field's
own criterion does not.

**The structural fact underneath it:** 8 of the 19 safe targets have a haematopoietic home cell
type. **0 of the 11 toxic targets do.** A perfect one-way separator across the labeled clinical
record.

**The sentence for the paper and the interview:** *"The field's greatest successes systematically
violate the field's own stated selection criterion — and the property that actually unites them
isn't tumor specificity, it's that the normal tissue they destroy is one the body can spare."*

---

## Category D — The direct head-to-head, with clinical outcomes as referee

Section 27.2, independent rapid-autopsy cohort, 18 shared targets (9 safe, 9 toxic), same labels, same statistic (numbers updated Sept 29, 2026; the Sept 26 version showed 17 targets, 0.528 / 0.750 / 0.495):

| Selection criterion | Separation | p |
| --- | --- | --- |
| Tumor-vs-normal expression difference (what the field uses, within-patient) | **0.556** — chance | 0.365 |
| Confinement (what the framework uses) | **0.778** (95% CI 0.52–0.96) | 0.025 |

This is the cleanest version of the whole theme: the field's operative criterion, tested against
real clinical outcomes on targets the field actually pursued, performs at chance. Not weakly —
at chance.

Supporting control (section 27.3): confinement predicts on-target toxicity (0.778) but does **not**
predict which targets reached phase 3 (0.448, p = 0.675). So it is tracking biology, not
popularity — checked against a trial registry compiled by an unrelated group.

---

## Category E — Failures possibly misattributed (weakest, flag as speculative)

Section 34's finding: the efficacy-fail group splits into a small immune-lineage subgroup scoring
like the safe targets and a larger group of broadly-expressed solid-tumor antigens scoring like the
toxic ones (ERBB3, AXL among them).

The implication, if it holds: some programs recorded as "failed for lack of efficacy" may actually
have failed because dose had to be capped for tolerability before efficacy could be demonstrated —
meaning the historical record misattributes the failure mode, and the framework predicts the real
one.

**This is a hypothesis, not a result.** It is consistent with the data and it is interesting, but
nothing here establishes it. Present it as a question the framework raises, never as a finding.

---

## Summary ranking for the write-up

1. **CD47 / magrolimab ($4.9B against a last-place ranking)** — lead with this. It is the single
   most vivid thing in the project and it needs no statistical background to land.
2. **Category C, the inverse** — the field's successes break the field's rule; 8/19 vs 0/11.
3. **Category D, the head-to-head** — 0.556 vs 0.778, same targets, clinical outcomes as referee.
4. **CEACAM5 and F3** — named the tissue, from RNA alone, with the honest 4-of-11 accounting.
5. **HER2 CAR-T** — the fatal case, and the clearest fold-change-versus-safety illustration.
6. **Category E** — mention as an open question only.

---

*Assembled September 26, 2026 by cross-referencing `target_labels.py` clinical evidence fields
against computed home cell types in `wholebody_v2.json`, plus verification of the Forty Seven
acquisition. The home-cell-type match rate (4 of 11, 2 exact) was computed for this document and
is not currently recorded anywhere else in the project.*
