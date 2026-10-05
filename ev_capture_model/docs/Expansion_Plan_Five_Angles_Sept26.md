# Expansion plan — five angles, scoped and sequenced

**Written September 26, 2026.** Companion to `Full_Assessment_Impact_ISEF_Weaknesses_NextSteps_Sept8.md`
and `Work_Completed_Sept9.md`. Those two cover weaknesses and what was already executed; this file
covers what is worth ADDING between now and the fairs, and in what order.

**Read this first, before picking anything up:** the standing recommendation has not changed — the
paper and poster are the highest-value use of remaining time, and nothing here outranks them. These
five angles exist because writing the paper well is roughly 6–8 weeks of part-time work, not four
months, so there is genuine room for one or two additions alongside drafting. If anything slips,
the paper wins. An unfinished expansion is worth nothing; an unwritten paper costs everything.

**The filter that produced this list:** an addition qualifies only if it DEEPENS the existing spine
(theorem -> the prediction it implies -> the evidence it holds). Breadth is actively harmful here —
see "What was considered and rejected" at the end. Dylan's own prior sessions flagged the
"don't look scattered" risk and that instinct was correct.

---

## Context a later session needs before starting any of these

**Timeline.** School fair end of January 2027; state fair (AzSEF) March 2027. Roughly four months
from this file's date.

**Assessment as of this file.** Qualification to ISEF ~85%. Any CBIO Grand Award ~70–75%.
Second-or-better ~40%. First Award ~20%. The dominant variable in all four numbers is the write-up
and the judging interview, NOT additional analysis. Say so plainly if Dylan asks again.

**Two findings from the September 26 review session that are not yet in any other document and
should be folded into the paper regardless of which expansions get built:**

1. **A label-sensitivity result that neutralizes the project's second-biggest weakness.** The
   whole-body test was re-run under four adversarial relabelings — the ones a hostile reviewer
   would actually propose. Computed directly from `wholebody_v2.json` + `target_labels.py`,
   20,000-iteration permutation test each:

   | Relabeling | One-factor AUC | Two-factor AUC |
   | --- | --- | --- |
   | As frozen (19 safe v 11 toxic) | 0.880 | 0.914 |
   | NECTIN4 -> toxic | 0.889 | 0.926 |
   | NECTIN4 + CD33 -> toxic | 0.855 | 0.873 |
   | EGFR -> safe (cetuximab is approved) | 0.835 | 0.870 |
   | NECTIN4 + CD33 + TACSTD2 -> toxic | 0.835 | 0.862 |

   Every permutation p stayed below 0.0012. **The headline finding does not depend on Dylan's
   label calls.** This belongs in the paper as a sensitivity table — it converts the softest input
   in the project into a demonstrated robustness result. The script used was ad hoc; a later
   session should rebuild it properly under `src/aim13_scoping/` with a test, so the numbers are
   reproducible rather than one-off.

2. **A structural fact worth stating as a sentence, not an AUC.** Of the 19 safe targets, 8 have a
   haematopoietic home cell type. Of the 11 toxic targets, **zero** do. Not one. That is a perfect
   one-way separator on this panel, and it explains why the +1.0 expendability bonus carries so
   much of the two-factor rule. Phrase for the interview: *"in the entire labeled clinical record,
   no antigen living in an expendable blood lineage has poisoned a patient."*

**Network access, because it will waste a session otherwise.** Both the cloud container and the
device shell are blocked by proxy allowlist from: NCBI/GEO, proteinatlas.org, nature.com,
**crisprsql.com, depmap.org, and ndownloader.figshare.com** (the last three confirmed by direct
`curl` on September 26 — all return connect_rejected). The reliable path for any sizeable file is
to give Dylan a click path and have him download it into the outputs folder himself. This worked
for GSE193531 and is the established precedent. Do not burn a session on blob-download tricks or
base64 round-trips; they have failed repeatedly.

---

## Angle 1 — Continuous expendability weights from an external source

**Highest value of the five. Target window: October–November, alongside paper drafting.**

**The problem it solves.** `w_c` is currently binary: +1.0 if the home cell type is haematopoietic,
0 otherwise, assigned by Dylan from cell-type identity. Two costs follow. First, it is a judgment
call sitting underneath the headline result — weakness W2/W6 in the September 8 assessment.
Second, and more interesting: the section 29 theorem is PROVED for continuous `w_c`, so the project
currently proves something more general than it actually uses. The generality is decorative rather
than load-bearing.

**What to build.** Replace the binary flag with a continuous, externally-sourced weight per normal
cell type, then report the headline result as a sensitivity curve across weighting schemes rather
than a single number.

**Candidate weight sources, in rough order of promise — none of these has been verified yet, and
verifying them is the first task:**
- **Measured tissue turnover / regeneration rates.** Blood and gut epithelium regenerate in days;
  neurons and cardiomyocytes essentially do not. There is a real literature on human cell turnover
  (the Frisen lab's carbon-dating work is the obvious starting point).
- **Clinical replaceability.** Is there an established replacement for losing this lineage —
  immunoglobulin replacement, transfusion, transplant, hormone replacement? This is closer to what
  "expendable" actually means clinically, but it is more judgment-laden than turnover rates.
- **DepMap gene essentiality**, which is also Angle 5's weight source — see the shared-infrastructure
  note below.

**Design notes.**
- Keep the binary version as the primary result and present continuous weights as the robustness
  analysis, at least initially. Do not swap the headline until the continuous version is shown to
  behave at least as well. The binary rule is what every existing number in sections 23–33 uses.
- The deliverable is a **curve or table across weighting schemes**, not a single replacement number.
  The scientific point is "the conclusion is stable across any reasonable way of quantifying
  expendability," which is stronger than any one scheme.
- Whatever source is chosen must be cited and must not have been selected by looking at which
  weighting made the result better. Fix the scheme before computing, and say so in the write-up.

**Why it is worth the most.** It removes the biggest remaining judgment call, makes the theorem's
generality load-bearing, and is itself a contribution — "here is how to put a number on
expendability" is something the field does not currently have.

---

## Angle 2 — Give the ceiling units

**Cheapest real win. Half a day. Do this BEFORE drafting Results, because it changes what gets
written.**

**The problem it solves.** Ceilings are currently reported as AUCs — rank-based, dimensionless,
and abstract to a judge with ten minutes. Nothing in the paper currently tells a reader what the
ceiling MEANS in units a clinician would recognize.

**What to build.** Read a different point off the ROC curves that already exist: for each target,
report the maximum achievable tumor kill at a fixed harm budget. Something of the form *"at most
X% of tumor cells can be killed while harming under 1% of normal cells of type Y."* Pick one or two
harm budgets (1% and 5% are natural) and produce a per-target table.

**Design notes.**
- No new data, no new method, no new assumptions. The curves exist; this is a different readout.
- Keep the harm budget explicit and fixed in advance.
- This likely becomes a headline figure. It is the single most vivid artifact available from work
  already done.
- Do NOT confuse this with absolute molecule-calibrated operating points — that remains blocked by
  the calibration gap documented in sections 11 and 11a. This readout is calibration-invariant
  because it is a fraction-of-population statement, not a molecules-per-cell statement. Say that
  explicitly in the write-up so a judge does not think the calibration gap was quietly ignored.

---

## Angle 3 — Confidence intervals instead of only p-values

**Cheap, purely defensive. One day. Also before drafting Results.**

**The problem it solves.** Small n (19v11, 8v10, 9v8) is the project's real statistical ceiling and
the thing a stats-literate judge will press on. A p-value asserts significance; a confidence
interval shows you know exactly how much uncertainty you have.

**What to build.** Bootstrap or exact confidence intervals on every headline AUC in sections 23,
26, 27, and 33. Report as e.g. "0.914, 95% CI [x, y]" throughout.

**Design notes.**
- Bootstrap over targets (the sampling unit is the target, not the cell type).
- Expect wide intervals. That is fine and it is the honest picture — the point is not to look
  precise, it is to demonstrate calibrated awareness of precision.
- Pair this with the label-sensitivity table above. Together they close the entire small-n line of
  attack: "here is how uncertain the estimate is, and here is how little it moves when you change
  my most questionable inputs."
- Wire the new numbers into `verify_documented_numbers.py` like every other headline figure.

---

## Angle 4 — One deep dossier on a genuinely novel screen candidate

**Medium cost, high upside. December, only if Angle 1 has finished clean.**

**The problem it solves.** The genome-wide screen (section 30) produced 214 candidates and was
validated with known targets. It has never committed to a NEW one. Right now the screen is a
method; a dossier makes it a discovery.

**What to build.** Pick the highest-ranked candidate that is genuinely not in clinical development,
and build a real case for it: independent expression evidence, whether tumors actually depend on it
(DepMap again — see shared infrastructure), surface accessibility, what is already published.
State it as a named, falsifiable prediction alongside the existing locked 52.

**Design notes.**
- **This can fail, and failing must be reported.** A careful look might weaken the candidate. If it
  does, say so — that is a real result about the screen's limits and it is more credible than
  quietly picking a different candidate until one survives. Pick the candidate BEFORE investigating
  it, and record which one was picked and why, before the investigation starts.
- Expression data alone is necessary and nowhere near sufficient — it says nothing about
  internalization, epitope accessibility, shedding, or tumor dependence. The dossier must be framed
  as a hypothesis with a case behind it, not a ready-to-build target.
- Judging value: a judge asking "so what should the field actually target?" and getting a specific
  named answer with evidence is a strong moment. That is the whole point of this angle.

---

## Angle 5 — The gene-editing domain-transfer section

**One Discussion section, ~400 words and one figure. December. Feasibility verified September 26.**

**What this is and — critically — what it is not.** This demonstrates that the theorem's structure
transfers to a second domain. It is **NOT** a second clinical validation, and must never be written
as one. See "the honest limits" below.

**Why it is worth doing.** The mark of a fundamental result is that it escapes its original domain.
Showing that escape once, concretely, is what the Horvitz-tier "contribution is a principle rather
than an application" framing rewards.

**The novelty trade, stated plainly so a later session does not lose it.** This addition carries a
REAL cost, not a hypothetical one: demonstrating the math is domain-general makes the "this is just
the Neyman-Pearson lemma applied" critique more available to a judge, because you cannot show the
theorem travels without also showing its engine is classical. The net is positive — reach gained
exceeds mathematical novelty conceded — **but only if the section stays bounded and the cancer work
stays the headline.** If this grows into a co-equal second half of the project, the trade flips
negative and the project's contribution starts reading as "applied a known theorem in two places"
instead of "found that expendability, not tumor specificity, governs safety." Hold the line at one
section.

**The datasets, verified to exist on September 26:**

- **crisprSQL** (Nucleic Acids Research 2021, database paper PMC7778913, site crisprsql.com):
  144 guide RNAs, 25,632 guide-target pairs, quantitative base-pair-resolved cleavage frequencies
  from GUIDE-seq / CIRCLE-seq / Digenome-seq / SITE-Seq and others. **Already attaches GENCODE gene
  names to off-target sites** — that is the join key to essentiality. Free CSV, no login.
  Caveat: last updated 2020, so five years stale.
- **DepMap Chronos gene-effect scores**: how much cells depend on each gene, measured across ~1,100
  cell lines. Freely downloadable. This is the continuous `w_c` for this domain, sitting in a
  database with no judgment call required. **Shared with Angle 1** — see below.
- **The dataset that actually matters** — Nucleic Acids Research 2023, "Guide-specific loss of
  efficiency and off-target reduction with Cas9 variants"
  (https://academic.oup.com/nar/article/51/18/9880/7249918): 328 sgRNAs across 1,753 target
  sequences, three Cas9 variants (wild-type, HiFi, LZ3), reporting **both on-target efficiency and
  off-target activity per guide per variant**. Supplementary tables S1–S8 plus code at
  github.com/MDhewei/GuideVar.

**Why that third dataset is the right test — this is the key insight of the whole angle.** The
theorem's signature claim is that engineering moves you ALONG a curve set by the target's
distribution and cannot move the curve itself. High-fidelity Cas9 variants are exactly that kind of
engineering. And the paper's own headline finding is that specificity gains come with
**guide-specific** efficiency costs rather than a uniform tradeoff — which is precisely what the
theorem predicts, since each guide sits on its own curve and identical engineering therefore costs
differently depending on which curve you are on. **The project would be explaining a published
empirical observation that its own authors reported without a framework for it.** Explaining someone
else's unexplained result is the strongest available form of external validation.

**The honest limits — both must be written into the section, not glossed:**

1. **The structure does not transplant cleanly.** In the cancer work, on-target and off-target are
   both distributions over cells. In gene editing the on-target is a single locus, so there is no
   distribution on one side and a literal ROC does not exist. What transfers is the theorem's
   substance: the bound on achievable trade-offs between on-target activity and essentiality-
   weighted off-target burden. State the adaptation explicitly. A sharp judge will find this gap if
   it is papered over, and it is entirely defensible if it is named.
2. **No clinical outcome axis exists in this domain.** Approved gene-editing therapies number in the
   single digits and off-target-driven adverse events are essentially unreported. There is no
   equivalent of the 19-versus-11 safety panel. Scope the section as "the ceiling exists and behaves
   as predicted in a second domain," full stop.

**Access:** all three sources are behind the proxy wall (confirmed). Dylan downloads them himself.

---

## Shared infrastructure — read before sequencing

**Angles 1 and 5 both need DepMap essentiality data.** Building the DepMap ingestion once serves
both, and Angle 4's "does the tumor actually depend on this gene" check is a third consumer. This
materially changes the cost calculation: Angle 5 is much cheaper if Angle 1 is built first, and
Angle 4 is nearly free on top of either. Sequence accordingly — do not build three separate
one-off loaders.

---

## Recommended sequencing

| When | What | Why then |
| --- | --- | --- |
| Next 2 weeks | Angles 2 and 3 | Days of work each, and both change what gets written — they must land before Results is drafted, not after |
| October onward | Paper drafting begins regardless | Non-negotiable; everything else runs alongside |
| October–November | Angle 1 | The substantive one; builds the DepMap infrastructure the others reuse |
| December | Angle 5, then Angle 4 if time | Both cheap once Angle 1's infrastructure exists; both are additive rather than load-bearing |

**Standing rule: if anything slips, the paper wins.**

Also still open and unrelated to these five, from the September 10 CURRENT file — do not let the
expansions crowd these out: the label review sheet still needs one outside signature (longest lead
time of anything, has been waiting since September 8), the OSF pre-registration and locked-prediction
timestamp are unsubmitted, the AI-use disclosure log still stops at August 26–27 and does not cover
the September work, and the full 195-test suite plus 22-check verifier have not been re-run
end-to-end since the panel-growth edits.

---

## What was considered and rejected

- **Repositioning the project as a general theory of selective intervention.** This is the failure
  mode Angle 5 is scoped to avoid. The biology — the expendability axis, the clinical validation,
  the counterintuitive result that tumor specificity does not predict tolerability — is what makes
  this a discovery rather than a restatement of classical detection theory. Generality goes in one
  section, never in the headline.
- **Other transfer domains.** Antimicrobial selectivity (pathogen vs. host homolog vs. commensal
  microbiome, with expendability transferring almost literally) and pesticide selectivity (target
  pest vs. pollinators, weighted heavily) are both real fits and are worth ONE SENTENCE each in the
  Discussion as domains satisfying the theorem's conditions. Do not build either. Radiotherapy was
  considered and rejected: that field already does formal multi-organ dose optimization with Pareto
  fronts, so a judge from it will say the work rediscovers standard practice. Diagnostic screening
  was rejected outright — ROC analysis originates there, so no novelty claim is available.
- **Another disease for breadth.** Two blood cancers, one solid tumor, a whole-body atlas and two
  protein references is already enough. A fifth dataset dilutes.
- **The EV-side Aim 3b co-occurrence model.** The project pivoted away; building it now adds a
  second spine.
- **Emailing for patient-level antigen-density + response data.** Already correctly scoped as a
  fifteen-minute lottery ticket. Send it if convenient; build nothing on it.
- **Chasing the label panel toward 74 targets.** Closed September 10 — the remaining registry genes
  are in unconcluded trials, not missed judgments.
- **Reopening the K_INTRA / z / k_geom calibration search.** Permanently closed by user decision.

---

*Written September 26, 2026, following a full review of the project's state (handoff v15, the
September 8 assessment, the September 9 tracking file, the September 10 CURRENT file, the paper
skeleton) plus direct re-execution of the section 23 whole-body test against the real codebase.
Saved in both `ev_capture_model/docs/` and the Claude project, per the convention established
September 9.*

---

## Appendix (added September 26) — the solid-tumor case, and why it makes Angle 1 urgent

Added after a direct query of `target_discovery.json` asking whether the framework predicts any
viable solid-tumor targets. It does — but the current screen ranking cannot surface them. This is
the most concrete argument for Angle 1 found so far, so it is recorded here rather than left in
conversation.

### The framework already has a validated theory of solid tumors

Expendability was never really about blood. Blood is just the most obvious instance of a more
general principle: **the normal tissue can be lost, replaced, or regenerated.** Every approved
solid-tumor target in the 44-target panel reaches its window by one of four narrower routes, and
the panel contains an example of each:

| Route | Approved example | Computed home |
| --- | --- | --- |
| Regenerating epithelium | CLDN18 (claudin-18.2, appr. 2024); TACSTD2; NECTIN4 | foveolar cells; esophageal/urothelial epithelium |
| Surgically dispensable organ | FOLH1 (PSMA) | prostatic glandular cells |
| Developmentally restricted tissue | GPC3 | fetal liver, absent in adults |
| Immune-privileged / limited tissue | DLL3 (tarlatamab, appr. 2024) | oocytes |

This also explains CEACAM5 cleanly: colonocytes *do* regenerate, but barrier function fails during
the gap, so regeneration speed alone is not sufficient — loss must be survivable *while* it
happens. Worth stating in the Discussion; it sharpens the expendability concept beyond
"haematopoietic yes/no."

### The screen produces solid-tumor candidates, but ranks them indefensibly

Of 214 candidates: **105 elevated in at least one solid tumor, 85 solid-only.** Volume is not the
problem. Ranking is.

Top of the solid-only list, all scoring >1.0 **purely on the haematopoietic bonus, not on
confinement**:

```
IL1RL1   score=1.549  conf=0.549  nT=3    home=mast cells
VCAM1    score=1.516  conf=0.516  nT=2    home=kupffer cells
KCNJ15   score=1.368  conf=0.368  nT=3    home=neutrophils
ADAM8    score=1.265  conf=0.265  nT=9    home=neutrophils
SIRPA    score=1.189  conf=0.189  nT=12   home=neutrophils
TGFBR1   score=1.055  conf=0.055  nT=59   home=microglia
AKAP10   score=1.026  conf=0.026  nT=137  home=neutrophils
```

AKAP10 sits at confinement 0.026 across 137 cell types — expressed essentially everywhere — and
outranks genuinely confined targets because its single highest-expressing cell type is a
neutrophil. The bonus is also substantively wrong for several: **microglia are CNS-resident and do
not regenerate**, and neutropenia is dangerous. The binary flag treats neutrophils, microglia and
B cells as equivalent. They are not.

### The real candidates are buried below 1.0 because they get no bonus

| Gene | Confinement | Cell types | Home | Note |
| --- | --- | --- | --- | --- |
| **DCT** | **0.948** | **1** | melanocytes | 361× in skin cancer. Melanocyte-directed immunotherapy is known to cause vitiligo, and vitiligo is tolerated — the framework's own logic predicts a real window here. Most confined solid-tumor candidate in the entire screen. |
| **ERVFRD-1** | 0.901 | 1 | syncytiotrophoblasts | placental; absent in non-pregnant adults (route 3) |
| **PSG1** | 0.845 | 1 | syncytiotrophoblasts | same |
| **GPR50** | 0.860 | 1 | lactotrophs | pituitary |

Correctly rejected by the framework: **ITGB3** (platelets — bleeding), **SLCO1B3** (hepatocytes —
liver failure).

### Implication for Angle 1

Continuous weights would invert this ranking: DCT and the placental antigens rise, AKAP10 and
TGFBR1 fall off. The four solid-tumor routes above become directly expressible as weights
(regeneration rate, surgical dispensability, developmental restriction) rather than being
invisible to a binary haematopoietic flag.

**Two cautions for whoever builds this.**

1. **Re-weighting is not risk-free.** Every headline number in sections 23–33 uses the binary rule.
   Continuous weights could move the 0.914 two-factor separation either direction. Decide *before
   running it* how a decrease gets reported, and report it either way. Keep the binary result as
   the primary until the continuous version is shown to perform at least as well.
2. **This fixes ranking artifacts, not viability prediction.** Re-weighting removes a known
   internal defect. It adds no information about internalization, shedding, epitope accessibility
   or tumor dependence — the screen still cannot predict whether a target yields a working drug,
   and no screen can. Do not let the cleaner list be described as a more *accurate* prediction of
   clinical viability; it is a more *defensible* prioritization.
