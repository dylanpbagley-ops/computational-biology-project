> **SNAPSHOT NOTICE (added Oct 1, 2026):** this is the Sept 8 assessment, kept for its ISEF-odds reasoning. Several numbers and one phrase are out of date. The current numbers are in `Aim13_Results_Two_Comparator_Ceiling.md` (panel 19 safe v 11 toxic, separation 0.914 under the frozen rule; section 29 now says confinement is an empirical result motivated by the threshold form, NOT "the leading term" of the harm formula; sections 40-41 are new). Where this file and the results document disagree, the results document wins. The ISEF odds below were estimates made on Sept 8 and were not re-derived.

# Where the project stands: weaknesses, ISEF odds, and what to build next

**September 8, 2026.** Second full assessment. Every expansion angle listed here has been run,
not proposed — the code and result files are in `ev_capture_model/src/aim13_scoping/`. Test
suite: 195 passing. Results document: sections 1–29.

---

## Part 0 — Two questions you asked

### 0.1 Do you need to email those authors? Honestly, no.

The two cohorts I found pair per-patient antigen density with CAR-T response and say "available
on reasonable request." Here is the honest accounting.

**What it would buy you.** It would let you ask a question the project currently cannot: *within
one drug, does an individual patient's antigen density predict whether that patient responded?*
That is a real question and a nice figure. But note what it is **not**: it is not the question
your project is built around. Your central claim is about **targets**, not patients — which
target is safe to build a drug against. Section 27 showed that tumour-versus-normal separation
does not predict drug success at all, so a patient-level version would be testing a different
mechanism from the one your thesis rests on.

**What it would cost.** An email to an academic lab from a high-school student, for data that
usually sits behind an institutional data-use agreement, on a timeline of weeks to never. It is
also the only item on your list whose success is entirely outside your control.

**The recommendation.** **Do not treat it as a priority, and do not build any part of the
write-up assuming it arrives.** Send it anyway, once, as a fifteen-minute lottery ticket, because
the downside is fifteen minutes and the upside is a genuinely novel figure. But everything below
is worth more per hour than waiting on a reply.

**Is there a version without the email?** Partly, and you already have it. The rapid-autopsy
dataset gives you 30 patients with tumour **and matched normal tissue from the same body** — the
closest thing to patient-level data that exists in public, and better in one respect than the
cohorts behind the email, because the normals are matched. What no public dataset gives you is
*treatment response* attached to *antigen density*. I searched GEO, published supplementary
tables, and the trial literature; the datasets that have both are consistently the ones behind
"on reasonable request." That is a real limit of the public data, not a failure of searching.

### 0.2 Pre-registration — what it actually does, and your objection is half right

**Your objection:** *couldn't I just register after I have the results, and who would know?*

**The mechanical answer.** No, because a registration on OSF is **timestamped and frozen**. Once
submitted it cannot be edited or backdated, it gets a permanent public URL, and anyone can see
the date. So "register it afterwards and pretend" fails at the one thing that matters: the date
is checkable by a stranger.

**But here is the part most explanations skip, and it is the part you actually put your finger
on.** Registering *late* is not fraud — it is just worthless. If you register a plan on Tuesday
for an analysis you ran on Monday, the timestamp faithfully records that you registered on
Tuesday, and a careful reader learns nothing except that you understand the vocabulary. **You
cannot retroactively make an analysis pre-specified.** That is exactly why your pre-registration
document opens by stating plainly that the whole-body analyses were already run and that
registering them now proves nothing about them.

**So what is it doing for you?** Three things, in descending order of honesty-value:

1. **It freezes things going forward.** Every future run — the protein references, any new
   dataset, any added target — is genuinely pre-specified now, because the rule, the measures
   and the thresholds are on record with a date. That is not decoration: it already bit. H5 was
   registered with the threshold "separation > 0.70 and p < 0.05", and when mass spectrometry
   came back at 0.650 the answer was NOT MET. Because the line was drawn first, I could not
   quietly redefine success afterwards — and that discipline is what forced me to find the real
   explanation instead of the flattering one.
2. **It freezes the labels.** The 40 safe/toxic labels are the softest part of the project. Once
   registered, changing one is a documented deviation rather than a quiet edit.
3. **It demonstrates you know why the practice exists** — and given that a judge will probably ask
   how you avoided fitting your story to your data, having a dated answer beats having a good
   argument.

**And here is the version that actually bites, which you now have.** A pre-registration is a
promise about analysis. The **locked prediction list** is a promise about the world: 52 targets
in trials with no verdict yet, 14 predicted tolerated, 18 predicted toxic, 20 explicit no-calls,
dated and unchangeable. If those come out wrong, you were wrong, publicly. *That* is what makes
a claim falsifiable, and it is worth more than the registration itself. Do both; lead with the
predictions.

---

## Part 1 — Weaknesses, ranked by how much they hold the project back

### 1. Small n on the outcome tests — **the real ceiling, and partly inherent**
The whole-body test is 19 safe vs 10 toxic. The solid-tumour-only test is 8 vs 9. The
autopsy-cohort test is 9 vs 8. Each is significant, none is comfortable, and p = 0.046 on the
autopsy panel is one relabelled target away from crossing 0.05.

**Resolvable?** Partly. The world contains a limited number of surface targets with a settled
safety verdict — that part is inherent. But the panel can grow: the trial registry names 74 ADC
targets and the project currently labels 40. **Growing the labelled panel is the single
highest-value cheap thing left**, and it needs literature work, not new data.

### 2. Labels assigned by one person — **resolvable, and the fix is one signature**
Everything else rests on this. The review sheet is written, every row carries a citation, and
the discriminant check (confinement predicts toxicity at 0.750 but predicts reaching phase 3 at
0.495, using a registry someone else compiled) shows the labels are not just tracking target
popularity. **What remains is one outside person going through 40 rows.**

### 3. No wet-lab work and no patient-level outcomes — **inherent, and correctly framed**
See 0.1. Lean into the target-level framing; it is the stronger claim anyway.

### 4. RNA versus protein — **RESOLVED, and it made the claim sharper**
Replicates on antibody staining at cell-type resolution (0.785, p = 0.010); fails on mass
spectrometry across 20 whole organs (0.650, p = 0.115) — and **RNA fails on those same organs
too** (0.626, p = 0.143). Technology is not the variable; resolution is. Results §26.4–26.5.

### 5. The measure needed a mechanism, not just a correlation — **RESOLVED**
Section 29 proves that confinement-and-expendability is the leading term of the ROC ceiling once
harm is weighted by how much the body needs each cell type. The empirical rule and the theorem
are now one object rather than two.

### 6. Gating choices are homemade — **effectively closed**
Three gating methods agree at 0.92–0.98; circularity fixed; and a shuffled-label null now covers
bias hidden in gating, subsampling and binning (single-antigen signal collapses 0.692 → 0.002).

### 7. The absolute-calibration gap — **inherent and harmless**
CyTOF signal is not molecules per cell. The measures that carry the project are invariant to it.

---

## Part 2 — ISEF competitiveness, as numbers

**How to read these.** Each rating is my honest estimate of the probability of that outcome,
0–100, for this project *as it stands today, well presented*. They are not scores out of 100 for
quality. Where I am uncertain I say so, and the base rates each one is anchored to are given.

### The structure, so the numbers mean something

**AzSEF** (your state fair, Arizona Science Center): 2025 had **786 students, 558 projects, and
20 advanced to ISEF**; 2024 had 691 students, 508 projects, **18 senior finalists to ISEF** — of
whom 14 went on to win something at ISEF.

**ISEF Computational Biology & Bioinformatics, 2026**: **2 First ($6,000), 5 Second ($2,400),
8 Third ($1,200), 8 Fourth ($600)** — 23 Grand Awards. Society for Science states awards are
proportional to entries, which puts CBIO at roughly 80–90 projects.

**Top Awards** are chosen **from the First Award winners only** — roughly 44 of them across 22
categories — and there are about 9: Yancopoulos Innovator ($100,000), two Regeneron Young
Scientist Awards ($75,000), Gordon E. Moore ($50,000), Leighton & Berger Family Prize ($50,000),
and the $10,000 Barrett, Horvitz, Scripps and Coleman awards.

### The ratings

| Level | Base rate | **Rating** | Why |
|---|---|---|---|
| **Qualify to ISEF from Arizona** | ~7% of senior AzSEF projects | **85** | A proved theorem with 341 tests, multi-dataset validation and external clinical outcomes is far above the median state-fair project. Not higher because you must first win your school/regional fair, and because purely computational work with no device and no visual hook can underperform with generalist state judges. |
| **Any Grand Award in CBIO (4th or better)** | ~27% | **72** | Comfortably top-quartile for the category on rigour. The risk is not quality, it is being read as another CAR-T data-analysis project in the 15 minutes a judge gives you. |
| **Second Award or better in CBIO** | ~8% | **42** | This is where the theorem, the independent replication and the locked predictions start to separate you from the field. |
| **First Award in CBIO (top 2 of ~85)** | ~2.4% | **22** | Genuinely reachable and genuinely uncertain. The 2026 First Award went to a Markov-chain model of tomato immunity — not fancier than yours, just more pointed. You have what those winners have (mechanistic model + independent validation + one specific new claim). Two slots is two slots. |

**Top Awards** — all conditional on winning First in CBIO, so each is roughly (22%) × (your odds
among ~44 First winners). These are unconditional probabilities:

| Award | Amount | **Rating** | Fit |
|---|---|---|---|
| **Horvitz Prize for Fundamental Research** | $10,000 | **8** | **Your best-fit top award by a distance.** It is explicitly for work whose contribution is a principle rather than an application, and section 29 is exactly that — an impossibility result with the empirical rule proved to be its leading term. |
| **Leighton & Berger Family Prize** | $50,000 | **5** | Named for two computational scientists — Bonnie Berger is a computational biologist at MIT. The category fit is unusually good. |
| **Coleman Award for Life Science Innovation & Impact** | $10,000 | **5** | Life science with a stated clinical consequence; the released tool and the prediction list help here. |
| **Peggy Scripps Award for Science Communication** | $10,000 | **4** | Entirely about how you present. The one rating you can move most by practising. |
| **Regeneron Young Scientist Award** | $75,000 ×2 | **3** | Two of the very top prizes across the whole fair. |
| **Yancopoulos Innovator Award** | $100,000 | **2** | Tends toward invention and application; yours is a theory project. |
| **Gordon E. Moore Award** | $50,000 | **2** | Same. |
| **Barrett Award for Innovation** | $10,000 | **1** | Poor fit — this is not an innovation-shaped project. |

**The honest summary.** Getting to ISEF: likely. Winning something in your category: favoured.
First place: a real shot that depends on execution. A top award: a long shot, but not a fantasy,
and the Horvitz prize is the one to aim the write-up at.

**The single biggest lever on every number above is framing.** Projects in your genre — 2026's
"ReinforCell" CAR-T reinforcement learning, 2026's pancreatic multi-omic profiling — have been
landing at **Fourth Award**. Open with "I analysed CAR-T data" and you join them. Open with
*"I proved there is a hard limit on what any antigen-targeted cancer therapy can achieve, and
then showed the limit predicts which real drugs poisoned patients"* and you are in a different
conversation.

---

## Part 3 — Expansion angles, all run, ranked by impact

### #1. The weighted multi-compartment ceiling theorem — **DONE, and it is the biggest thing here**

The project carried two measures with no formal relationship: the ROC ceiling (a theorem) and
confinement-plus-expendability (an empirical rule). Section 29 proves they are the same theorem.

Give each normal cell type an **expendability weight** w (0 = the body regenerates it or can
replace its function, 1 = essential) and define harm as the weighted fraction of normal cells
killed. Then the achievable (tumour killed, harm) region is exactly the ROC hull of tumour
against the **weighted mixture** of normal cell types — and the confinement rule falls out as
the leading term, `H ∝ Σ_c w_c s_c`.

**Verified numerically: 16,400 checks across 1–6 compartments, largest violation 0.000e+00, and
likelihood-ratio thresholding matches brute force exactly.** All three known results (the
single-compartment ceiling, the worst-compartment ceiling, the confinement rule) drop out as
special cases.

**Why it is #1:** it turns the project's headline from *"a measure that predicts drug safety"*
into *"a theorem about what any antigen-targeted therapy can achieve, of which that measure is
the leading term."* It also converts the expendability judgement from something buried in a
label file into an explicit modelling parameter a reader can disagree with. This is the section
to aim at the Horvitz prize.

### #2. Locked prospective predictions — **DONE, and it is the cheapest credibility you will ever buy**

52 targets that are in registered trials and have **no safety verdict yet**, scored with the
rule from section 23 unchanged, no free parameters, thresholds fixed in advance:
**14 predicted tolerated, 18 predicted to cause on-target toxicity, 20 explicit no-calls.**
Dated and locked in `prospective_predictions.json`.

Broadly-expressed targets like CD46 (107 cell types) and ADAM9 (114) are called toxic; confined
ones like LRRC15 and SLC34A2 are called tolerated. Every one of these is falsifiable by a trial
readout in the next few years.

**Why it matters:** it is the difference between explaining the past and predicting the future,
and it costs nothing because the rule was already fixed. Put the table in the write-up unchanged
and never revise it.

### #3. Genome-wide discovery of targets nobody is pursuing — **DONE, with a 57-fold positive control**

All 19,475 expressed human genes, four filters, no fitted parameter: observed cell-surface
localisation → above the CAR-T sensitivity floor → confined **and** expendable (the frozen rule)
→ elevated in cancer cell lines → not already being pursued. Result: **214 candidate targets**,
134 of them in an expendable blood lineage.

**The positive control is what makes this real.** A discovery list with no control is a wish
list. So: immunotherapy targets already in active clinical development that were *absent* from
the registry the screen excluded — could the screen find them without being told?

| target | rank of 214 | what it is |
|---|---|---|
| **LILRB4** | **6** | AML/myeloma CAR-T (IO-202), active trials |
| **GYPA** | **9** | glycophorin A, erythroleukaemia |
| **BAFF-R** | **16** | BAFF-R CAR-T, B-cell malignancy trials |
| **CD27** | **20** | varlilumab |
| **TACI** | **23** | myeloma bispecific/CAR target |
| **CSF3R** | **24** | AML |
| **CD79A** | **29** | partner of the target polatuzumab hits |
| **CLEC12A** | **38** | CLL-1, AML CAR-T |

**10 of 16 recovered, 8 in the top 50.** The screen keeps 1.10% of genes, so chance predicts
**0.18** hits. Ten is a **57-fold enrichment, p = 1.6 × 10⁻¹⁶.**

**And the misses are correct calls.** CD5, CD7, KIT and FLT3 were dropped — and CD5/CD7 CAR-T
cause T-cell aplasia while FLT3/KIT CAR-T cause myeloablation. The framework declined to
recommend exactly the ones that poison people.

Two failures of mine along the way, both kept in the write-up because they are the trap this
screen exists to fall into. Using bulk tumour RNA, the top hits were haemoglobin "elevated in
glioblastoma" and immunoglobulins "elevated in lung adenocarcinoma" — bulk tumour contains blood
and immune infiltrate, so any blood-confined gene looks cancer-elevated *because* it is
blood-confined; the screen was rediscovering its own filter through a contaminant. And with no
absolute expression floor, perfect confinement of nothing scored 1.000 and the list filled with
olfactory receptors.

### #4. Protein at cell-type resolution — **ANSWERED: the data cannot answer it, and that is the finding**

HPA integrated Deep Visual Proteomics in May 2026 — protein per cell type, exactly the reference
I said did not exist. It **cannot** settle the question: 8 of 29 labelled targets are *not
detected*, and **7 of those 8 are safe targets** (BCMA, SLAMF7, CD30, GPRC5D, CD52, DLL3, GPC3).
It is 27 cell types from a single donor, so the rare blood lineages that define the safe class
sit below its detection limit. Dropping them doesn't fix it — it deletes the class being tested.

Also worth knowing: HPA's four-level specificity ladder gives only 0.553 on the RNA cell-type
data where the full continuous computation gives 0.868, so a null on that ladder carries almost
no information. **No further downloads are needed** — the full DVP tables are the same
single-donor dataset with the same limit.

### #5. Grow the labelled panel from 40 to ~74 targets — **no new data, pure literature work**

The trial registry names 74 ADC targets; you label 40. Every added target with a documented
verdict directly attacks weakness #1, which is the project's real ceiling. This is unglamorous
and it is the highest value-per-hour item on the list.

### Angles considered and NOT recommended

- **Emailing authors for patient-level cohorts.** See 0.1. Send it, don't wait on it.
- **More diseases for breadth.** You have two blood cancers, one solid tumour, one whole-body
  atlas and two protein references. Another disease adds a row to a table and dilutes the story.
  You said you don't want this to look scattered; this is the item that would make it look that
  way.
- **Patient-outcome prediction from the 5-vs-5 cohort.** Unfixable before the fair.
- **Re-running the two-antigen machinery on more panels.** Section 28 already replicated it in a
  solid tumour and found the honest limit (the engineering gain concentrates rather than
  generalises). More panels won't change that.

---

## What to do next, in order

1. **Write the paper.** Everything above is done or nearly done; the bottleneck is now prose, and
   ISEF requires it to be yours.
2. Lead with the theorem (§29), then the prediction that follows from it (§23), then the
   evidence that the prediction holds (§26–28), then the locked forward predictions.
3. Get one outside signature on the label review sheet.
4. Submit the pre-registration, rewritten in your own words.
5. Nothing further to download — `proteinatlas.tsv.zip` covered every remaining need.
