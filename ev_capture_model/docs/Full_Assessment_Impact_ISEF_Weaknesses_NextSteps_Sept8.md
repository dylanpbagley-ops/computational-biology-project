# Full Project Assessment — September 8, 2026

A synthesis of everything on record (sections 1–32 of the results document, both prior status
memos, and the Aim13-15 proposal), pushed further per your instructions: nothing pre-narrowed,
every weakness surfaced, every viable expansion kept on the table.

**Status update, September 9, 2026: see `Work_Completed_Sept9.md` for exactly what has been done
from the Phase 0-2 plan below and the weakness list, what's blocked and why, and a walkthrough
for submitting the locked predictions. Read that file before re-doing anything listed here.**

---

## Part 1 — How impactful is this project, honestly

Strip away the ISEF framing for a moment and ask what this project *actually found*, as
science.

**The core contribution is a theorem, not a dataset analysis.** Section 29 proves that for any
therapy whose kill probability increases monotonically with antigen density — which covers
essentially every CAR-T, bispecific, and ADC in clinical use — the achievable trade-off between
tumor cells killed and normal cells harmed is bounded by the ROC curve of tumor-vs-normal antigen
distributions, weighted by how expendable each normal cell type is. That means: no amount of
affinity engineering, avidity tuning, or costimulatory redesign can push a bad target past its
ceiling. The ceiling is a property of *where the target lives in the body*, fixed before anyone
builds a construct. This is a genuine, falsifiable, mathematically proven claim, not a
correlation dressed up as one — and it was verified against brute-force search 16,400 times with
zero violations.

**The theorem predicts something checkable, and the check worked.** Sections 22, 23, and 27 show
that the ceiling (tumor-vs-normal separability) does *not* predict which real drugs succeeded
clinically — but restriction to an expendable normal lineage *does* (0.868 separation, p =
0.0003 on 154 cell types; replicated at 0.792, p = 0.023 in solid tumors alone; replicated again
in an independent breast-cancer autopsy cohort at 0.750, p = 0.046). That is three independent
replications of a genuinely counterintuitive finding: **tumor selectivity is not what makes a
cancer drug safe. Hitting a lineage the body can live without is.** This overturns the implicit
assumption behind most target-selection practice in the field (which optimizes for
tumor-vs-normal contrast) and explains, mechanistically, several well-known but previously
unexplained clinical facts: why CD19 CAR-T causes B-cell aplasia despite being a "great" target,
why CD38-directed therapy causes T-cell fratricide, why CEACAM5 causes colitis despite being
confined, why HER2 CAR-T was fatal.

**The genome-wide discovery screen (§30) is real hypothesis-generation with a real positive
control.** Run blind across all 19,475 human genes with no fitted parameters, it recovered 10 of
16 known-but-excluded active clinical targets (57-fold enrichment over chance, p = 1.6×10⁻¹⁶),
and its misses were mechanistically correct rejections (CD5/CD7 cause T-cell aplasia, FLT3/KIT
cause myeloablation — exactly what the framework is built to flag). That is not just a
retrodiction engine; it is a working discovery tool that produced 214 candidate targets, 134 of
them following the confined-and-expendable pattern that defines every approved blood-cancer drug.

**The locked prospective predictions (in Status_v2, referencing `prospective_predictions.json`)
make the whole thing falsifiable going forward**, not just backward-fitting: 52 targets
currently in trials with no safety verdict yet, scored with the frozen rule, no free parameters,
dated. That is the single strongest thing a high-school project can do to earn credibility —
most science fair projects (and a fair number of published papers) never commit to a prediction
that could embarrass them later.

**How this compares to real published work in the space.** The novelty search (§7 of the results
doc) found the closest prior art: SCAN-ACT, Dannenfelser et al. (*Cell Systems* 2020), a 2025
bioRxiv logic-gate paper, and a 2025 *Frontiers in Immunology* review of CAR-T modeling that
directly confirms no reviewed model does what this one does — treats antigen density as a
distribution yielding an escape fraction, derives a selectivity bound from it, or uses
Neyman-Pearson/ROC reasoning for target ranking. That's a real, checked novelty claim, not an
assumed one. If this were submitted as a paper (independent of ISEF), the theorem plus the
three-cohort clinical replication plus the discovery-screen positive control would be a
legitimate, publishable computational/systems-biology contribution — the kind of thing that
lands in a mid-tier journal (Cancer Immunology Research, npj Precision Oncology, or similar) as a
methods/framework paper. That is a genuinely high bar for a project built entirely from public
data by one student, without a lab.

**Where the impact is currently overstated or fragile, in one paragraph (full detail in Part 3):**
the outcome tests that anchor the "predicts clinical success" claim run on n = 5-19 vs 5-10 per
test — real, significant, but small; the safety labels are one person's literature judgment,
unreviewed; nothing here is validated in a wet lab or in a patient; and the protein-vs-RNA
question (whether the whole thing is an RNA artifact) is answered "no" but only partially — it
holds at matched resolution in one protein technology (IHC) and is untestable in another
(DVP, due to detection-limit gaps in exactly the safe class). None of these breaks the core
claim, but a sharp judge will find every one of them within fifteen minutes, so they need
answers ready, not hopes that they go unasked.

**Bottom line on impact:** this is not a "did some data analysis on CAR-T" project. It is a
theorem with a mechanism, three independent clinical replications, a working discovery tool with
a positive control, and a locked set of falsifiable forward predictions. That combination —
proved mechanism + independent validation + one new specific claim — is, per your own prior
session's research into ISEF's own recent winners, exactly the shape of project that wins First
Award in Computational Biology & Bioinformatics. It is a legitimately good project. The work
remaining is overwhelmingly about hardening what exists and writing it up, not about finding
something new to make it good.

---

## Part 2 — Realistic ISEF projection, by level

These are recalibrated from the September 8 status memo (which itself already did real
groundwork — pulled real AzSEF and ISEF CBIO base rates). I'm not discarding that work; I'm
stress-testing it and flagging where I'd move a number and why, then giving a single reconciled
table.

### Sanity-checking the prior estimates

The prior estimates (Qualify 85, Any Grand Award 72, Second+ 42, First 22) are grounded in real,
cited base rates (AzSEF: ~786 students, 20 advance; ISEF CBIO: 2 First / 5 Second / 8 Third / 8
Fourth out of ~85 entries) and a real comparison set (2026 CBIO First Award went to a Markov-chain
model of tomato immunity; recent Fourth Awards include CAR-T RL and multi-omic PDAC profiling
projects). That comparison set is doing real work — it's not vibes, it's "here are the actual
winners in your category in the actual recent years, here's how sophisticated they were."

I'd adjust three things:

1. **The "purely computational, no device, no visual hook" risk deserves more weight, not less.**
   Regional and state fairs (AzSEF) are judged by a broader panel than ISEF's category-specialist
   judges, and a project whose evidence is 32 sections of statistics, theorems, and code is a
   genuinely harder sell to a generalist judge in a 10-minute interaction than a poster with a
   physical prototype or a clear before/after visual. This project's single best asset for that
   audience — the confinement web tool (an interactive, live, gene-lookup tool) — needs to be the
   thing on the table at judging, not a screenshot. I'd keep AzSEF qualification around 82-85 but
   flag this explicitly as the thing to manage, not assume away.

2. **The n-size weakness is real enough to shave a few points off "Second or better," not because
   the science is wrong, but because ISEF judges specifically probe statistical power, and a
   panel of five separate outcome tests each individually modest (n=5v5, n=8v9, n=9v10, n=9v8,
   n=19v10) is a pattern a stats-literate judge will notice and press on. The weighted
   multi-compartment theorem (§29) is the strongest possible answer to that pressure — it reframes
   the small-n empirical tests as *confirmations of a proven theorem* rather than as the whole
   evidentiary base — but only if the write-up leads with the theorem exactly as the memo already
   recommends. I'd keep Second-or-better around 40, contingent on that framing actually landing in
   the final product.

3. **First Award at 22 is defensible but I'd flag the specific risk that moves it up or down by a
   lot: how many of the ~85 CBIO entries this year are AI/genomics-adjacent "big model, big
   claim" projects competing for the same judges' attention.** That's not knowable in advance.
   I'd widen the honest range to 18-27 rather than pin 22, and say plainly that execution and
   framing at the table will move this more than any further analysis will.

### Reconciled table

| Level | Rating (0-100) | Confidence in the rating itself |
|---|---|---|
| Qualify to ISEF from Arizona (via AzSEF) | **84** | High — grounded in real base rates, project is well above median rigor |
| Any Grand Award in CBIO (4th or better) | **74** | High — the theorem + 3-cohort replication clears the bar recent 4th-place projects met |
| Second Award or better in CBIO | **40** | Medium — depends on framing execution and whether n-size objections are handled live |
| First Award in CBIO (top 2 of ~85) | **20** (range 15-27) | Medium-low — real shot, genuinely depends on the field that year and on presentation |
| Any ISEF Top/Special Award (Horvitz-tier) | **8** (Horvitz specifically) | Low-medium — best-fit award exists and is real, but these are long shots by design |

**What would move these numbers the most, in order:** (1) an outside signature on the label
review sheet — removes the single most-repeatable judge objection; (2) rehearsing the "why do
your outcome tests have small n" answer until it's delivered in one confident breath pointing at
§29; (3) having the confinement web tool live and clickable at the table; (4) not opening with
"I analyzed CAR-T data."

---

## Part 3 — Every weakness, exhaustively, with two ratings each

**Rating scale, as you specified.** *Resolvability* (0-10): 0 = impossible to fix before
competition, 10 = trivial, nearly free fix. *Current impact* (0-10): 0 = no effect on the
project's strength right now, 10 = severely undermining the project as it stands. A weakness can
be highly resolvable and low-impact (cheap insurance) or low-resolvability and high-impact (a
real, inherent limitation to frame honestly) — both matter, and I've tried not to let one rating
bleed into the other.

I pulled from all three status documents plus the results document's own self-reported issues
(recorded bugs, wrong diagnoses, honest negatives) and did not filter any of them out, per your
instruction not to prematurely narrow.

### Statistical / evidentiary weaknesses

**W1. Small n on every outcome test.**
The tests validating "confinement/lineage predicts clinical outcome" run at n=5v5 (blood outcome
flip, §18), n=19v10 (whole-body, §23), n=8v9 (solid-tumor-only, §23.6), n=9v8 (autopsy
replication, §27), n=9v10 (efficacy-fail specificity, §23.4, itself not significant). Several of
these are one relabeled target from crossing p=0.05.
*Resolvability: 6* — the trial registry names 74 ADC targets against the current 40 labeled;
growing the panel via literature work is real, bounded, unglamorous effort, not new data
collection. Full elimination is impossible (the world only contains so many targets with a
settled clinical verdict), but meaningful growth is achievable before competition.
*Current impact: 7* — this is the single most likely thing a stats-savvy judge presses on, and
several individual tests are fragile in isolation. Mitigated significantly, but not eliminated,
by the fact that three independent tests all point the same direction and by the theorem (§29)
existing as a non-statistical backstop.

**W2. Safety labels assigned by one person (you/Claude), unreviewed.**
Every outcome test rests on 40 safe/toxic/efficacy-fail labels assigned from literature review,
with no outside check. This is arguably the softest input in the entire pipeline, since it feeds
every downstream statistic.
*Resolvability: 9* — the review sheet already exists (`Label_Review_Sheet.md`) with a citation
and tick-box per row, written specifically for a non-computational reviewer. What's needed is
one qualified adult (a teacher, Professor Wideman, a physician, a knowledgeable relative) to sign
off. This is close to a solved problem mechanically — it just hasn't happened yet.
*Current impact: 6* — not catastrophic (the discriminant check in §27.3 shows confinement
predicts toxicity but not trial-phase-advancement, which argues the labels aren't just tracking
popularity), but "who checked your labels?" is a certain question and "nobody yet" is a weak
answer to walk in with.

**W3. No wet-lab validation, no patient-level outcome data.**
Every result is a reanalysis of existing public data. Nothing here was tested at the bench or
against a real patient's actual drug response.
*Resolvability: 1* — genuinely close to inherent for a high schooler with no lab access and
~6 months. The one lever (emailing authors of patient-level cohorts) is assessed elsewhere in the
project notes as a 15-minute lottery ticket, correctly not worth planning around.
*Current impact: 4* — lower than it sounds, because the project's central claim is explicitly
about *targets*, not patients, and is framed that way already. This is a real, inherent
limitation to state plainly rather than a design flaw — judges respect an honestly-scoped claim
far more than an overreaching one.

**W4. The absolute molecule-count calibration gap.**
CyTOF and CITE-seq signal are not confirmed molecules-per-cell. Section 11 quantifies this
directly: applying real measured activation thresholds (n* = 81, 843, 8,690 molecules/cell) to
the data, only the most sensitive agent lands inside the dynamic range at all — the other two
fall off the top entirely. So *specific drug operating points* are not currently computable, only
*ceilings* (AUCs), which are calibration-invariant by construction.
*Resolvability: 2* — would require wet-lab antibody-binding-capacity bead calibration, which is
out of reach. Not resolvable in the available time.
*Current impact: 3* — low, because the project's own framework anticipates and correctly scopes
around this: the ceiling claim (AUC-based) doesn't need calibration; only the "where does a real
drug sit on the ROC" claim does, and that's already labeled as calibration-dependent and reported
as a negative result rather than glossed over (§11). This is a good example of a limitation that's
already been converted into a strength (an honest negative result) rather than left exposed.

**W5. RNA-vs-protein: partially resolved, but with a real, documented failure mode.**
The headline result replicates on one protein reference (antibody IHC, 0.785, p=0.010) but the
question is explicitly *unanswerable* on the other available protein-at-cell-type reference
(Deep Visual Proteomics) because exactly the safe/confined class (BCMA, SLAMF7, CD30, GPRC5D,
CD52, DLL3, GPC3 — 7 of 8 undetected targets) falls below DVP's single-donor detection limit
(§31). Separately, §32 found that RNA-derived and antibody-derived confinement *values*
correlate weakly (ρ=0.167) even when both predict the same clinical outcome — meaning the two
technologies agree on the conclusion but are not measuring numerically the same underlying
quantity.
*Resolvability: 2* — no existing public dataset covers protein across both blood and solid
tissue at cell-type resolution in the same donors (this was searched for directly and confirmed
absent, §32.1). This is a field-level gap, not a you-level gap.
*Current impact: 4* — moderate. The core outcome claim survives (validated in two technologies),
but the *stronger* claim ("confinement is a property of the protein and RNA is a proxy for it")
cannot be made, and needs to be explicitly walked back in the write-up rather than implied.

**W6. The "does confinement alone explain everything" gap — required the two-factor patch.**
Confinement alone misses CEACAM5 (confined to colonocytes, still caused severe colitis because
colon isn't expendable). The two-factor rule (confined AND expendable) fixes this, but
"expendable" is currently assigned by the researcher from cell-type identity, which is a
judgment call, even if it's a defensible, no-fitted-parameter one.
*Resolvability: 7* — the weighted ceiling theorem (§29) already generalizes this from a binary
judgment call into a continuous, theorem-grounded parameter (w_c) that a reader can set
themselves and disagree with explicitly. That's mostly done; what remains is stating clearly in
the write-up which w_c values were used and why, ideally with a short justification per cell type
rather than an implicit "trust me."
*Current impact: 3* — low, because this has already been substantially addressed by the theorem
work, and the theorem framing converts what would otherwise look like an ad hoc patch into a
principled generalization.

### Methodological / integrity weaknesses

**W7. Aim 12's headline number has never been reproducible from the codebase itself.**
Per §10a, `aim12_single_cell_cytof_reordering.py` contains hardcoded result constants
transcribed from an earlier browser session, with no function in the repo that actually computes
them, and no test checking them. The parser and underlying computation were independently
re-verified this session and matched — but the *original* number's provenance, strictly, was
never checkable from the code alone until now.
*Resolvability: 9* — already substantially addressed (the worked example and full ranking were
independently reproduced and matched this session). What remains is making the reproduction
itself a permanent, tested artifact rather than a one-time check, which is cheap.
*Current impact: 3* — low now that it's been checked and confirmed, but worth explicitly noting
in the write-up as "verified independently this session" rather than silently trusting the old
number.

**W8. A real inconsistency between the documented counting criterion and what the code actually
computed (§10b).**
The write-up's prose describes the 59.3%/81.8pp result one way; the actual code implements a
different (also defensible, but different) convention. Both are internally consistent, but the
number-and-description pairing as currently written is not accurate.
*Resolvability: 10* — this is purely a writing fix. Report 48.1%/53.1pp under the strict
criterion, or keep 59.3%/81.8pp and describe the actual criterion accurately. Trivial to execute,
just needs to actually be done in the poster/paper text and captions (§12 lists exactly which
figures need updating).
*Current impact: 5* — moderate right now only because it hasn't been fixed yet; a judge cross-
checking a caption against a number could catch this, and "we found this and are fixing it" is a
much better story than a judge finding it first. High leverage, cheap fix — do this early.

**W9. Figure 5's hero example inverts under gating.**
The P7676-vs-Q7915 comparison, the project's headline "distribution beats mean" illustration, is
computed on ungated whole-biopsy means and literally reverses direction once gated to pure B
cells (§12). The figure is accurate to what Aim 12 as-run actually computed, but the specific
sentence "even though Q7915 has higher average density" is only true pre-gating.
*Resolvability: 8* — three concrete options are already scoped in the results doc itself,
ranging from a one-clause caption fix to a full figure rebuild on gated data. This is bounded,
known work.
*Current impact: 6* — this is a flagship figure, so a reversal a judge could catch by asking "was
this gated?" is a real risk to project credibility if not addressed before the figure is finalized
for the poster.

**W10. Two-antigen "engineering gain" section was overstated by ~10x in an earlier draft due to a
real bug (only testing "kill where marker is high," missing markers that are informatively
low).**
*Resolvability: 10* — already found and corrected, with the corrected numbers documented
alongside the original mistake rather than silently replacing it (§17). This is functionally
resolved; the only remaining action is making sure the corrected numbers, not the original
inflated ones, are what ends up in the poster/paper.
*Current impact: 2* — low, since it's caught and fixed. The residual risk is purely
administrative (an old draft figure or number surviving into final materials by accident).

**W11. A genuinely wrong self-diagnosis of a "third bug" that turned out not to be a bug (§20).**
Worth listing not because it currently hurts the project, but because it's exactly the kind of
thing a rigorous judge respects seeing documented rather than hidden — and it's a template for how
this project handles being wrong.
*Resolvability: 10* — already resolved and documented.
*Current impact: 1* — essentially none; if anything this strengthens the project's credibility
narrative when included in the write-up's methods/limitations section.

**W12. Gating is homemade percentile-based, not expert manual gating or clustering.**
*Resolvability: 8* — already substantially addressed: three independent methods (percentile,
Otsu, Gaussian mixture) agree at 0.92-0.98 rank correlation on matched files (§14). What remains
is stating this comparison clearly in the write-up as the answer to "how do we know your gating
is right," which is a writing task, not new analysis.
*Current impact: 2* — low; this has effectively been closed out already.

**W13. A substring-matching bug inflated a headline number (0.916 vs the correct 0.905) by
misclassifying "conjunctival goblet cells" as blood lineage because the string contains "t
cell."**
*Resolvability: 10* — already found, fixed, and guarded against regression (word-boundary
matching plus an assertion-based guard table that fails at import if it recurs).
*Current impact: 1* — resolved; only relevant as a documented example of self-correction for the
write-up's integrity narrative.

**W14. Circularity in early gating for the two-antigen work (CD19/CD20 partly defined the B-cell
gate they were then tested within), later fixed by switching comparator or excluding
gate-defining markers (§19).**
*Resolvability: 10* — already resolved.
*Current impact: 1* — resolved; note for completeness only.

**W15. The likelihood-ratio (theoretically optimal) classifier overfits badly and cannot be used
as a practical result — measured overfit gaps of 0.020-0.267 versus ~0.02 for simple
two-threshold gates.**
*Resolvability: 9* — already correctly identified and handled by reporting the practical
(two-threshold gate) number instead of the theoretical optimum, with both explained.
*Current impact: 2* — low; this is already a strength of the write-up (an honest methodological
finding) rather than an unaddressed problem.

**W16. The "exclusion beats requirement" (NOT-gates outperform AND-gates) finding did not
replicate across diseases — 79-86% of pairs in DLBCL vs 51-55% (a coin flip) in myeloma.**
*Resolvability: 3* — this is a genuine biological/statistical finding about generalization limits,
not a bug; "more data" might sharpen it but won't necessarily resolve the disagreement, since it
may be a real difference between diseases.
*Current impact: 3* — low-moderate; already correctly reported as disease-specific rather than
general, which is the honest and defensible framing. The risk is only in accidentally overstating
it as general in the final write-up.

**W17. The pancreatic 6-target panel result ("top 2 ranked targets = the 2 independently
calibrated ones") was falsified when expanded to 14 targets — mesothelin dropped from
effectively-first to 5th.**
*Resolvability: 10* — already caught, corrected, and the weaker/correct claim substituted (§24.3).
*Current impact: 2* — low, resolved; worth keeping visible in the write-up as evidence of
rigor rather than hiding the correction.

**W18. CPTAC pancreatic survival analysis is a clean negative result (mesothelin level does not
predict survival, p=0.35).**
This isn't really a "weakness" in the sense of a flaw — it's an honest null result — but it needs
to be handled carefully in framing since a naive read could look like "the model failed."
*Resolvability: 10* (in the sense that the correct framing is already known and written) — the
result doc already correctly notes this tests a *prognostic* question the framework never claimed
to answer, distinct from the *response* question it does address.
*Current impact: 2* — low, provided the write-up keeps this distinction explicit; if blurred, a
judge could wrongly read it as a failure of the core thesis.

### Framing / presentation weaknesses

**W19. High risk of being mistaken for a generic "I analyzed CAR-T data" project** in a judge's
first 15-20 seconds, which the comparison set shows tends to land at Fourth Award regardless of
underlying rigor.
*Resolvability: 9* — purely a presentation/framing choice, already correctly identified in both
prior memos: lead with the theorem (§29) and the locked predictions, never with the datasets.
*Current impact: 5* — currently unresolved only because the actual poster/paper/pitch haven't
been written yet; this is the single highest-leverage framing fix available and costs nothing but
discipline in the writing.

**W20. Purely computational project with no physical device or visual hook, competing partly
against generalist judges** (especially at the AzSEF/regional level before reaching CBIO
category specialists).
*Resolvability: 6* — cannot manufacture a physical device, but the confinement web tool is a real,
interactive, clickable asset that substitutes reasonably well as a "hook" if it's actually present
and demoed live at the table, not just mentioned.
*Current impact: 4* — moderate; manageable with deliberate staging of the judging conversation
around the tool.

**W21. AI-use disclosure and authorship boundary.**
ISEF bars generative AI from writing the submitted abstract, research plan, or poster text. Every
piece of prose that goes on the poster or in the paper needs to be Dylan's own words, even where
Claude did structure, analysis, and figures. The `AI_Use_Disclosure_Log.md` needs to stay
accurate and current.
*Resolvability: 10* — entirely a matter of discipline in the writing workflow (draft numbers/
structure with Claude, write final prose yourself, keep the log honest as you go).
*Current impact: 3* — low currently since nothing has been improperly submitted, but this is a
disqualification-class risk if mishandled, so it belongs on the list even though it's "just"
a process rule.

**W22. Novelty check is described as "preliminary, not complete" in the earliest proposal
document**, though substantially expanded by results-doc §7 with specific DOIs and specific gaps
identified in each piece of prior art.
*Resolvability: 8* — the hard part (finding and reading the closest prior art) is largely done;
what remains is a final systematic pass closer to competition, ideally with a second search angle
(e.g., searching for "ROC" or "Neyman-Pearson" combined with "immunotherapy target selection" in
a second database) to close out residual risk.
*Current impact: 3* — low-moderate; the existing novelty case is already reasonably strong and
specific, but "we checked thoroughly" needs to be true, not just asserted, by the time this is
defended in front of a judge who may know the literature.

---

## Weakness summary table

| # | Weakness | Resolvability | Current Impact |
|---|---|---|---|
| W1 | Small n on outcome tests | 6 | 7 |
| W2 | Labels unreviewed by outside person | 9 | 6 |
| W3 | No wet-lab / patient-level validation | 1 | 4 |
| W4 | Absolute molecule-calibration gap | 2 | 3 |
| W5 | RNA-vs-protein partially unresolved (DVP detection limit) | 2 | 4 |
| W6 | Expendability weight is a judgment call | 7 | 3 |
| W7 | Aim 12 headline number's reproducibility | 9 | 3 |
| W8 | Documented counting-criterion mismatch (59.3% vs 48.1%) | 10 | 5 |
| W9 | Figure 5 hero example inverts under gating | 8 | 6 |
| W10 | Two-antigen gain overstated (fixed) | 10 | 2 |
| W11 | Wrong self-diagnosed "third bug" (documented) | 10 | 1 |
| W12 | Homemade gating (largely validated) | 8 | 2 |
| W13 | Substring-match bug (fixed) | 10 | 1 |
| W14 | Gating circularity (fixed) | 10 | 1 |
| W15 | Likelihood-ratio overfitting (handled) | 9 | 2 |
| W16 | NOT-gate finding doesn't generalize across diseases | 3 | 3 |
| W17 | Pancreatic 6→14 target correction (fixed) | 10 | 2 |
| W18 | CPTAC survival null result (needs careful framing) | 10 | 2 |
| W19 | Risk of "generic CAR-T project" first impression | 9 | 5 |
| W20 | No physical device / visual hook | 6 | 4 |
| W21 | AI authorship boundary discipline | 10 | 3 |
| W22 | Novelty check needs a final pass | 8 | 3 |

**Reading this table:** nothing here scores a 9-10 on impact — which is itself a meaningful
finding: there is no single fatal flaw. The two highest-impact items (W1, small-n outcome tests;
W2, unreviewed labels) are both real but only moderately resolvable-to-fully-resolvable, and both
already have a known fix in motion. The cluster of "already fixed, low residual impact" items
(W7, W10, W11, W13, W14, W15, W17) is worth noting for a different reason: it's evidence of how
much of this project's rigor comes from catching and correcting its own mistakes — which is a
genuine asset for the write-up's credibility, not just cleanup work.

---

## Part 4 — Next steps, ambitious, sequenced to competition

You asked me not to default to "just write the paper" if there's more expansion worth doing
first. I looked hard for that. Two real candidates for further expansion surfaced beyond what's
already in the prior memos, plus I'm keeping every already-identified high-value expansion. Here
is the full ambitious plan, phased.

### Phase 0 — Cheap, high-leverage fixes (do these first, ~1-2 weeks)

1. **Fix the counting-criterion mismatch (W8).** Pick one convention, update every figure and
   caption that cites 59.3%/81.8pp or the alternative. This is pure writing/figure-regeneration
   work using numbers you already have.
2. **Resolve Figure 5's inversion (W9).** Either add the ungated-data caveat, swap in a
   gating-robust example from the B-gated re-run, or rebuild on gated data. Given this is your
   flagship illustration figure, I'd lean toward the full rebuild if time allows — it removes the
   objection entirely rather than managing it.
3. **Send the label review sheet to an outside reviewer now, not later (W2).** This has the
   longest lead time of any cheap fix because it depends on someone else's schedule. Professor
   Wideman is an obvious candidate given the existing relationship. Send it today or this week;
   let it work in the background while everything else proceeds.
4. **Start growing the labeled panel from 40 toward 74 (W1).** This is literature work you (or
   Claude, checked by you) can do incrementally in parallel with everything else — it doesn't
   block on anything and directly strengthens the project's biggest real weakness.

### Phase 1 — Close remaining scientific gaps (this is the "push further" phase, ~3-5 weeks)

These are the expansions worth doing *before* settling into pure write-up mode, because each one
could materially change or strengthen a claim, and doing them now means the paper is written once
against final results rather than revised after the fact.

5. **Express the ceiling as the achievable operating point of real agents (already flagged in
   §8/§11 as a "what to do next" item, not yet done).** Section 11 shows only one real calibrated
   agent (CLN-978, n*=81) currently lands inside the CyTOF dynamic range; the other two fall off
   the top. This is worth one more attempt: are there *other* published, calibrated CAR-T/
   bispecific activation thresholds in the literature — not just the three currently on hand —
   that would land inside the measurable range for other target/disease pairs in this project
   (myeloma, pancreatic)? Finding even 2-3 more real calibrated agents to place on ROC curves
   would convert "we have a ceiling" into "here's where real drugs sit relative to the ceiling,"
   which is a more concrete, more vivid claim for a judge.
6. **A genuinely new angle: test the theorem's central prediction directly against the locked
   forward-prediction list as it starts to resolve.** The 52 locked predictions were dated this
   week. They won't fully resolve before competition, but it's worth checking now — and then
   again in January — whether *any* of the 52 targets have had a safety readout reported publicly
   since trials started (some Phase 1 CAR-T/ADC trials report safety data relatively fast, within
   months). Even one or two early confirmations before the fair would be a genuinely compelling,
   dated, can't-be-faked addition. This is low-cost (a literature search) and high-upside (a real
   forward validation, not just backward-fitting) — worth checking every few weeks between now
   and the fair rather than only at the end.
7. **A second, independent novelty-search pass (W22), from a different angle than the first.**
   The first pass searched CAR-T modeling reviews and specific keyword combinations. A second
   pass specifically targeting the ROC/Neyman-Pearson framing — search terms like "receiver
   operating characteristic" + "target selection" + "immunotherapy," and separately "achievable
   region" or "Pareto" + "on-target off-tumor" — run across PubMed and a general web/preprint
   search, would close out the residual risk that the core theorem framing (not just the
   CAR-T application) has independent prior art in a different subfield (e.g., diagnostic test
   theory, biomarker panel design) that hasn't been checked yet. This is cheap and directly
   protects the project's single strongest claim (genuine novelty of the theorem).
8. **Revisit the myeloid/solid-tissue protein gap (W5) one more time before declaring it closed.**
   The write-up currently and correctly states no dataset exists that covers protein across both
   blood and solid tissue at cell-type resolution in the same donors. That was checked as of this
   session. It's worth one targeted re-check closer to competition (datasets do appear — DVP
   itself only became available in May 2026) rather than assuming the landscape is frozen for the
   next several months. If nothing new appears, that's fine — it becomes "checked again in
   [month], still nothing," which is itself a stronger statement than checking once.
9. **Consider the myeloma bone-marrow healthy reference gap explicitly flagged in §16.** Section
   16 notes plasma-cell/myeloma targets can't be scored against the healthy PBMC (blood) reference
   because plasma cells live in marrow, not blood — so BCMA, SLAMF7, FcRH5, CD138 all read as
   uninformative there, not as failures. If a public healthy bone-marrow single-cell or CITE-seq
   reference exists (worth a real search — several exist for MDS/AML research, e.g. the Human Cell
   Atlas bone marrow reference, van Galen et al., or Triana et al. *Nature Immunology* 2021), this
   would let the myeloma replication (§15) be extended with the same rigor the DLBCL/lymphoma work
   already has (healthy-reference restriction ranking + forward predictions), which is currently
   asymmetric across the two diseases. This is a real, bounded, worthwhile expansion if a suitable
   dataset exists — genuinely worth a dedicated search session rather than assuming it's the same
   dead end as the DVP protein search.

### Phase 2 — Panel growth and label hardening (runs in parallel with Phase 1, ongoing)

10. **Continue growing the 40→74 target label panel (from item 4), aiming to complete before
    Phase 3 begins.** Every added target with a documented verdict directly strengthens W1, the
    project's single most legitimate weakness.
11. **Get the outside label review signed** (from item 3) — chase this actively if it stalls;
    don't let it become the thing that quietly never happens.
12. **Submit the OSF pre-registration**, rewritten in your own words per the ISEF AI-authorship
    rule, once any Phase 1 expansions that would count as new pre-specifiable analyses are
    decided (so they can be registered before running, preserving their evidentiary strength).

### Phase 3 — Write-up (once Phase 1 findings are in hand, ~6-8 weeks, can start structurally
earlier)

13. **Structure and outline can start now** — the section order is already effectively decided
    (theorem first, then the prediction it makes, then the evidence, then the locked forward
    predictions, then discovery/tool, then honest limitations) — but full prose drafting should
    wait until Phase 1's findings (especially item 5 and item 9, if pursued) are final, so nothing
    gets rewritten from a changed result.
14. **Build/finalize the poster and paper**, written entirely in your own words (Claude's role:
    structure, figure production, feedback on drafts you've already written, not authoring
    submitted text).
15. **Rehearse the judge-facing defense of the known soft spots** — specifically W1 (why small n
    is inherent to some tests but the panel is being actively grown), W5 (why the RNA/protein
    question is honestly scoped rather than dodged), and W19 (the 30-second opening framing). This
    is not a "next step" in the analytical sense but it is a real, schedulable task — a mock judge
    Q&A session, ideally with someone unfamiliar with the project pressing on exactly these points,
    is worth scheduling deliberately in the last 2-3 weeks before competition.

### Phase 4 — Final verification pass (last 1-2 weeks before competition)

16. **Run `python -m pytest -q` (195 tests) and `verify_documented_numbers.py` one final time**
    against the fully assembled poster/paper numbers — not just against the code — to make sure
    nothing that made it onto the final poster drifted from what the code actually produces. This
    is the single cheapest insurance available and should not be skipped even under deadline
    pressure.
17. **Confirm the AI-use disclosure log is complete and accurate** against everything actually
    done across the whole project, not just the most recent sessions.

---

## What I'd explicitly *not* do, and why (kept here rather than silently dropped)

- **Adding more diseases for breadth.** Already correctly identified in Status_v2 as
  counterproductive — the project covers two blood cancers, one solid tumor, one whole-body atlas,
  two protein references; another disease would dilute the story rather than strengthen it. This
  is one place where narrowing scope is the right call, not premature — it's a considered decision
  already backed by the "don't look scattered" concern you raised in the prior session.
- **Waiting on emailed patient-level cohort data as a load-bearing plan.** Already correctly
  scoped as a low-probability, out-of-your-control lottery ticket — send it, don't build the
  schedule around it.
- **Re-running the two-antigen machinery on more panels just for more panels' sake.** Already
  tested in a second disease (myeloma) and a second technology (bulk RNA-seq, breast autopsy
  cohort) — the finding that generalizes (statistical gain near zero) and the finding that doesn't
  (NOT-gate advantage is lymphoma-specific) are both already known. A third replication is
  unlikely to change either conclusion enough to be worth the time versus Phase 1's items.

---

*This document synthesizes `Status_ISEF_Ratings_and_Expansion_v2.md`,
`Project_Status_Weaknesses_and_Expansion.md`, `Aim13_Therapeutic_Index_Ceiling_PROPOSAL.md`, and
`Aim13_Results_Two_Comparator_Ceiling.md` (sections 1-32), all current as of September 8, 2026.*
