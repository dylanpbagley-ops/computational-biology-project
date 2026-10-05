# Work completed — September 9, 2026

Tracking file for what has and has not been done from the Phase 0-2 plan in
`Full_Assessment_Impact_ISEF_Weaknesses_NextSteps_Sept8.md`, plus the 22 weaknesses list in that
same file. Read this before re-doing anything below — check here first so a later session
doesn't repeat work.

**How to read the status marks:** DONE = verified complete and re-checked working. PARTIAL =
real progress made, real gap remains, explained. BLOCKED = could not be done from this session,
reason given. NOT STARTED = correctly sequenced for later, not attempted.

**Update, later on September 9:** a second round of work happened the same day, going deeper on
the three items flagged above as too time-consuming (panel growth, item 5, item 9). See the new
section "Round 2 — deeper work on panel growth, item 5, item 9" near the end of this file for
what changed. The table rows below are left as originally written (that's what was true after
round 1) — do not edit them; read the round-2 section for the current state of those three
items.

---

## Phase 0 — cheap, high-leverage fixes

| # | Item | Status | Detail |
|---|---|---|---|
| 1 | Fix counting-criterion mismatch (W8) | **DONE** | Investigated the actual source of the mismatch: it turned out the figure scripts (fig0, fig5, fig5_poster, fig6, fig8, fig8_poster) had *already* been corrected in a prior session — all currently use 48.1%/53.1pp consistently, and the PNGs postdate their scripts (confirmed via file mtimes, so they're current renders, not stale). The one remaining bug was in `aim12_single_cell_cytof_reordering.py`'s top docstring, which compared the 59.3%/81.8pp full-cohort number against a **stale 9-sample-subset** number (48.1%/59.4pp) as if it were a different *counting convention* on the same cohort — it wasn't; it was a different cohort size under the same convention. Rewrote that docstring block to correctly present both real conventions (A: any-flip + full-cohort-spread, giving 16/27=59.3%/81.8pp; B: strict pairwise, giving 13/27=48.1%/53.1pp on the full 41-patient cohort) with an explicit instruction never to cross-pair a number from one convention with the other's description. Verified: file parses, `test_aim12_single_cell_cytof_reordering.py` (11 tests) still passes, and a repo-wide grep confirms no remaining stale pairing anywhere in active source or docs (only `reproduce_aim12_headline.py`, which legitimately cites the old number as what it's testing against, and historical archive files, which are correctly left untouched per `PROCESS.md`). |
| 2 | Resolve Figure 5 gating-inversion (W9) | **DONE (found already resolved)** | Checked `fig5_aim12_cytof_reordering_example.py` directly: the caption and worked example were already rebuilt on B-gated data using patients whose reordering direction survives gating (Q7998 and Q5971 entering the top 5, displacing patients with genuinely higher gated means, not the old P7676/Q7915 pair that inverted). This is option 3 from the results doc (full rebuild), the strongest of the three options listed — already done, already re-rendered (PNG postdates the script). No further action needed. |
| 3 | Send label review sheet to an outside reviewer | **BLOCKED — needs you** | I cannot send this myself; it needs to go from you to a real person (Professor Wideman was suggested as an obvious candidate given the existing relationship). The sheet itself (`docs/Label_Review_Sheet.md`, 40 rows, one citation and tick-box each) is ready to send as-is — I read it in full this session and it's genuinely usable by a non-computational reviewer with no changes needed. **Action needed from you: send it this week.** |
| 4 | Start growing labeled panel 40→74 | **NOT STARTED** | Confirmed current count is still 40 (`target_labels.py`, `TARGETS` has 40 entries). This is real literature work — finding the clinical outcome (approved/toxic/efficacy-fail) for each of the ~34 additional ADC targets named in the trial registry, with a citation each — and I did not want to rush it in the same session as everything else, since bad labels actively hurt the project (they feed every downstream statistic) more than no labels do. Flagged here as the single most valuable thing to pick up next, in a dedicated session with time to do it carefully. |

## Phase 1 — scientific expansion ("push further" items)

| # | Item | Status | Detail |
|---|---|---|---|
| 5 | Place more real calibrated agents on the ROC | **PARTIAL** | Confirmed the current codebase has exactly two calibrated molecules on file (CLN-978 and blinatumomab, both CD19-directed, giving the three n* thresholds 81/843/8,690). Searched the literature for additional real, quantitative antigen-density activation thresholds for other constructs. Found one genuinely relevant paper: Nature Communications 2025 (doi 10.1038/s41467-025-61427-4), a B7H3 CAR-T functional-avidity study using a QuantiBright flow-cytometry assay, reporting that high-avidity B7H3 CAR-T binders activate against targets carrying **roughly 1,000 molecules/cell**. **This is not a drop-in fit**, and I want to flag that honestly rather than overstate it: unlike CLN-978/blinatumomab (which were fit to a full dose-response curve, giving a precise k_intra), the B7H3 number is an approximate threshold from a different assay methodology, and the paper explicitly states no single universal molecules/cell threshold exists — it's binder-dependent, not construct-independent. **Decision needed from you:** is an approximate, differently-derived threshold like this worth adding to the ROC placement (with the honest caveat that it's not calibrated the same rigorous way), or should this wait until a more directly comparable dose-response calibration is found for a non-CD19 target? I did not add it to any code or figure pending that call. |
| 6 | Check whether any of the 52 locked predictions have resolved | **DONE (checked, none have — expected)** | Pulled the actual predictions from `prospective_predictions.json`: 14 TOLERATED, 18 ON-TARGET TOXICITY, 20 no-call (I want to flag that my Sept 8 assessment doc's summary numbers were correct, but I made an error *this session* while first re-deriving them — I initially filtered on the wrong label string and got 0 toxic before catching and fixing it; noting this so it doesn't look like I'm hiding a mistake). Searched for any public safety readout on any of the 52 genes since the list was dated (Sept 8, one day before this check). Found nothing — expected, since Phase 1 trial safety readouts take months, not a day. **This is not a one-time task — it's an ongoing check.** Recommend re-running this search every few weeks between now and the fair; a single early hit would be a genuinely compelling addition. |
| 7 | Second novelty-search pass (W22) | **DONE** | Searched specifically for the ROC/Neyman-Pearson framing (rather than the CAR-T-modeling framing the first pass used) combined with target-selection and biomarker-panel language. Results were all generic ROC/Neyman-Pearson methodology papers (diagnostic test theory, a SARS-CoV-2 vaccine biomarker paper, general references) — nothing that applies the ROC-as-achievability-bound idea to immunotherapy target selection. This is a genuine negative result in the useful direction: it did not surface any new prior art that threatens the theorem's novelty claim. Recommend treating this as done rather than repeating it verbatim, though a differently-worded pass closer to submission remains cheap insurance. |
| 8 | Re-check protein-at-cell-type landscape (W5) | **NOT STARTED THIS SESSION** | Correctly sequenced for closer to competition per the original plan (the point was to check again later, not now, since DVP itself only appeared in May 2026 and the landscape needs time to move). No action taken, none needed yet. |
| 9 | Healthy bone-marrow reference search (myeloma asymmetry) | **PARTIAL** | Searched for a public healthy human bone-marrow single-cell/CITE-seq reference with surface protein measurements covering plasma cells. Found several real candidates, but none reached the bar of "confirmed usable" in the time available: a 2023 PubMed paper ("Human Bone Marrow Plasma Cell Atlas") looked promising by title but needs direct checking for protein (not just RNA) data and a real accession number; a Nature Immunology 2026 pediatric bone-marrow multiomic atlas surfaced but wasn't checked for adult applicability or protein-level data; and a "novel CITE-seq workflow" bone-marrow atlas (132-plex antibody panel, >80 cell types, described as covering stem/progenitor/immune/stromal cells) turned out on closer inspection to be a **conference webcast announcement, not a paper** — no GEO accession, no confirmed plasma-cell coverage, no data availability statement. **This needs a dedicated follow-up session**, not because it's hard in principle but because properly vetting each candidate (checking it's public, checking it actually has plasma cells and protein-level data, checking sample size) takes real time I didn't have left in this session after items 5-7. |

## Phase 2 — panel growth and label hardening

| # | Item | Status | Detail |
|---|---|---|---|
| 10 | Continue growing 40→74 panel | **NOT STARTED** | Same item as Phase 0 #4 above — listed once, not duplicated work. |
| 11 | Get outside label review signed | **BLOCKED — needs you** | Same as Phase 0 #3 — needs a real person's action, not mine. |
| 12 | Submit OSF pre-registration | **BLOCKED — needs you, and see the walkthrough below** | This is a different document from the "locked predictions" — see the explanation and step-by-step walkthrough at the end of this file. I did not submit anything on your behalf (I have no way to, and this needs to be in your own words per ISEF's AI-authorship rule regardless). |

## Weakness fixes not already covered by Phases 0-2 above

Cross-checking the full 22-item weakness list from the assessment doc against the phase plan:
every weakness that had a concrete, executable-now fix is covered by one of the phase items
above. The remainder are either already resolved (documented as such in the original
assessment — W7, W10, W11, W12, W13, W14, W15, W17, W18, all confirmed still correctly resolved
by this session's fresh test run and verifier pass), inherently blocked for this session (W3,
W4 — need wet-lab or wet-lab-adjacent resources I don't have), or genuinely require your
authorship/discipline going forward rather than a one-time fix (W21, the AI-authorship
boundary — see the flag below; W19, the framing risk — addressed by *how* the paper/poster get
written, not by a script I can run; W16, the disease-specific NOT-gate finding — already
correctly scoped as disease-specific rather than general, nothing further to do; W6 — already
substantially addressed by the existing weighted-ceiling theorem, only the write-up needs to
state the w_c values used, a writing task for Phase 3).

**One additional gap found this session, not on the original 22-item list:** the AI Use
Disclosure Log (`docs/AI_Use_Disclosure_Log.md`) stops at the August 26-27 AND-gate/novelty
session and does not cover anything from the massive September 6-9 Aim13-32 work, or this
session's fixes. This is a real, currently-open gap in an ISEF-required document. I did not
write new entries myself, because this log is specifically your accountability record of how AI
was used, and I did not want to unilaterally decide how to characterize several days of
extensive work without your review of the framing. **Recommend this be a Phase 3 task:
sit down together (or you alone) and write the missing entries before submission**, using the
existing entries' level of detail as a template.

---

## Verification performed this session (baseline confirmed clean)

- `python -m pytest -q` in `ev_capture_model/`: **195 passed**, 2m14s. (Initially failed on
  missing packages — `scipy`/`numpy`/`pandas`/`scikit-learn`/`matplotlib` were not installed in
  this session's environment on your computer; installed them, then all tests passed. This was
  an environment gap, not a code problem — worth knowing in case a future session hits the same
  wall.)
- `python verify_documented_numbers.py` in `src/aim13_scoping/`: **22/22 checks passed**,
  document and code agree.
- Both re-run again after my docstring edit to `aim12_single_cell_cytof_reordering.py`, to
  confirm the edit didn't break anything: still 195/195 and 22/22.

---

## How to submit the locked predictions — walkthrough

You asked specifically about "that locked prediction thingy" — this is `prospective_predictions.json`
(52 targets, dated 2026-09-08, scored TOLERATED/ON-TARGET TOXICITY/no-call with no free
parameters). **Important distinction: this is different from the OSF pre-registration.** They
serve different purposes and both need to happen, but neither has been submitted yet. Here's
what each one is and exactly what to do:

### The locked predictions (`prospective_predictions.json`) — what "submitting" it means

This isn't something you upload to a registry the way a pre-registration is. Its entire value
comes from having a **public, timestamped, unchangeable record that existed before the trials
it predicts resolved** — so "submitting" it means making it exist publicly, dated, before you
or anyone touches it again. The cheapest, most standard way a solo researcher does this:

1. **Put the predictions on OSF as a component of your project, not the pre-registration
   itself.** Go to osf.io, and if you haven't already made a project for this work, make one
   (free, ~5 minutes). Upload `prospective_predictions.json` as a file (or a clean CSV/table
   version of it — up to you, but keep the raw JSON too since it's the actual generated
   artifact). **Make the OSF project or that specific file component public** — OSF stamps
   files with an upload timestamp the moment you make them public, and that timestamp is
   externally verifiable by anyone (including a judge) forever after.
2. **Alternative/addition: register it as a formal OSF pre-registration** using their
   "registries" feature specifically for the predictions (a short, separate registration from
   the main analysis-plan one) — this is stronger than just an uploaded file because OSF
   registrations are cryptographically timestamped and literally cannot be edited after
   registration, not even by you. This is probably worth doing for this specific artifact given
   how much weight the write-up puts on it being genuinely locked.
3. **Either way, the two things a judge needs to be able to verify are: (a) a date, and (b) that
   the content hasn't changed since that date.** OSF gives you both automatically once a file or
   registration is public.
4. **Keep a note of the date and the OSF link** — this is what goes in the paper/poster as the
   citation for "these predictions were locked on [date], see [link]."

### The OSF pre-registration (`OSF_Preregistration_DRAFT.md`) — separate, also not done

This is the *analysis plan* pre-registration — the document that says "here's how we're going to
test whether our ranking predicts clinical outcomes, decided before running the remaining
analyses." It's a genuine OSF Registry submission (not just an uploaded file), and it needs to
be rewritten in your own words before submitting — ISEF explicitly bars AI-written text in
anything submitted, and a pre-registration counts. Steps:

1. Read `docs/OSF_Preregistration_DRAFT.md` in full (I read it this session — it's well-organized:
   the question, four hypotheses, labeling rules, the primary analysis, what would falsify it,
   fixed label lists for both diseases, and an honest "analyses already run before registration"
   section that keeps the record accurate rather than backdating anything).
2. **Rewrite it in your own words.** Not a light edit — actually write it the way you'd explain
   it, using the draft as your reference for content and structure, not as text to copy.
3. Go to osf.io/registries, click "New Registration," and pick a general/open-ended template
   (OSF has several; a generic one works fine for this — you don't need a clinical-trial-specific
   template).
4. Paste in your rewritten version, fill in the standard fields (title, description, your name).
5. **Submit it before running any of the "remaining analyses" the draft refers to** (the healthy-
   PBMC comparator work and beyond) — though per the draft's own honest accounting, most of that
   already ran in the Sept 6-8 sessions. Re-read section 7 of the draft ("Analyses already run
   before registration") carefully before submitting — it already correctly discloses what was
   run before registration versus what the registration will actually govern going forward. If
   more analyses have run since the draft was written (they have — sections 13-32 of the results
   doc happened after this draft), that section needs updating to reflect the true current state
   before you submit, so the registration doesn't misrepresent what it's actually pre-specifying.
   **I'd flag this to you directly: as currently drafted, this pre-registration would be
   submitted quite late relative to how much analysis has already run. That doesn't make it
   worthless — section 7's honesty about that is exactly the right move — but it does mean its
   evidentiary value is smaller than a pre-registration submitted before section 13 would have
   had. Worth deciding deliberately whether to submit it now (with an honest, updated disclosure
   of everything already run) or treat it as primarily insurance for whatever analysis happens
   between now and the fair.**

Both take about 15-30 minutes each once you've written your own version of the text. Neither
requires anything from me beyond what's already prepared.

---

---

## Round 2 — deeper work on panel growth, item 5, item 9

You asked me to go deeper on the three things flagged above as too time-consuming: grow the
panel from 40 toward 74, find a better-calibrated agent for item 5 (you said the B7H3 paper
"looks a little too approximate"), and continue the bone-marrow search for item 9. Here is what
actually happened with each, done honestly rather than padded to look complete.

### Panel growth: 40 → 42, not 74 — and why that's the right number

I went and checked the real 72-gene ADC clinical-trial registry (De Schepper et al. 2025) that
the panel is meant to be drawn from, gene by gene, against the actual clinical literature. The
honest finding: the overwhelming majority of the ~32 unlabeled registry genes **do not have a
settled clinical verdict yet** — they're in ongoing trials with no published approved/toxic/
efficacy-fail outcome. Those already sit correctly in `prospective_predictions.json` as forward
predictions, which is exactly where an undecided target belongs — giving them a fake label just
to hit "74" would have been actively harmful, since every downstream statistic in the results
document is fit against these labels.

I found exactly two genuinely well-documented, verifiable additions:

- **F3 (tissue factor) — toxic.** Tisotumab vedotin (Tivdak) carries an FDA **boxed warning**
  for ocular toxicity (up to ~60% of patients: conjunctivitis, dry eye, corneal
  ulceration/damage, blepharitis), dose-limiting, requiring a dedicated eye-care REMS program.
  Approved despite this toxicity, not because it's absent — the same "approved-but-toxic"
  pattern already used for FOLR1/mirvetuximab in the panel.
- **GPNMB — efficacy_fail.** Glembatumumab vedotin (CDX-011) failed its phase 2b METRIC trial
  (randomized, gpNMB-overexpressing triple-negative breast cancer) on its primary endpoint
  (PFS). Development did not advance. No dose-limiting on-target toxicity signal — a clean
  efficacy failure, not a safety failure.

Both were added through the full pipeline, not just dropped into a spreadsheet:
`target_labels.py` (the frozen label table) → `Label_Review_Sheet.md` (the outside-reviewer
sheet, with full citations and a dated note explaining the panel-growth search) →
`wholebody_v2.py` re-run fresh → the deliberate tripwire test in
`test_wholebody_restriction.py` updated (19/11/10, was 19/10/9) → `verify_documented_numbers.py`
updated → all of Section 23 in `Aim13_Results_Two_Comparator_Ceiling.md` (the 108KB results
document) updated in 8 places with the new numbers, each with a "prior panel, 40 targets"
comparison line so nothing is silently changed. Full pytest suite: 195 passed. Verifier: 22/22
checks pass, document and code agree.

**The genuinely good news:** adding these two targets *strengthened* every headline separation
statistic in Section 23 (one-factor separation 0.868→0.880, two-factor 0.905→0.914, etc.) — not
because I picked favorable targets, but because F3's own independently-computed "home cell type"
(ocular epithelial cells) correctly matches its real toxicity mechanism, and that kind of
independent agreement is exactly what a validity check should show. One honest caveat: the
specificity check in §23.4 crossed from p=0.063 to p=0.043 (borderline-significant) purely from
small-n fragility of adding one target to an 11-target toxic group — I added a paragraph
explaining this is not a claim the underlying separation changed in kind.

**Bottom line: 42 real, defensible, cited targets beats 74 targets where 32 of them would have
been guesses.** This is a case where doing it honestly capped the number lower than you asked
for, and I think that's the right call — happy to discuss if you disagree.

### Item 5: a better calibrated agent for the ROC — no clean winner found, here's the honest tradeoff

You were right that the B7H3 paper was too approximate. I looked hard for something closer to
the CLN-978/blinatumomab bar (a real engineered dose-response titration system, absolute
molecule calibration via a Quantibrite-type bead kit, and an explicit numeric threshold) and did
not find a perfect match. Two real candidates, each falling short in a different, specific way:

**Candidate 1 — Bluhm et al. 2018, BCMA CAR-T** (*Molecular Therapy*, DOI:
[10.1016/j.ymthe.2018.06.012](https://doi.org/10.1016/j.ymthe.2018.06.012), PMID 30078440).
According to PubMed, the peer-reviewed abstract states plainly: "BCMA CAR T cells triggered
target cell lysis with an activation threshold in the range of 100 BCMA molecules." That's a
clean, single, citable number — and BCMA is a **target already in your panel**, which is exactly
the kind of comparability you asked for. The problem: I could not get past the paywall to verify
*how* that number was calibrated. I tried the PubMed full-text API (PMC returns an empty
full-text field for this article — it's Elsevier-licensed, not in the open PMC corpus),
ScienceDirect (403), Cell/Molecular Therapy's own site (403), and the companion ASH abstract on
a related BCMA-density paper (403). Every route to the methods section is blocked. So this
number is real and peer-reviewed, but I cannot personally verify it was derived with
Quantibrite-equivalent absolute-molecule rigor rather than a cruder relative-titration approach
— I only have the abstract's word for "activation threshold."

**Candidate 2 — Majzner et al., "Tuning the Antigen Density Requirement for CAR T Cell
Activity"** (open-access PDF, eScholarship). I got full methods access to this one. It's
genuinely rigorous — CRISPR-KO of native CD19 in NALM6 cells, followed by FACS single-cell
cloning to build a library of clones at defined antigen densities, calibrated in absolute
molecules/cell via the **BD Quantibrite kit** (the same calibration standard your CLN-978/
blinatumomab numbers use). But two things keep it from being a drop-in: (a) it titrated **CD19
and HER2**, neither of which is a currently-uncalibrated target you need — CD19 is already your
best-calibrated target — and (b) critically, **it does not report a single clean threshold
number**. It reports full dose-response curves, and — this is the more important finding —
**the threshold itself is construct-dependent even for the same antigen**: CD19-CD28ζ shows
superior low-density activity compared to CD19-4-1BBζ. That means "the CD19 threshold" isn't
even a well-defined single number in this system; it depends on which CAR construct you pair it
with. That's a genuine scientific finding worth noting for the discussion section, but it makes
this paper unsuitable as a source of one more clean ROC point.

**My honest recommendation:** neither source clears the CLN-978 bar cleanly. Rather than force
an approximate number into the ROC (which is the same mistake you correctly flagged for B7H3), I
recommend **not adding a third calibrated point right now**, and instead treating this as an
open item to revisit — a well-calibrated non-CD19 threshold paper may simply not exist yet in
the literature, which is itself informative (it's part of why the field doesn't already have the
kind of ceiling analysis this project is proposing). If you want, the Bluhm ~100-molecule BCMA
number could be cited in the discussion as *directionally consistent, low-confidence supporting
evidence* (with the paywall caveat stated explicitly) without being plotted as a calibrated ROC
point — that's a real option if you want at least a mention, short of a bad drop-in. **Decision
needed from you** on whether to do that or leave it out entirely. Nothing was added to the
codebase or results document for item 5 this round.

### Item 9: healthy bone-marrow reference — one strong new lead, verified

Went back to the myeloma/plasma-cell asymmetry gap. The most promising thing from round 1 (a
"132-plex antibody panel, >80 cell types" bone-marrow atlas) turned out on closer inspection
*not* to be a mere webcast as I'd flagged — it's a real, peer-reviewed paper I'd simply lost
track of:

**Zhang et al. 2024, "An immunophenotype-coupled transcriptomic atlas of human hematopoietic
progenitors"** (*Nature Immunology*, DOI:
[10.1038/s41590-024-01782-4](https://doi.org/10.1038/s41590-024-01782-4), PMID 38514887,
[PMC11003869](https://pmc.ncbi.nlm.nih.gov/articles/PMC11003869/)). According to PubMed: CITE-seq
performed on primary human bone marrow cells using a machine-learning-optimized panel of **132
antibody-derived tags** (from an initial titration of 266), resolving **>80 stem, progenitor,
immune, stromal and transitional cell states**, explicitly including **plasma cells** as one of
the annotated populations. Donors are healthy adults (the paper's methods describe healthy
bone-marrow aspirate donors, HIV/hepatitis-negative, consented for open data sharing). I found a
public data accession for it: **Zenodo record 15065910**, indexed by the NIAID Data Discovery
Portal as the companion dataset to this exact paper. I could not get a JavaScript-rendered
listing of the individual files at that Zenodo record in this session (the fetch tool can't
execute the page's script), so the exact file-by-file breakdown of what's downloadable (raw
counts vs. processed, whether plasma cells have enough cells to be usable rather than just
present as a label) still needs a five-minute manual check — but the paper itself, the
133-antibody protein panel, the plasma-cell annotation, the healthy-donor status, and the public
accession number are all now independently confirmed, not just a lead.

Second thing I re-confirmed rather than newly found: **Duan et al., "Human Bone Marrow Plasma
Cell Atlas: Maturation and Survival Pathways Unraveled by Single Cell Analyses"** (bioRxiv
preprint, DOI: [10.1101/2023.01.18.524601](https://doi.org/10.1101/2023.01.18.524601), PMID
36711623, [PMC9882341](https://pmc.ncbi.nlm.nih.gov/articles/PMC9882341)). According to PubMed
and a direct read of the preprint text: 5 healthy adult donors (Emory University, 2018-2019),
17,347 bone-marrow antibody-secreting cells profiled by single-cell RNA-seq, with **flow
cytometry validation that directly measured BCMA, CD19, CD138, TACI, OX40, and GITR** — BCMA
being a target already in your panel, which makes this directly relevant to the myeloma
asymmetry question. The honest gap: I checked specifically for a data-availability statement
with a GEO/SRA/EGA accession and found **none** — I read through the preprint's methods and
results and there is no downloadable-data statement in this bioRxiv version, and I could not
find evidence it was ever published in a peer-reviewed journal with an updated data statement
(it still appears everywhere as the 2023 preprint). That means this source is usable as a
**citable published finding** (the reported BCMA/CD19/CD138 flow levels on healthy plasma cells)
but not as a **downloadable raw dataset** to re-analyze yourself.

**Recommendation:** the Zhang et al. Nature Immunology atlas (Zenodo 15065910) is the strongest
candidate found across both rounds — real journal publication, real protein panel, real
accession, plasma cells explicitly annotated. Next concrete step is a direct look at the Zenodo
files (or the paper's supplementary tables) to confirm plasma-cell sample size and whether
protein-level (not just RNA) data is broken out per cell type in the deposited files, rather than
only in the aggregate figures.

---

## Round 3 — item 5 implemented; item 9 lead checked, count too thin, still searching

You confirmed the Bluhm et al. BCMA paper via a working PMC link (mine had been hitting
reCAPTCHA every time; the browser pane got past it). Read the actual methods this time: it
uses **QuantiBRITE PE calibration beads**, the same absolute-molecule standard as CLN-978/
blinatumomab — a real methodological match, not a downgrade like B7H3 was. You said to add it
to the ROC if I judged the data solid enough, and I did.

**What I actually did, and why it's not a simple "add a row":** BCMA isn't measured in the
CD19/DLBCL CyTOF panel that Section 11's table uses, so Bluhm's n\*=100 couldn't go in that
table directly. The correct place is the myeloma CITE-seq cohort (Section 16+), where CD269/
BCMA is a real measured column. I wrote `bcma_operating_point.py`, which applies the identical
test Section 11 already ran for CD19 — checking whether the antigen-density scale in this
dataset can support placing a real molecule-calibrated threshold — to BCMA in the myeloma data
instead. **Honest result: it's the same absolute-calibration gap, found again, independently, in
a second disease and a second technology (CITE-seq vs. mass cytometry).** Raw CD269 ADT counts
across all 11 myeloma samples top out at 357 (one sample, one cell), with per-sample tumor-
compartment medians mostly in the single digits — nowhere near the 100-molecule threshold, and
in some samples the tumor compartment doesn't even clearly exceed T cells on this raw scale,
which is itself a sign the raw counts aren't tracking real molecule numbers here. This is now
written up as new Section 11a in `Aim13_Results_Two_Comparator_Ceiling.md`, right after Section
11, with the same rigor and honesty as the original. All 195 tests still pass.

**Why I'm calling this a genuine positive result even though no operating point got plotted:**
the same gap showing up independently in two diseases, two technologies, two labs is a
stronger, more defensible statement than either instance alone — it says the gap is structural
(how these single-cell datasets are generated), not a fluke of the DLBCL CyTOF panel
specifically. That's a better addition to the paper's honesty section than a forced third ROC
point would have been.

**Item 9:** you said to keep searching rather than build on the 135-plasma-cell Zhang et al.
dataset yet. I pulled the actual GEO accession (GSE245108 — the Zenodo record from round 2 was
only the raw pre-processing flow files, not the finished dataset) and downloaded just the cell-
type annotation file to check numbers before committing to anything larger. Confirmed: 135
plasma cells out of 72,179 total cells across 4 healthy donors — real, labeled, with BCMA/
SLAMF7/CD138/CD38/FcRH5 all in the antibody panel, but thin because this atlas targets CD34+
progenitors, not mature plasma cells, which are rare passengers in it. Per your direction, I
did not download the full matrices or build anything on this yet. Continuing the search for a
larger, plasma-cell-focused healthy bone-marrow reference is the next step — see below.

**Continued the search per your instruction — found a real second candidate.** Checked whether
Duan et al. (the bioRxiv-only "Human Bone Marrow Plasma Cell Atlas" from round 1/2, no
accession found at the time) had a formal peer-reviewed publication yet, since bioRxiv
preprints often gain a proper data statement at journal acceptance. It has: **Cell Reports,
July 2023** (DOI:
[10.1016/j.celrep.2023.112682](https://doi.org/10.1016/j.celrep.2023.112682), PMID 37355988),
and the published version does have a real accession: **GEO GSE230705** (30 samples, 5 healthy
donors × 3 sorted plasma-cell-maturation populations, scRNA-seq + scVDJ-seq, only 156.7MB
total — much smaller and more tractable than Zhang et al.'s multi-GB matrices).

**Honest, direct comparison between the two real candidates now on file:**

| | Zhang et al. 2024 (GSE245108) | Duan et al. 2023 (GSE230705) |
|---|---|---|
| Plasma cells | 135 of 72,179 (0.19%) — incidental, atlas targets CD34+ progenitors | 17,347, all plasma cells — this is a plasma-cell-focused study |
| Healthy donors | 4 | 5 |
| Data type | scRNA-seq **+ CITE-seq surface protein (ADT)**, 132-antibody panel including BCMA/CD138/CD38 | scRNA-seq only in the GEO deposit; BCMA/CD138/CD38/TACI/OX40/GITR were measured by **flow cytometry** (confirmed via the paper's key resources table) but that flow data is not in a downloadable per-cell matrix — only the paper's own summary figures |
| What you'd actually get | A small but real per-cell surface-protein reading (only 135 cells) on plasma cells | A large, real RNA-level reading on 17,347 plasma cells, but protein/molecule-density information only as a citable published summary, not raw per-cell data |

**Neither one is a clean win, and that's worth being upfront about rather than picking one to look decisive:** Zhang et al. has the right data type (protein, not just RNA) but a thin plasma-cell sample; Duan et al. has a large, purpose-built plasma-cell sample but the surface-protein layer isn't in the downloadable data. Kept searching per your instruction and found two more candidates, one of which is now the strongest overall:

**Da Vià et al. 2025 (*Blood*)**, DOI:
[10.1182/blood.2024025643](https://doi.org/10.1182/blood.2024025643), PMID 40009503 — 21 healthy
donors, 234,789 total plasma cells, directly compares malignant vs. genotypically-normal
(polyclonal) plasma cells and explicitly reports BCMA/CD38/CD138 expression differences. Sounds
ideal, but two honest catches: only 7 of the 21 "healthy donors" were newly collected for this
study — the other 14 were pulled from two prior public datasets (one of which is the Duan et al.
data already found). And the accession is **EGA** (European Genome-phenome Archive,
EGAC50000000232), which is normally controlled-access — a data-access-committee application,
not an instant download, unlike the GEO datasets found so far. Also scRNA-seq only, no CITE-seq
protein layer.

**Boiarsky et al. 2022 (*Nature Communications*)**, DOI:
[10.1038/s41467-022-33944-z](https://doi.org/10.1038/s41467-022-33944-z), PMID 36396631 — this
is the strongest candidate found across all three rounds. **GEO GSE193531**, freely and
instantly downloadable (a 52.1MB processed UMI count matrix CSV + a 516.9KB cell-level metadata
CSV — no multi-GB matrices, no controlled access), from the Ghobrial lab at Dana-Farber. **9
healthy/normal bone marrow donors** (labeled NBM-1 through NBM-9 in the data, so healthy cells
are trivially filterable), 29,387 total plasma cells across the full cohort (26 myeloma-spectrum
patients + the 9 healthy donors), with normal and transformed plasma cells compared **from the
same bone marrow biopsy** in the patient samples — an even cleaner within-sample design than
this project's own myeloma cohort. Same honest limitation as Duan et al.: scRNA-seq only, so it
gives RNA-level BCMA/CD38/CD138/etc. expression on real healthy plasma cells, not surface
protein counts — but it is by far the easiest of the four candidates to actually implement, and
the healthy-donor design is purpose-built rather than incidental.

**Where this leaves the search:** four real candidates now on file, no single one has everything
(large healthy-plasma-cell sample size + easy public access + surface protein data, all three at
once). Boiarsky et al. is the best balance of sample size, ease of access, and honest healthy-
donor design; Zhang et al. remains the only one with real per-cell surface protein (ADT) data,
just on very few plasma cells. Continuing to search for something that clears all three bars at
once is the next step if you want it — otherwise Boiarsky et al. is ready to implement now as an
RNA-level healthy-marrow restriction test, structurally parallel to the existing PBMC one.

---

## Round 4 — Da Vià verified (EGA-blocked, confirmed), Boiarsky implemented and written up

You found the Da Vià et al. PDF yourself and put it in the outputs folder, and asked me to check
whether it's the right article and then decide between it and Boiarsky. Then, rather than pick
one of the three options I offered, you asked the more fundamental question directly: is Boiarsky
worth adding at all, and if so, go ahead. Here's what happened with both.

### Da Vià et al. — confirmed correct, confirmed EGA-blocked, not implementable right now

Read the PDF you found in full (DOI 10.1182/blood.2024025643, filename matches — this is the
right paper). Two corrections to my own round-3 summary of it, found on actually reading the
figures rather than relying on the abstract: the real healthy-donor plasma-cell count is
**54,419** (not the 234,789 "total plasma cells" figure I'd cited before, which was the sum
across all groups — HD + MGUS + SMM + MM — not the healthy-donor count alone), and the paper's
own data-availability section confirms, in its own words, that the raw data is deposited at the
**European Genome-phenome Archive, accession EGAC50000000232**, with access "directed to and
fulfilled by" the corresponding author on request — i.e., a real controlled-access gate, not a
formality. This closes the question cleanly: scientifically the strongest of the four candidates
by sample size, but not implementable in this session or probably before the fair, since a data-
access-committee request has its own timeline outside anyone's control here. Filed as closed,
not pursued further.

### Boiarsky et al. — judged worth adding, implemented, written up as new Section 33

Before touching any code, I judged whether Boiarsky was actually worth adding given its one real
limitation (scRNA-seq only, no surface-protein layer) — since you asked me to make that call, not
just implement on request. My reasoning: the entire myeloma side of this project has never had an
**independent, externally-collected** healthy reference — every existing myeloma comparison
(sections 15, 19) is within-patient, tumor vs. normal compartments from the *same* biopsy, which
can't rule out a batch- or cohort-specific artifact the way outside data can. Boiarsky is real,
independent, easy to verify, and directly testable, even though it's RNA-level rather than
protein-level. That's a genuine, if modest, gap worth closing, so I went ahead.

**Downloading it turned out to be the hard part, not the analysis.** This session's tools
couldn't reach any NCBI/GEO host directly (a network block on this side), so getting the two GEO
files (a 516.9KB metadata file and a 52.1MB matrix) onto your computer took several failed
approaches — browser-based streaming, in-browser filtering — before you downloaded both files
yourself and dropped them in the outputs folder. That's what actually unblocked this. Once they
were there, the rest went quickly:

1. **Moved both files** into `ev_capture_data/healthy_bone_marrow_boiarsky/`, verified gzip
   integrity on both (clean).
2. **Inspected the real structure directly** with `zcat`/`awk` instead of guessing: 22,274 gene
   rows (not the full transcriptome — a filtered/HVG set), gene-symbol row names (not Ensembl
   IDs), 29,387 total cells matching the metadata file exactly. All 11 target genes (BCMA/
   TNFRSF17, CD38, SDC1/CD138, SLAMF7, GPRC5D, CD19, CD3D, CD3E, CD14, NCAM1, MS4A1) found on a
   direct grep in under a minute — the earlier in-browser streaming approach had scanned 1,084
   of these 22,274 rows without a hit, which in hindsight was just bad luck with row order
   (genes aren't alphabetically sorted), not a real problem with the data or the gene list.
3. **Wrote `extract_boiarsky_markers.py`**: one-time extraction of the 11 target-gene rows plus
   full cell metadata into a compact 3.4MB JSON, so no script ever needs to re-touch the 52MB
   raw matrix again. Verified matrix column order exactly matches metadata row order (an
   explicit `assert`, not an assumption).
4. **Checked what the dataset actually is before analyzing it — and found an unplanned scope
   limit.** Canonical T-cell/myeloid/B-cell/NK markers are detected in under 12% of all 29,387
   cells (most under 4%), while plasma-cell markers CD38/SDC1 sit at 60-67%. **This file is
   already CD138-sorted/plasma-cell-enriched, not whole bone marrow.** That's consistent with
   the study's own design, but it rules out repeating Section 16's exact "restriction against
   other lineages" AUC test — there aren't enough other-lineage cells left in the file to test
   against. Flagged this honestly rather than forcing the wrong test onto the data.
5. **Ran the test the data can actually support** (`boiarsky_healthy_reference.py`): restricted
   to the 9,329 genuinely healthy NBM donor cells, BCMA (TNFRSF17) is detected in **82.87%** of
   them — more often than the plasma-cell identity markers CD38 (73.78%) and CD138 (57.44%)
   themselves. Checked whether the 15% of cells negative for both CD38 and CD138 were a real
   separate population (they're not — their sequencing depth is roughly half the rest, the
   signature of dropout, and BCMA is still detected in 72.6% of them). For contrast, GPRC5D — a
   real myeloma target known to be lower-abundance — comes back sharply lower (5.67%), which is
   reassurance the pipeline can tell a genuinely sparse marker from a genuinely common one within
   the same sorted population, i.e. the 82.87% BCMA number isn't just a ceiling effect.
6. **Wrote up as new Section 33** in `Aim13_Results_Two_Comparator_Ceiling.md` ("Fifth addendum"),
   matching the document's existing honesty conventions exactly: what the data can and can't show
   stated before the result, the RNA-vs-protein limitation stated up front, and an explicit line
   that this closes the "is BCMA an artifact of this project's own cohort" question, not the "is
   BCMA safe to target" question — the within-patient comparisons and Section 11a's calibration-
   gap finding remain the load-bearing analyses.
7. **Re-ran the full verification suite**: `pytest -q` → 195/195 passed (this addition is
   documentation plus a standalone analysis script, not part of the tracked pipeline, so nothing
   else was expected to move — confirmed nothing did). `verify_documented_numbers.py` → 22/22
   still passing, unchanged (Section 33's numbers aren't wired into that checker, matching how
   Section 11a was also handled — both are narrative findings from standalone scripts, not
   numbers persisted to a JSON file the checker reads back).

**Bottom line: the myeloma-side asymmetry flagged in item 9 is now closed as far as public,
freely-downloadable data allows.** BCMA's behavior in this project's own cohort is now backed by
an independent healthy-donor check, the same way CD19 already was via Section 16. The honest
residual gap — no myeloma-side reference has both a real healthy-donor design *and* surface-
protein (not just RNA) data — is stated plainly in the new section, and Da Vià (the one candidate
that might have closed it) is blocked by EGA controlled access, not by anything fixable from
here.

---

*Companion file to `Full_Assessment_Impact_ISEF_Weaknesses_NextSteps_Sept8.md`. Written
September 9, 2026, after executing everything in Phases 0-2 that could be done without an
outside reviewer, wet-lab access, or your own authorship. Round 2 added the same day after
going deeper on panel growth, item 5, and item 9. Round 4 added the same day after Da Vià
verification and the completed Boiarsky implementation (Round 3 covers item 5's implementation
and the initial Boiarsky candidate search).*
