# Complete Implementation Plan (Third Session, Full Re-Derivation)

**Purpose of this document:** You asked me to fully read both diagnosis documents, go through everything in them (not just the narrow "next steps" list from the handoff doc), do whatever research was needed, and build a complete plan — then show it to you before implementing anything else. This is that plan. Nothing below has been implemented yet except where explicitly marked "ALREADY DONE" (meaning it was done in the first pass, before your correction, and is included here only so you can see the complete picture in one place).

I built this by making a full inventory of every numbered weak point, every brainstorm category, every restructuring option, and every followthrough-document finding in both source documents — not just the items the handoff document happened to single out.

---

## Part 1: What's already done (from the first pass, for context only)

These were implemented before your correction. Including them here so this document is a complete picture, not because I'm asking you to approve them again.

1. Fixed the "no CRC healthy-vs-cancer data" documentation error.
2. Added Spitzberg et al.'s real single-marker-vs-cocktail capture-loss numbers to Aim 6.
3. Added Guo et al.'s real 4-marker AUC values and cluster11 healthy→CRC decline numbers to Aim 7.
4. Built Aim 9 (the marker-ranking-order check) — two honestly-scoped checks using real data only, including a caught-and-fixed impossible-percentage bug.
5. Added citations for the superselectivity, CAR-T antigen-density, and capture-recapture ecology literatures to `binding_model.py`.
6. Updated the handoff document and ran the full test suite (149 tests passing).

---

## Part 2: New research completed this session (informs everything below)

- **SPR/QCM Langmuir paper** (the single most-hyped lead in both documents): confirmed real via PubMed, but **has no accessible full text anywhere** (not open access, no free version exists). Two related papers checked the same way: also paywalled. One related paper (Sugiyama et al., dual-targeting liposomes) **was** fully readable — real evidence that halving ligand density reduces single-target binding, but it's about a different question (two-ligand cooperative binding) and doesn't hand over a usable density-vs-affinity table.
- **PDAC integrin paper** (PMC11189911): fully read. Real data, but it's Western blot relative-expression data (bulk, not per-particle), so it can't plug into the project's per-vesicle model the way Guo's data does. Usable only as a citable supporting fact.
- **CRC supermeres paper** (PMC12902025, published this January): fully read. A brand-new nanoparticle class study with real head-to-head AUC comparisons across three particle types for the same five markers, in a real 22-CRC-vs-25-healthy cohort. Same limitation as the PDAC paper — different particle type, bulk/SPR measurement, not per-vesicle. Citable as independent confirmation that marker choice changes diagnostic power a lot.
- **Shi et al. 2024 CRC systematic review** (DOI 10.1186/s12885-024-12312-8): confirmed real — a meta-analysis of 48 studies, 7,727 patients. Key number: EV protein panels reach AUC 0.92 vs. 0.87-0.90 for RNA-based markers, and panels beat individual markers. Directly supports the project's "marker choice and combination matter" thesis from an independent, large, aggregated source.

---

## Part 3: The five restructuring options — my evaluation of each

The diagnosis document listed five options (A through E) for changing the project's overall shape, not just adding data. Here's an honest read on each, now that the research is in:

**Option A — reframe the whole project around ONE new prediction, not a list of validations.**
This is a *writing/framing* change, not a code change. Aim 9 already function as this "one new claim" — it just needs to be positioned as the headline of the project's story, with everything else (Aims 1-8) explicitly framed as supporting infrastructure that built toward it. **Recommended: do this** — it's cheap (a rewrite of the summary/narrative sections, not new science) and it's exactly what separated top ISEF finishers from the rest per the research.

**Option B — use the SPR/QCM paper to anchor K_INTRA to a real number instead of a swept range.**
Dead end. The paper has no accessible full text. **Not implementable** — I did not attempt to bypass the paywall, consistent with the restrictions I operate under. This should be marked closed in the handoff doc so a future session doesn't re-attempt it, while noting a truly different cross-domain search (not the same Langmuir paper) could still be tried someday if you want.

**Option C — build a "prediction registry": state 2-3 falsifiable predictions before searching for whether they can be checked.**
Genuinely new and cheap to do. **Recommended: do this.** It's a short, honest document (or a section of an existing one) written BEFORE any further data search, listing specific predictions the model makes that haven't been checked yet, stating plainly which ones currently have no known way to check them. This is a structural strength independent of whether every prediction can be resolved.

**Option D — actively search for a real marker-ranking comparison before assuming none exists.**
This is what this session's research already did. The Shi et al. review and the supermeres paper both count as real (if imperfect) evidence in this direction. **Already substantially addressed** — no further action needed beyond citing what was found (see Part 4 below).

**Option E — small first-hand evidence-gathering (e.g., a publicly time-stamped prospective prediction).**
The diagnosis document itself flags this as your call, not mine to decide. **I'm bringing it to you as a decision point, not a recommendation.** It would mean picking one falsifiable prediction from the registry (Option C) and publicly time-stamping it (e.g., posting it somewhere dated) before a specific upcoming paper's data becomes available to check it against. This is a real, if modest, way to add an original contribution — but it takes deliberate timing and a bit of nerve (a public prediction that could be wrong). I'll ask about this directly below.

---

## Part 4: Concrete, itemized implementation plan

Everything in this section is new work I'm proposing to do, organized by document location so you can see nothing was skipped.

### From the Diagnosis Followthrough document, Part 6.3's priority list:
1. ~~Build the ranking-order check~~ — DONE (Aim 9, first pass).
2. **Cite superselectivity + CAR-T together** — DONE (first pass).
3. **Cite capture-recapture ecology stats** — DONE (first pass).
4. Contact Ma/Wu authors about their five-cancer dataset — the followthrough doc itself says this "should not be treated as a priority." My plan: draft a short, ready-to-send email (I will not send it — sending an email requires your explicit go-ahead per my operating rules) so it's a five-second decision for you later, not a research task.

### From the Novelty Diagnosis document, Part 2 (the "how to strengthen each weak point" section):
5. **Weak Point 1 fix (Option A + C above):** rewrite the project's headline framing to lead with Aim 9 as "the one new claim," and add the prediction registry (Option C) as new, standalone content.
6. **Weak Point 3 documentation fix:** update the handoff doc to say the SPR/QCM K_INTRA search is closed due to a real, confirmed access barrier (paywall), not because the idea failed — an important distinction for a judge or future session to understand.
7. **Weak Point 4 additional citations:** add the Shi et al. systematic review numbers (AUC 0.92 for protein panels) and the supermeres paper's head-to-head AUC comparison as new supporting citations — genuinely new evidence for the "marker choice/combination changes diagnostic power" thesis, found this session, not yet in the codebase.
8. **PDAC integrin paper citation:** add as a supporting citation (not model input) alongside the existing Aim 6/7/8 real-data citations, since it's real evidence of per-marker, per-cell-line ranking differences (from a different measurement method) supporting the same general thesis.

### From Part 2's appendix (the 20-category brainstorm) — items NOT yet acted on:
9. **Categories 12/13 (lateral flow assay math, CTC-chip capture-density literature)** — flagged in the diagnosis doc as "not yet scoped, high priority." I have not scoped these yet. **Proposing to do a first scoping pass** (search + read abstracts, not full implementation) to see if either turns up something as good as the SPR/QCM lead did, since neither has been touched at all yet.
10. Categories left genuinely deprioritized with no new action proposed (already honestly labeled weak/speculative in the diagnosis doc itself, and nothing found this session changes that): clinical trials (1), FDA submissions (3), NCI biobanks (7), population biobanks (9), non-English sources (10, partially folds into Ma/Wu), antibody vendor Kd sheets (17, confirmed structurally dead), investor materials (20).

### Housekeeping (applies regardless of which options above you approve):
11. Update `EV-Capture-CURRENT.md` again to reflect this entire second research pass and whatever you approve below, written the way you asked — specific enough to guide a future session, not so narrow it causes tunnel vision again.
12. Run the full test suite again after any code/citation changes to confirm nothing broke.

---

## Part 5: Decisions I need from you before proceeding

**Question 1 — scope of what to implement.** Everything in items 5-9 above is new work I think is worth doing, but it varies in size (some are quick citation additions, item 9 is a genuine new research pass). Do you want all of it, or should I skip anything?

**Question 2 — Option E (the public prospective-prediction idea).** This is the one item that's structurally different from everything else — it would mean publicly committing to a specific, dated prediction before knowing if it's right. Do you want to explore this at all, or leave it out given how close you are to the fair?

**Question 3 — Ma/Wu author-contact email.** Should I draft it (not send) so it's ready if you want it, or skip it entirely since the followthrough doc itself says it's not a priority?

I'll wait for your answers on these before writing or editing any more code.
