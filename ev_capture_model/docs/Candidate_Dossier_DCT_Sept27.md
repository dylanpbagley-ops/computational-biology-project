# Candidate dossier: DCT (dopachrome tautomerase / TYRP2)

**Angle 4, executed 2026-09-27.** Candidate fixed in
`src/aim13_scoping/dossier_precommitment.md` BEFORE any investigation, with the selection
rule and the failure conditions written down first. This is the result of that
investigation, reported as it came out.

**Verdict: the candidate survives, with one substantial risk and a narrower novelty claim
than the screen's output implied.**

---

## Why DCT was selected

From the re-ranked screen (`screen_rerank_weighted.json`, continuous weights gated at
top_share >= 0.25): rank 4 of 214, behind two placental antigens and one immunoglobulin
variable gene, both categories excluded by the pre-registered selection rule.

| property | value |
| --- | --- |
| confinement (`top_share`) | **0.948** -- the most confined candidate in the screen |
| normal cell types above threshold | **1** (melanocytes) |
| cancer elevation | **361x** in skin cancer cell lines |
| rank under binary rule | 135 of 214 |
| rank under continuous weights | **4 of 214** |

The rank change from 135 to 4 is the entire case for Angle 1: under the binary flag DCT was
buried beneath broadly-expressed genes whose only qualification was a spurious haematopoietic
home.

---

## 1. Expendability -- the strongest part of the case

Melanocyte destruction has a named, well-characterised clinical phenotype: **vitiligo**.

Immune-checkpoint-inhibitor-induced vitiligo is extensively documented in melanoma patients,
is managed rather than dose-limiting, and is repeatedly associated with **favourable**
prognosis and treatment response. It is studied as a *biomarker of response*, not as a
toxicity that stops therapy.

This is a stronger expendability anchor than the framework normally gets. For most targets
expendability is an argument from tissue biology. Here there is a human clinical phenotype
showing the exact normal-tissue cost, at scale, already accepted in oncology practice.

## 2. Tumour dependence -- unusually strong, and it addresses antigen escape

Dependence was the one of the four downstream properties that public data can actually
speak to, and DCT does well on it:

- DCT upregulation drives resistance to chemotherapy (cisplatin, carboplatin, methotrexate)
  and to radiation (UV, X-ray).
- Ectopic DCT overexpression **abrogates UVB-induced apoptosis** -- a direct survival
  function, not a passenger marker.
- It remains well expressed in **amelanotic** melanoma cells where tyrosinase is inactive,
  so its role is not merely pigmentation.
- **Its expression stays constant across melanoma progression, unlike other melanoma
  antigens that downregulate.**

That last point matters more than it looks. Antigen escape -- the tumour dropping the target
and regrowing -- is the most common cause of CAR-T relapse, and it is one of the four
properties the framework is blind to. A target the tumour actively depends on, and which is
documented not to downregulate, is structurally resistant to that failure mode.

## 3. Surface accessibility -- THE RISK, and it is real

DCT is a tyrosinase-related protein and its dominant localisation is **intracellular**:
perinuclear, co-localising with the trans-Golgi marker syntaxin 6, with weak overlap with
late melanosome markers.

Plasma-membrane presence IS documented -- DCT is detected in substructures at the plasma
membrane in both murine and human melanoma lines, which is why it passed the screen's
surface filter. But the review reporting this frames anti-DCT antibody approaches as a
*speculative future direction* ("if studies will confirm that DCT is internalized from the
PM, this will open interesting perspectives"), not as established.

**Honest reading: surface presence is real but its abundance is unquantified, and this is
exactly the epitope-accessibility gap the framework cannot see.** If DCT turns out to be
present at the surface only sparsely, a surface-directed agent fails regardless of how
confined the target is. This is the single most likely way the candidate dies, and it is
also precisely what the proposed flow-cytometry experiment would measure.

### Strong supporting precedent from the sibling protein

**TYRP1 -- the same protein family, the same melanosomal biology -- is in active CAR-T
development explicitly targeting its surface expression** (CIRM-funded programme, cutaneous
and rare melanoma subtypes).

This materially de-risks the accessibility concern: it demonstrates that a melanosomal
tyrosinase-related protein can be surface-accessible enough to support a CAR-T programme.
DCT is the same family with substantially higher confinement.

## 4. Novelty -- genuine, but narrower than "untouched"

DCT is **not** an unknown protein. It is a well-established melanoma antigen with a long
immunotherapy history: peptide vaccines, dendritic-cell vaccines loaded with DCT peptide,
gene-modified cellular vaccines, lentiviral and liposomal delivery.

Every one of those routes targets DCT as an **MHC-presented intracellular antigen** -- a
fundamentally different modality from surface-directed binding.

What has not been done is the surface route. The review states directly: *"there are no
reports so far, to our knowledge, about melanoma therapies targeting directly the DCT gene
or protein."* Only a patent claiming antisense targeting of DCT mRNA exists.

**So the defensible claim is: "surface-directed targeting of DCT is unexplored, despite DCT
being the most confined melanocyte-restricted antigen in a genome-wide screen, with
documented tumour dependence and a tolerated normal-tissue phenotype."** It is NOT "nobody
has thought about DCT." Stating the second would be caught immediately and would discredit
the rest.

---

## A limitation this investigation found IN THE SCREEN ITSELF

The screen's "already pursued" filter is `taken = {t[0] for t in TL.TARGETS}` -- it excludes
only the 44 genes in the clinical label panel.

**It therefore has no knowledge of any target in development that is not in that panel.**
TYRP1 CAR-T is in active development; TYRP1 is not in the panel; the screen had no mechanism
to flag it. TYRP1 did not reach the candidate list for some earlier reason (not determined
here), but it could have, and nothing would have caught it.

**Consequence: "214 NEW candidates" overstates novelty.** Some unknown fraction of those 214
are in development and the screen cannot see it. This needs stating plainly in the write-up,
and it is a cheap fix for a future session -- widen the exclusion set using a trial registry
rather than the label panel.

This is a real finding about the screen's limits, surfaced only because the dossier forced a
literature check on one specific candidate. That is an argument for doing dossiers at all.

---

## Summary

| dimension | assessment |
| --- | --- |
| Confinement | **strongest in the screen** (0.948, one cell type) |
| Expendability | **strong, empirically anchored** -- vitiligo is tolerated and marks response |
| Tumour dependence | **strong** -- drives resistance, blocks apoptosis, does not downregulate |
| Surface accessibility | **UNRESOLVED -- the principal risk.** Documented but unquantified |
| Novelty (surface route) | **genuine**, but DCT is a known antigen by other modalities |
| Precedent | **TYRP1 CAR-T in development** -- same family, de-risks accessibility |

**What would settle it:** quantitative flow cytometry for surface DCT on a melanoma line
versus normal cells. That is the wet-lab experiment already scoped in the expansion plan,
and this dossier gives it a specific, pre-committed question rather than a general one.

*Sources consulted: a review of TRP-2/DCT in melanoma (IntechOpen, chapter 56861); the
checkpoint-inhibitor vitiligo literature; a CIRM-funded TYRP1 CAR-T programme description.
The DCT review is not recent, so a future session should re-check whether a surface-directed
DCT programme has appeared since.*
