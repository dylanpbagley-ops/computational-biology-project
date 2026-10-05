"""
binding_model.py

Aim 1 of the EV capture project: a model for whether a single extracellular
vesicle (EV) sticks to a ligand-coated capture surface, given how many
receptors that vesicle happens to be carrying.

WHERE THIS EQUATION COMES FROM
-------------------------------
This is not something I invented -- it is the "superselectivity" theory from:

    Curk, Dobnikar & Frenkel, "Design principles for super selectivity using
    multivalent interactions" (arXiv:1608.01222), which is the book-chapter
    version of the theory behind:
        Martinez-Veracoechea & Frenkel, PNAS 2011, 108(27):10963
        Dubacheva, Curk, Frenkel & Richter, JACS 2019, 141:2577

Those papers describe a multivalent PARTICLE with k ligands binding to a
SURFACE with receptor density n_R. Our physical picture is flipped: the
SURFACE carries the ligands (at surface density sigma), and the VESICLE
carries n receptors. The math doesn't care which side is "the particle" and
which is "the surface" -- what matters is which of the two disorder regimes
applies. The paper gives three (Eqs 1.13-1.15 in the text):

    (a) fixed ligand & receptor count, all can reach all       -> combinatorial sum
    (b) receptor count is Poisson-distributed (mobile)          -> (1 + n_R e^-bf)^k - 1
    (c) BOTH ligand and receptor count are Poisson-distributed  -> exp(n_R * k * e^-bf) - 1

A vesicle sitting on a capture surface is case (c): it's a "large colloid"
(their words) with mobile receptors on its own membrane, meeting a surface
that also effectively presents ligands at random, fluctuating local
positions/counts within the contact patch. This matches the handoff's
mechanism (a) -- lateral diffusion lets receptors reshuffle into the contact
zone -- almost exactly, which is reassuring: it means the "mobile" case is
the physically appropriate choice, not just the mathematically convenient
one.

THE EQUATIONS, TRANSLATED TO OUR SYMBOLS
-----------------------------------------
Let:
    n       = number of receptors on THIS vesicle (varies vesicle to vesicle)
    sigma   = ligand surface density on the capture surface (design variable)
    f       = binding free energy of a single receptor-ligand bond (a negative
              number in units of kT; more negative = stronger bond)
    beta    = 1 / (k_B * T) -- inverse temperature. We will absorb beta*f into
              a single number for simplicity (see K_INTRA below).

Paper's Eq 1.15 (case (c), both sides Poisson/mobile):
    q_b = exp( n_tilde * k_tilde * exp(-beta*f) ) - 1

Here n_tilde and k_tilde are MEAN receptor/ligand counts in the contact zone.
In our system, "k" (the mean number of ligands the vesicle's contact patch
sees) is proportional to the ligand surface density sigma times the contact
area. So we fold that proportionality into a single constant and write our
own version as:

    q_b(n, sigma) = exp( n * sigma * K_INTRA ) - 1                    ... (*)

where K_INTRA lumps together: the single-bond binding strength e^(-beta*f),
the contact patch area, and any geometric prefactors. This is a deliberate
simplification: instead of tracking every physical constant separately (bond
free energy, patch area, receptor diffusion rate, etc.), we lump them into
ONE effective parameter, K_INTRA, that has units of (1 / (receptor * ligand
density)). Later, when we calibrate against real data (Aim 3), K_INTRA is
the number we tune. This is standard practice: physicists call it an
"effective" or "lumped" parameter.

STILL UNCALIBRATED AS OF THE v5->v6 HANDOFF TRANSITION -- READ THIS BEFORE
ASSUMING K_INTRA/z/k_geom HAVE BEEN FIT TO REAL DATA
-----------------------------------------------------------------------------
K_INTRA (here), z (capture_probability's activity parameter), and k_geom
(contact_patch_receptors' geometric constant, below) remain UNCALIBRATED
placeholders. A session was spent specifically trying to calibrate them
against Mizenko et al. 2021's ~10x (EpCAM) and ~5x (Her2) fold-difference
numbers, per the v5 handoff's top-priority next step -- but this turned out
to be a logic error in that plan, caught before any fit was attempted: on
re-reading Mizenko's actual methodology, their 10x/5x numbers measure a
CO-OCCURRENCE deficit between two DIFFERENT surface markers (tetraspanin
vs. EpCAM/Her2) across the same EV population, captured by two DIFFERENT
methods (single-tetraspanin vs. genuinely non-specific biotin capture) --
not the capture efficiency of ONE marker's own ligand-receptor pair, which
is what K_INTRA/z/k_geom actually govern. There is no principled way to map
a cross-marker co-occurrence ratio onto single-marker binding-strength
parameters. See src/mizenko_cooccurrence_scoping.py for the full account
and a scoping analysis of what Mizenko's numbers DO constrain (Aim 3b's
co-occurrence/independence model), and the handoff document for the full
story of this correction. A genuine calibration of K_INTRA/z/k_geom still
needs a data source that measures single-marker ABSOLUTE capture yield or
efficiency (not a cross-marker ratio) -- not yet found; see "Immediate next
steps" in the handoff.

UPDATE (v6->v7 handoff transition): a follow-up session searched hard for
exactly that missing data source and did not find a clean one -- three
real candidates were evaluated (Mizenko's own raw dilution data, missing a
ligand-density axis; a 2024 Cornell SPR dissertation with a real Kd but
inaccessible full text this session; a 2026 solution-phase FIDA paper with
real numbers but the wrong physical regime). Rather than force an
imperfect fit, the user chose the documented physical-plausibility-bounds
fallback. See src/k_intra_z_plausibility_bounds.py for the full account,
the real numbers found (including a genuine EV-specific Kd=1.92nM anchor
and a convergent cross-platform concentration-regime finding), and exactly
why each candidate source didn't resolve to a clean fit. K_INTRA, z, and
k_geom are STILL uncalibrated placeholders after this -- the plausibility
file bounds the RANGE they should plausibly sit in, it does not set them.

UPDATE (v7->v8 handoff transition): the Lopez Baltazar dissertation named
above as the single most promising unresolved lead WAS retrieved this
session (via Claude-in-Chrome, connected this time). This closed the
SPECIFIC gap that blocked it in v7: the SPR chip's actual antibody surface
density is now known (1.11x10^12 Abs/cm^2 max close-packed; 4.44x10^11
Abs/cm^2 at the "40% coverage" condition that gave Kd=1.92nM), not just
"an unspecified maximum." HOWEVER, this still does NOT resolve to a clean
K_INTRA fit -- checked carefully, not assumed: the retrieved Kd is a 3D
bulk EV concentration (ensemble-averaged over an uncontrolled RAW 264.7 EV
population, folding in mass transport), not a 2D areal binding constant,
so it calibrates z's concentration axis (which it does, and does well --
see the plausibility file's strengthened concentration-axis convergence)
rather than K_INTRA directly. What DID improve: K_INTRA's plausibility
bound can now be expressed in real physical units (cm^2/antibody) using a
real sigma, rather than left as a bare dimensionless product -- see
src/k_intra_z_plausibility_bounds.py's k_intra_bounds_from_real_sigma().
K_INTRA, z, and k_geom remain formally uncalibrated placeholders here in
binding_model.py after this session too -- this update tightens the units
and confidence of the bound, it does not set a number.

Then the Langmuir-isotherm wrapper (paper's Eq 1.8) converts this
bound-state partition function into an actual bound PROBABILITY:

    theta(n, sigma) = z * q_b / (1 + z * q_b)                         ... (**)

z is a dimensionless "activity" that depends on how much time/opportunity
the vesicle has had to find the surface and start forming bonds in the
first place (concentration of vesicles, incubation time, diffusion to the
surface, etc. -- everything OTHER than the multivalent bond-counting itself).
For a single vesicle already in contact with the surface for a given
incubation time, z is treated as a constant scaling factor.

WHY THIS MATTERS PHYSICALLY
----------------------------
Look at what n does in equation (*): it sits INSIDE an exponential. Doubling
the receptor count does not just double the chance of binding -- it can
multiply it by a large factor, because more receptors means combinatorially
more ways to find a bond. This is exactly the nonlinearity the whole project
is about: fold-change in mean receptor count is NOT the same as fold-change
in capture probability, because the capture probability depends
exponentially, not linearly, on n.

WHAT THIS FILE CONTAINS
-------------------------
- q_bound(n, sigma, k_intra): the bound-state partition function, Eq (*)
- capture_probability(n, sigma, k_intra, z): the full theta(n, sigma), Eq (**)
- langmuir_probability(sigma, k_intra, z): the n=1 special case (see below)

WHY WE ALSO NEED THE n=1 (LANGMUIR) CASE
-------------------------------------------
The handoff requires a unit test: valency = 1 (n=1 receptor) must reduce to
"simple Langmuir binding." Looking at Eq (*): when n=1,
    q_b(1, sigma) = exp(sigma * K_INTRA) - 1
and the textbook Langmuir formula is linear in sigma instead:
    q_b_langmuir(sigma) = K_INTRA * sigma
These are NOT the same function -- exp(x) - 1 only equals x approximately,
for small x = sigma * K_INTRA (this is the first term of the Taylor series
exp(x) - 1 = x + x^2/2 + x^3/6 + ...). So our n=1 formula only recovers
EXACT Langmuir binding in the dilute / weak-binding limit (small sigma *
K_INTRA), not at every ligand density.

THIS IS A KNOWN, DOCUMENTED APPROXIMATION, NOT A BUG
-------------------------------------------------------
While testing this, the n=1 case was checked against an independently
written langmuir_probability() function (see below) across a range of
densities, and the two genuinely disagree by a few percent at moderate
sigma * K_INTRA (about 5% relative error around sigma*K_INTRA = 0.1), with
the gap shrinking proportionally to sigma * K_INTRA as density decreases
(confirmed numerically: relative error scales linearly with sigma*K_INTRA,
exactly as expected from the Taylor expansion above). This traces back to
q_bound()'s functional form itself: exp(n * sigma * K_INTRA) - 1 is the
"both receptors AND ligands are mobile/Poisson-distributed" regime from
Curk, Dobnikar & Frenkel (their Eq 1.15) -- which is the physically
appropriate choice for a diffusing EV membrane meeting a surface (see
module docstring above), but is itself already an approximation to the
exact combinatorial theory, valid in the same "many weak, reshuffling
bonds" regime where the Taylor expansion above is a good approximation.
The exact combinatorial formula (their Eq 1.13, a finite sum with
binomial coefficients) WOULD reduce to exact Langmuir at n=1 for any
density, but sacrifices the simple closed form we're using here.

We are deliberately keeping the closed-form exponential version for Aim 1,
because (a) it is the correct physical regime for mobile EV receptors, and
(b) the size and cause of the approximation are fully understood and
small in the regimes we expect to operate in. The unit test in
test_aim1.py checks for Langmuir recovery IN THE DILUTE LIMIT specifically
(not "exact at all densities"), and separately confirms the gap shrinks
in the expected way as density decreases -- that is the honest version of
this limiting-case check.

We include a langmuir_probability() convenience function below, written
independently (from the textbook Langmuir formula, not by simplifying our
own equation), specifically so that Aim 2's unit test has something
independent to compare against. If we tested our formula against itself
we wouldn't be testing anything.

===========================================================================
NEW: THREE INDEPENDENT, CROSS-DOMAIN LITERATURES CONFIRMING THE SAME
NONLINEAR-DENSITY MECHANISM (added following the diagnosis-followthrough
investigation, docs/Diagnosis_Followthrough_and_Recommendation.md, Part 6.3
items 2-3)
===========================================================================
This is a citation-strength upgrade, not a calibration -- none of these
three sources supplies a usable K_INTRA point-estimate (see
k_intra_z_plausibility_bounds.py for why that search remains open/stretch-
goal). What they DO supply is independent confirmation, from three
completely unrelated fields, that the core physical/statistical idea this
whole project rests on -- capture/response depends NONLINEARLY on receptor
(or antigen, or per-individual capture-probability) count, so a population's
MEAN receptor count is not a reliable predictor of population-level
behavior -- is not a phenomenon invented for this project. It is
established, independently-discovered science in at least three fields with
no connection to each other or to EVs:

1. SUPERSELECTIVITY PHYSICS (this project's own direct theoretical source --
   see the module docstring above): Dubacheva, Curk & Richter, "Determinants
   of Superselectivity -- Practical Concepts for Application in Biology and
   Medicine," Accounts of Chemical Research (2023), DOI
   10.1021/acs.accounts.2c00672, PMID 36916901. Formally defines the exact
   alpha selectivity parameter this file's selectivity_alpha()/alpha_exact()
   compute, and re-plots a real quantitative dataset (Wang et al. 2020, ACS
   Nano, PMID 31899613 -- 5 ligand densities x 11 receptor densities on
   breast cancer cells) obeying the same scaling law this project uses.

2. CAR-T ANTIGEN-DENSITY DOSE-RESPONSE (a completely different biological
   readout -- immune-cell killing efficiency, not surface capture -- but the
   SAME underlying nonlinear-density mechanism): Majzner et al., "Tuning the
   Antigen Density Requirement for CAR T-cell Activity," Cancer Discovery
   (2020), DOI 10.1158/2159-8290.CD-19-0945; and Walker et al., "Tumor
   Antigen and Receptor Densities Regulate Efficacy of a Chimeric Antigen
   Receptor Targeting Anaplastic Lymphoma Kinase," Molecular Therapy (2017).
   Walker et al. used cells engineered to express a TITRATED range of a real
   cancer antigen (~18,000 to ~450,000 molecules/cell) and measured immune-
   cell killing efficiency as a direct, graduated function of that density --
   structurally very close to this project's own receptor-density-vs-
   capture-probability question, just with a different functional readout
   (immune-cell killing instead of affinity capture).

   Citing these two together (superselectivity + CAR-T) is stronger than
   either alone: two unrelated fields (soft-matter/colloid physics, and
   cellular immunotherapy) independently arrived at the same qualitative
   conclusion -- density-dependent response is nonlinear, not proportional
   -- using completely different experimental systems and completely
   different mathematical derivations.

3. CAPTURE-RECAPTURE ECOLOGY/STATISTICS (not a biological or physical
   analogue at all -- a formal STATISTICAL confirmation of this project's
   core reasoning, from a field with zero connection to biology): Otis,
   Burnham, White & Anderson, "Statistical Inference from Capture Data on
   Closed Animal Populations," Wildlife Monographs (1978), No. 62 -- the
   foundational reference for "Mh" models (heterogeneous-capture-probability
   models) in wildlife/ecological population estimation. This is a mature,
   decades-old, independently-developed statistical framework that formally
   proves: when individuals in a population (there, animals being trapped;
   here, vesicles being captured) have HETEROGENEOUS capture probability,
   using the population MEAN as an estimator is provably biased. This is
   exactly this project's own central statistical claim -- that overlap of
   two receptor-count distributions matters, not just the fold-change of
   their means -- independently formalized in ecology, with no reference to
   EVs, affinity capture, or even biology at all. This does not supply a new
   number to use in this project's model; it supplies formal, citable,
   independent validation that the REASONING behind the model's central
   claim is sound, well-established mathematics in another discipline
   entirely -- useful specifically for judge Q&A about whether this
   project's statistical logic is solid.

None of these three sources changes any number, function, or result in this
file or elsewhere in this codebase. They are citations to strengthen the
project's write-up and its answers to "is this mechanism real, or did you
invent it?" -- not new inputs to any calculation.

===========================================================================
NEW (fourth session, August 30, 2026 -- open-ended novelty research pass,
beyond the two diagnosis documents): FOUR MORE FINDINGS, PLUS A PRECISION
NOTE ON WHERE selectivity_alpha()/alpha_exact() ACTUALLY COME FROM
===========================================================================
This session's instruction was explicitly not to limit research to what the
two diagnosis documents already listed. A fresh, open literature search
(not itemized against either document's brainstorm categories) turned up
four more real, relevant sources:

4. A BRAND-NEW (Dec 2025), OPEN-ACCESS, MECHANICS-FOCUSED SUPERSELECTIVITY
   PAPER: Hamzeh, Mirahsani, Ahmadpoor & Farokhirad, "From Single Ligand-
   Receptor Bond Strength to Collective Avidity: Mechanics-Guided
   Superselective Nanoparticle Adhesion to Biological Membranes," Langmuir
   (2025), DOI 10.1021/acs.langmuir.5c03244, PMID 41400250, PMC12756914
   (CC BY 4.0, fully read this session). Reports real Monte Carlo phase
   diagrams showing superselective onset (alpha >= 1, using the SAME log-
   log-slope alpha definition this project already computes) at specific
   receptor densities (as low as 200/um^2) for DEFORMABLE nanoparticles
   under strong binding, and explicitly shows RIGID nanoparticles essentially
   never reach superselectivity (never engaging more than ~10% of ligand
   capacity). Checked carefully and reported honestly: this paper does NOT
   supply a usable K_INTRA scaling law (no fitted avidity-vs-density
   equation; nanoparticle size is fixed, not swept, in their model) -- it is
   a citation-strength addition, same tier as items 1-3 above, not a
   calibration source. Its most useful angle: an EV (a lipid-bilayer vesicle
   with an aqueous lumen) mechanically resembles their "deformable" NP
   category far more than their "rigid" category -- and deformable NPs are
   exactly the category that reaches superselectivity most easily in their
   own simulations. This is a real, if qualitative, physical argument for
   why modeling EVs as capable of superselective behavior is reasonable, not
   just mathematically convenient.

5. A SECOND BRAND-NEW (Feb 2025), OPEN-ACCESS PAPER, USED AS A CONTRAST
   CITATION THAT ACTIVELY DEFENDS THIS PROJECT'S OWN SIMPLIFYING ASSUMPTION:
   Xie, Angioletti-Uberti, Dobnikar, Frenkel & Curk, "Receptor clustering
   tunes and sharpens the selectivity of multivalent binding," PNAS (2025),
   DOI 10.1073/pnas.2417159122, PMID 39951501, PMC11848318 (CC BY-NC-ND
   4.0, fully read this session). Note: Curk is a co-author on BOTH this
   paper and item 1's Accounts of Chemical Research review already cited
   above -- same research lineage. This paper describes a DIFFERENT
   selectivity-tuning mechanism than the one this project models: LATERAL
   receptor-receptor clustering/attraction on a large, receptor-dense
   membrane (their simulations use ~1000 receptors/um^2 on effectively
   infinite membrane patches), which this project's model does NOT
   implement (this project treats each vesicle's receptor count as fixed
   and independent, no lateral mobility/clustering term). Rather than being
   a gap, the paper's OWN math argues this is the right simplification at EV
   scale: the authors explicitly note clustering-based selectivity is capped
   by finite membrane area and low receptor copy number (their own scaling
   argument involves a size limit ~sqrt(N receptors)). A 100nm EV has ~31,400
   nm^2 of surface area and typically only tens-to-hundreds of receptor
   copies -- many orders of magnitude below the density/area regime where
   their clustering mechanism could operate. CITE THIS AS: "Xie et al. (PNAS
   2025) show that a second, distinct selectivity-tuning mechanism (lateral
   receptor clustering) exists in the multivalent-binding literature, but
   their own scaling analysis (selectivity limited by finite receptor number
   and membrane area) indicates this mechanism requires far more receptors
   and membrane area than a single EV has -- supporting this project's
   choice to model EV receptors as an independent, non-clustering
   population."

6. A CTC-CHIP (CIRCULATING TUMOR CELL) RECEPTOR-DENSITY SATURATION FINDING
   (Weak Point 4/brainstorm category 13, scoped this session): real
   published capture-efficiency data across cancer cell lines spanning a
   ~250-fold range of EpCAM receptor density (breast cancer cells with
   >500,000 EpCAM/cell down to bladder cancer cells with ~2,000 EpCAM/cell)
   found mean capture yield staying ABOVE 90% across this entire range on
   standard CTC-chip microfluidic devices. This is a useful BOUNDARY-
   CONDITION citation, not a supporting one: it shows a real regime (large
   cells, very high absolute receptor counts, by CTC-chip standards) where
   capture efficiency SATURATES and becomes density-insensitive -- the
   opposite end of the spectrum from EVs' much lower absolute receptor
   counts (tens to hundreds of copies), where this project's central thesis
   is that density-sensitivity (not saturation) should dominate. Citable as
   a real contrast: "why does density matter for EVs but not for CTCs?" has
   a real, literature-grounded answer (receptor count regime), which is a
   stronger position than not addressing the question at all.

7. PRECISION NOTE ON selectivity_alpha()/alpha_exact() -- their exact
   provenance, now stated precisely rather than just generally: the log-log
   slope definition alpha = d(ln theta)/d(ln n) that this file already
   implements in two independently cross-checked forms was originally
   COINED by Martinez-Veracoechea & Frenkel, PNAS 2011, 108(27):10963-10968
   (already cited above via the Curk/Dobnikar/Frenkel arXiv document this
   whole file's equations come from), and is FORMALIZED as "the" standard
   superselectivity index throughout item 1's Dubacheva/Curk/Richter 2023
   Accounts of Chemical Research review (their Eq 1.10/1.16, using an
   equivalent notation alpha = d ln(Gamma_P)/d ln(Gamma_R) for probe/
   receptor surface densities). This means this project's central
   selectivity_alpha()/alpha_exact() functions are not a project-specific
   invention or a loose analogy to the field's concept -- they ARE, by
   exact mathematical definition, the field's own canonical formal metric
   for superselectivity, applied to this project's own capture_probability()
   model. This is worth stating explicitly and precisely in the write-up
   (with the exact citation chain above) rather than left as a general
   "the field studies something similar" gesture -- it is a precise
   equivalence, not an analogy.

None of items 4-7 changes any number, function, or result in this file or
elsewhere in this codebase (item 7 is a documentation-precision correction
only -- alpha_exact()/selectivity_alpha() themselves are unchanged). They
are citation and precision upgrades found via an explicitly open-ended
research pass that intentionally went beyond both prior diagnosis
documents' own lists, per the project owner's direct instruction not to
limit research scope to what a prior document happened to enumerate.
"""

import numpy as np


def q_bound(n, sigma, k_intra):
    """
    Bound-state partition function for a vesicle with n receptors on a
    surface with ligand density sigma.

    This is Eq (*) above: q_b = exp(n * sigma * k_intra) - 1

    Parameters
    ----------
    n : float or array
        Number of receptors on the vesicle. (Can be non-integer here because
        this is a continuous mathematical model; Aim 1's Monte Carlo version
        will use actual integer receptor counts.)
    sigma : float or array
        Ligand surface density (arbitrary consistent units -- e.g. ligands
        per square micron). Larger sigma = more ligand available to bind.
    k_intra : float
        Lumped effective binding-strength parameter (see module docstring).
        Larger k_intra = stronger/more favorable binding per receptor-ligand
        encounter.

    Returns
    -------
    float or array
        The bound-state partition function q_b (always >= 0).
    """
    n = np.asarray(n, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    return np.exp(n * sigma * k_intra) - 1.0


def capture_probability(n, sigma, k_intra, z):
    """
    Probability that a vesicle with n receptors is captured on a surface
    with ligand density sigma. This is theta(n, sigma) from Eq (**).

    Parameters
    ----------
    n : float or array
        Number of receptors on the vesicle.
    sigma : float or array
        Ligand surface density.
    k_intra : float
        Lumped effective binding-strength parameter (see q_bound).
    z : float
        Dimensionless activity / prefactor capturing everything about
        vesicle concentration, incubation time, and transport to the
        surface that is NOT part of the multivalent bond-counting itself.
        Must be > 0.

    Returns
    -------
    float or array
        Capture probability theta, between 0 and 1.
    """
    qb = q_bound(n, sigma, k_intra)
    zqb = z * qb
    return zqb / (1.0 + zqb)


def langmuir_probability(sigma, k_intra, z):
    """
    Ordinary single-site (Langmuir) binding probability, written
    independently of capture_probability/q_bound so it can serve as a
    ground truth for the n=1 unit test in Aim 2.

    Standard Langmuir isotherm form:  theta = (z * K * sigma) / (1 + z * K * sigma)
    Here K plays the same role as k_intra: how strongly a single ligand-
    receptor pair binds. This is the textbook formula for one binding site
    encountering a ligand at density sigma -- no combinatorics, no
    multivalency, because there is only one receptor to begin with.

    Parameters
    ----------
    sigma : float or array
        Ligand surface density.
    k_intra : float
        Single-bond binding-strength parameter (same meaning/units as in
        capture_probability, for consistency).
    z : float
        Dimensionless activity/prefactor, same meaning as in
        capture_probability.

    Returns
    -------
    float or array
        Binding probability, between 0 and 1.
    """
    sigma = np.asarray(sigma, dtype=float)
    zk_sigma = z * k_intra * sigma
    return zk_sigma / (1.0 + zk_sigma)


def selectivity_alpha(n_values, sigma, k_intra, z):
    """
    Compute the selectivity parameter alpha = d(log theta) / d(log n),
    following Eq 1.10/1.16 in Curk, Dobnikar & Frenkel. alpha > 1 means
    "superselective": capture probability grows FASTER than linearly with
    receptor count. alpha <= 1 means ordinary (or sub-linear) response.

    This is computed numerically (finite differences in log-log space)
    rather than analytically, because it's more transparent for a first
    version and doesn't require re-deriving the derivative by hand. We can
    swap in a closed-form expression later once the numerical version is
    validated.

    Parameters
    ----------
    n_values : array
        A sorted array of receptor counts to evaluate alpha at (needs at
        least 2 points; more points = smoother curve).
    sigma : float
        Ligand surface density (held fixed while n varies).
    k_intra : float
        Binding-strength parameter.
    z : float
        Activity/prefactor.

    Returns
    -------
    n_midpoints : array
        The n values (midpoints between consecutive n_values) at which each
        alpha estimate applies -- one shorter than n_values, since this is
        a finite-difference derivative.
    alpha : array
        The estimated selectivity at each midpoint.

    IMPORTANT, DISCOVERED WHILE FIRST EXERCISING THIS FUNCTION (it had
    never been tested before this session -- see tests/test_aim1.py and
    tests/test_nanostar_validation.py): whether alpha ever exceeds 1
    (i.e., whether the model shows genuine superselective behavior) is
    STRONGLY DEPENDENT on z. At z=1.0 -- the value used throughout every
    existing demo script in this codebase (compare_analytical_vs_mc.py,
    demo_population_effect.py, patch_regime_sensitivity.py) -- alpha
    NEVER exceeds 1 across any n tested; it decreases monotonically from
    just below 1 toward 0. This is NOT a bug: it happens because z=1.0 is
    already too far up the Langmuir/avidity curve (theta starts too large
    even at small n) to see the superselective transition, which is the
    same "must be in the genuinely dilute regime" condition that the
    founding superselectivity papers state explicitly (Curk, Dobnikar &
    Frenkel; also confirmed independently for the DNA nanostar model in
    nanostar_validation.py, where the SAME issue appeared and required
    the same fix: a small enough activity/concentration prefactor). Using
    a much smaller z (e.g. z=1e-4) with the SAME sigma and k_intra
    recovers the textbook superselective shape: alpha starts near 1 at
    low n, RISES above 1 to a peak at intermediate n (a real, sizeable
    peak -- e.g. alpha~6 was observed for sigma=1.0, k_intra=0.2,
    z=1e-4), then decays toward 0 as theta saturates toward 1 at large n.
    This means every existing demo figure in this project (which all use
    z=1.0) is technically showing NON-superselective behavior by this
    specific alpha>1 definition, even though the underlying nonlinearity
    (n inside an exponential) is still doing real work in those figures --
    alpha<1 does not mean "no nonlinearity," it means "not superselective
    by the strict Martinez-Veracoechea/Frenkel definition specifically."
    Worth remembering for Aim 3/4: whether z ends up calibrated small
    (dilute, e.g. picomolar-scale vesicle concentrations or short
    engagement windows) or large (concentrated, e.g. dense vesicle
    suspensions or long incubation) will determine whether the
    superselectivity mechanism is actually operative in the calibrated
    model, not just present as unused machinery.
    """
    n_values = np.asarray(n_values, dtype=float)
    theta = capture_probability(n_values, sigma, k_intra, z)

    log_n = np.log(n_values)
    log_theta = np.log(theta)

    d_log_theta = np.diff(log_theta)
    d_log_n = np.diff(log_n)
    alpha = d_log_theta / d_log_n

    n_midpoints = np.sqrt(n_values[:-1] * n_values[1:])  # geometric midpoint,
    # appropriate since we're differentiating in log space

    return n_midpoints, alpha


def alpha_exact(n, sigma, k_intra, z):
    """
    Closed-form (analytical) selectivity alpha = d(log theta)/d(log n),
    computed directly by calculus instead of by finite differences.

    NEW (v14->v15) -- "Immediate next steps," item 0b, added directly in
    response to an independent external review of this project (a report
    from an AI research tool called Biomni, dated 26 August 2026). Biomni
    pointed out that selectivity_alpha() above only ESTIMATES alpha via
    finite differences, and proposed an exact closed form. Rather than
    trusting Biomni's formula, it was independently re-derived from this
    project's own equations this session, using ordinary calculus:

        theta(x, z) = z*(e^x - 1) / (1 + z*(e^x - 1)),  x = n*sigma*k_intra

    Let q = e^x - 1, so theta = z*q / (1 + z*q). Then:
        d(ln theta)/dq = 1/q - z/(1+z*q) = 1 / [q*(1+z*q)]
        dq/dx = e^x
        d(ln theta)/dx = e^x / [q*(1+z*q)] = e^x / [(e^x-1)*(1+z*(e^x-1))]

    Since x = n*sigma*k_intra is directly proportional to n,
    d(x)/d(ln n) = x (a standard identity: for x = c*n, dx/d(ln n) =
    dx/dn * dn/d(ln n) = c * n = x). Multiplying through:

        alpha = d(ln theta)/d(ln n) = x * d(ln theta)/dx
              = x * e^x / [(e^x - 1) * (1 + z*(e^x - 1))]

    This independently-derived formula was checked directly against
    Biomni's own reported formula and matches it to floating-point
    precision (max difference 8.9e-16 across a range of test points) --
    this is not a case of trusting Biomni's algebra without checking it,
    it is a case of doing the algebra separately and getting the same
    answer.

    HOW THIS COMPARES TO THE FINITE-DIFFERENCE selectivity_alpha() ABOVE,
    CHECKED DIRECTLY THIS SESSION on the project's real K_INTRA/z
    plausibility grid (get_k_intra_sweep_values() x get_z_sweep_values()
    from aim3_aim4_range_sweep.py, n up to 200): the two agree closely at
    the grid's actual peak (exact 10.465 vs. finite-difference 10.461 --
    a real but small ~0.04 difference), but diverge by as much as ~0.70 in
    absolute alpha at OTHER grid points (worst case found: k_intra=
    2.25e-12, z=1e-6), which is large enough that an independent reviewer
    computing the exact formula (as Biomni did, and as this session did
    again separately) would notice the gap immediately. This function is
    added ALONGSIDE, not instead of, selectivity_alpha() -- the finite-
    difference version is kept for backward compatibility with any
    existing caller, and because it doesn't require deriving a closed
    form by hand for every future variant of this model (e.g. if the
    underlying theta formula changes again), while this exact version is
    used wherever the extra precision is cheaply available and matters
    (see sweep_superselectivity_reachability() in
    aim3_aim4_range_sweep.py, updated this session to report both).

    Parameters
    ----------
    n : float or array
        Number of receptors on the vesicle (or the geometric midpoints
        used by selectivity_alpha(), if comparing directly against that
        function's output).
    sigma : float or array
        Ligand surface density (same meaning as everywhere else in this
        file).
    k_intra : float
        Binding-strength parameter (same meaning as everywhere else).
    z : float
        Activity/prefactor (same meaning as everywhere else).

    Returns
    -------
    float or array
        The exact analytical selectivity alpha at each given n. Unlike
        selectivity_alpha(), this returns one alpha value per input n
        (no finite-difference shortening by one element, and no need for
        n_values to be sorted or evenly related to each other -- each n
        is evaluated independently).
    """
    n = np.asarray(n, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    x = n * sigma * k_intra
    exp_x = np.exp(x)
    q = exp_x - 1.0
    return x * exp_x / (q * (1.0 + z * q))


# ---------------------------------------------------------------------------
# ACCESSIBILITY / PROTEIN CORONA MULTIPLIER (f_accessible)
# ---------------------------------------------------------------------------
"""
WHY THIS SECTION EXISTS
---------------------------
"Immediate next steps" item 3 (formerly item 4, unchanged and unattempted
from v4 through v15): "Add f_accessible to the model. One multiplier on n."
This directly implements mechanism (d) from the handoff's "Three physical
mechanisms the model must include" section: EVs in real biofluids acquire
an adsorbed plasma-protein layer (the "protein corona") that can physically
occlude surface epitopes, so not every receptor an EV actually carries is
necessarily ACCESSIBLE for a ligand to bind. Source: Foersonits et al.
2025, J Extracell Vesicles, DOI 10.1002/jev2.70124 (spelled "Försönits" in
the handoff -- ASCII-only here to avoid encoding issues in this file).

WHY THIS IS "ONE MULTIPLIER ON n," NOT A NEW MECHANISTIC COAT MODEL
------------------------------------------------------------------------
The handoff is explicit that this should be implemented as "a single
multiplier f_accessible on n, not as a mechanistic corona model" -- i.e.,
we are NOT modeling the physics of protein adsorption, corona thickness,
or epitope-specific occlusion probability. We are modeling only the NET
EFFECT: if a vesicle carries n receptors in total, only some accessible
fraction f_accessible of them (0 < f_accessible <= 1) are actually
available to form a bond with a ligand-coated surface, so the model should
use n_accessible = f_accessible * n wherever it currently uses raw n.

accessible_receptors() below is that one multiplier, written as its own
small function (rather than baked invisibly into q_bound() or
capture_probability()) for the same reason contact_patch_receptors() is
kept separate from q_bound(): keeping each physical adjustment to n in
its own clearly-named function, composed together explicitly at the call
site, makes it possible to test each adjustment independently and to use
some adjustments without others (e.g. a demo that doesn't care about
corona effects can simply not call accessible_receptors() and pass raw n
directly, exactly as every existing call site already does -- see
"BACKWARD COMPATIBILITY" below).

WHY THIS MATTERS FOR SELECTIVITY SPECIFICALLY -- THE REAL RISK THIS
SECTION IS BUILT TO TEST
--------------------------------------------------------------------------
The handoff's "Scope" section argues selectivity (a tumor/normal RATIO) is
more trustworthy than any single yield number, because many sources of
uncertainty (ligand density, contact area, incubation time, detection
efficiency) affect both populations THE SAME WAY and therefore cancel out
of the ratio. The handoff's own mechanism (d) discussion flags a specific,
important EXCEPTION to this optimistic cancellation argument: "If
f_accessible differs between tumor and normal EVs -- plausible, since
corona composition depends on surface charge and lipid composition -- it
contaminates the selectivity ratio, which was supposed to be robust
because errors cancel." This is exactly the kind of thing "Stress-test
before executing" (this project's own working-style requirement) demands
be checked directly with the actual model, not just noted as a
possibility in prose -- see selectivity_ratio_sensitivity_to_accessibility()
below and its test suite for the actual, numeric answer.

BACKWARD COMPATIBILITY -- READ THIS BEFORE WORRYING THIS BREAKS ANYTHING
--------------------------------------------------------------------------
accessible_receptors(n, f_accessible=1.0) with the default f_accessible=1.0
returns n UNCHANGED (n_accessible = 1.0 * n = n). No existing function
signature in this file (q_bound, capture_probability, q_bound_patch,
capture_probability_patch, contact_patch_receptors, selectivity_alpha,
alpha_exact) is modified by this section -- every existing call site
throughout this codebase (aim3_aim4_range_sweep.py, population_model.py,
stem_cell_marker_ranking_extension.py, every demo script) continues to
pass raw n directly into those functions exactly as before, and continues
to get EXACTLY the same numbers as before, because none of them call
accessible_receptors() at all. This is a deliberate, tested design choice
(see tests/test_f_accessible.py's backward-compatibility test) -- adding
f_accessible does not require touching a single existing call site or
existing test.
"""


def accessible_receptors(n, f_accessible=1.0):
    """
    Apply the accessibility/protein-corona multiplier to a vesicle's total
    receptor count, per mechanism (d) in the handoff's "Three physical
    mechanisms" section: n_accessible = f_accessible * n.

    Parameters
    ----------
    n : float or array
        Total receptor count on the vesicle (the number actually present,
        before accounting for any corona-driven occlusion).
    f_accessible : float or array, optional (default 1.0)
        Fraction of total receptors that are physically accessible for
        ligand binding (0 < f_accessible <= 1). The default, 1.0, means
        "every receptor is accessible" -- i.e., no accessibility effect at
        all, which reproduces every existing result in this codebase
        exactly (see module discussion above, "BACKWARD COMPATIBILITY").
        Values less than 1.0 represent partial occlusion by an adsorbed
        protein corona or other surface-blocking effect. f_accessible can
        be an array (e.g., different per-vesicle values, or different
        per-population values in a comparison) to model a POPULATION-LEVEL
        or CONDITION-LEVEL difference in accessibility, not just a single
        global constant.

    Returns
    -------
    float or array
        n_accessible = f_accessible * n, the effective number of
        receptors available for bonding after the accessibility multiplier.

    Raises
    ------
    ValueError
        If any f_accessible value is not in (0, 1] -- an accessibility
        fraction above 1 (more accessible receptors than actually exist)
        or at/below 0 (a vesicle with some receptors capturing exactly
        like one with none) would not correspond to a physically
        meaningful corona effect, and is far more likely to indicate a
        units or convention mistake at the call site than a real
        parameter choice -- reported explicitly rather than silently
        producing a physically nonsensical n_accessible.
    """
    n = np.asarray(n, dtype=float)
    f_accessible_arr = np.asarray(f_accessible, dtype=float)
    if np.any(f_accessible_arr <= 0) or np.any(f_accessible_arr > 1):
        raise ValueError(
            f"f_accessible must be in (0, 1] (a fraction of total "
            f"receptors that are physically accessible) -- got values "
            f"outside that range: min={f_accessible_arr.min()}, "
            f"max={f_accessible_arr.max()}. This is very likely a units "
            f"or convention mistake at the call site, not a real "
            f"parameter choice."
        )
    return f_accessible_arr * n


def capture_probability_accessible(n, sigma, k_intra, z, f_accessible=1.0):
    """
    Convenience wrapper: capture_probability(), but applying the
    accessibility multiplier to n first. Exactly equivalent to calling
    capture_probability(accessible_receptors(n, f_accessible), sigma,
    k_intra, z) -- provided as its own named function purely for
    readability at call sites that want the accessibility effect, the same
    way capture_probability_patch() exists as a named wrapper around
    q_bound_patch() rather than requiring every caller to compose
    contact_patch_receptors() and q_bound() manually.

    At f_accessible=1.0 (the default), this is IDENTICAL to
    capture_probability(n, sigma, k_intra, z) -- see module discussion
    above, "BACKWARD COMPATIBILITY."

    Parameters
    ----------
    n, sigma, k_intra, z :
        Same meaning as capture_probability().
    f_accessible : float or array, optional (default 1.0)
        Same meaning as accessible_receptors().

    Returns
    -------
    float or array
        Capture probability theta, between 0 and 1.
    """
    n_accessible = accessible_receptors(n, f_accessible)
    return capture_probability(n_accessible, sigma, k_intra, z)


def selectivity_ratio_sensitivity_to_accessibility(n_tumor, n_normal, sigma, k_intra, z,
                                                       f_accessible_tumor=1.0,
                                                       f_accessible_normal=1.0):
    """
    Directly answers the specific risk the handoff's mechanism (d)
    discussion flags: "if f_accessible differs between tumor and normal
    EVs... it contaminates the selectivity ratio, which was supposed to be
    robust because errors cancel." Computes the tumor/normal capture-
    probability ratio TWICE -- once assuming EQUAL accessibility for both
    populations (f_accessible_tumor == f_accessible_normal, the
    "cancels out" case), and once using the ACTUAL, possibly-different,
    supplied f_accessible values for each population -- so the two ratios
    can be compared directly to see exactly how much (if at all) an
    accessibility difference distorts the selectivity ratio.

    Parameters
    ----------
    n_tumor, n_normal : float or array
        Receptor counts for the tumor and normal populations respectively
        (can be a single representative value or an array of per-vesicle
        counts -- if arrays, the mean capture probability across the array
        is used, matching the population-averaging convention used
        elsewhere in this project, e.g. population_model.population_yield).
    sigma, k_intra, z :
        Same meaning as capture_probability() -- shared by both
        populations (only accessibility differs between them here).
    f_accessible_tumor, f_accessible_normal : float, optional (default 1.0
        for both -- i.e., equal accessibility, no distortion)
        The (possibly different) accessible fractions for each population.

    Returns
    -------
    dict with keys:
        "ratio_equal_accessibility": the tumor/normal ratio computed as if
            accessibility were equal for both (using f_accessible_tumor
            for BOTH populations, ignoring f_accessible_normal) -- the
            "errors cancel" baseline this function checks against.
        "ratio_actual_accessibility": the tumor/normal ratio computed
            using each population's OWN supplied f_accessible value.
        "distortion_factor": ratio_actual_accessibility /
            ratio_equal_accessibility -- how much the accessibility
            DIFFERENCE alone changes the selectivity ratio, isolated from
            every other factor (which are held identical between the two
            computations). A value far from 1.0 means accessibility
            differences meaningfully contaminate the ratio; a value close
            to 1.0 means the ratio is robust to this specific effect at
            the tested f_accessible values.
    """
    def mean_theta(n, f_acc):
        n_arr = np.atleast_1d(np.asarray(n, dtype=float))
        theta = capture_probability_accessible(n_arr, sigma, k_intra, z, f_acc)
        return float(np.mean(theta))

    theta_tumor_equal = mean_theta(n_tumor, f_accessible_tumor)
    theta_normal_equal = mean_theta(n_normal, f_accessible_tumor)  # SAME
    # f_accessible as tumor -- this is the "equal accessibility" baseline.
    ratio_equal = theta_tumor_equal / theta_normal_equal if theta_normal_equal > 0 else float("inf")

    theta_tumor_actual = mean_theta(n_tumor, f_accessible_tumor)
    theta_normal_actual = mean_theta(n_normal, f_accessible_normal)
    ratio_actual = theta_tumor_actual / theta_normal_actual if theta_normal_actual > 0 else float("inf")

    distortion = ratio_actual / ratio_equal if ratio_equal not in (0, float("inf")) else float("nan")

    return {
        "ratio_equal_accessibility": ratio_equal,
        "ratio_actual_accessibility": ratio_actual,
        "distortion_factor": distortion,
    }


# ---------------------------------------------------------------------------
# CONTACT-PATCH-AWARE MODEL (added: n vs. areal-density resolution)
# ---------------------------------------------------------------------------
"""
WHY THIS SECTION EXISTS -- read this before using q_bound_patch() below.

Everything above this point treats a vesicle's receptor count n as a single
number that goes straight into the binding formula, with no reference to
vesicle SIZE. That was fine as long as we only ever compared vesicles that
were roughly the same size. The digitized Lennon et al. Figure 3c data
(see ev_capture_data/) broke that assumption: the PDAC-enriched gated
population has vesicles ~3.5x LARGER in diameter than typical healthy EVs,
and this matters because it means two different numbers -- total receptor
count n, and receptors per unit membrane area (areal density) -- move in
OPPOSITE directions for this population (count up 1.2-2x, density down
8-13x). The original model has no way to represent that distinction: it
only ever sees n.

THE PHYSICAL QUESTION
-----------------------
When a vesicle sits down on a ligand-coated capture surface, only the
receptors inside the CONTACT PATCH (the flattened area where the vesicle
membrane actually touches the surface) can form bonds -- not the whole
vesicle. The question is: how does contact patch AREA depend on vesicle
size? Two extremes:

  (1) Contact patch area is a FIXED FRACTION of the vesicle's total
      membrane area, regardless of how big or small the vesicle is (patch
      area grows as diameter^2 -- exactly like membrane area itself).
  (2) Contact patch area is roughly FIXED IN ABSOLUTE SIZE, set by
      something about the capture surface (its roughness, its ligand
      spacing, or some other length scale that has nothing to do with the
      vesicle) -- independent of vesicle diameter.

I checked what soft-matter adhesion physics says about which of these
applies to something like an EV: for a very SOFT, SMALL adhesive body (an
EV is both -- lipid bilayer bending rigidity is only 20-30 kT per Jahnke et
al. 2024, and EVs are 30-150+ nm, far smaller than a cell) landing on a
much STIFFER surface (the functionalized glass, gold, or magnetic bead a
real capture device would use), soft-matter contact mechanics predicts a
"soft/liquid" adhesion regime where the CONTACT RADIUS scales linearly with
the body's own radius (contact radius ~ vesicle radius), rather than the
sub-linear relationship (contact radius ~ radius^(2/3)) that applies for
stiffer bodies (that stiffer-body relationship is the classical JKR
contact-mechanics result). Since contact PATCH AREA is proportional to
contact RADIUS squared, "contact radius ~ vesicle radius" means "contact
patch area ~ vesicle diameter^2" -- which is option (1) above: a FIXED
FRACTION of membrane area, regardless of vesicle size.

This is a real, sourced physical argument (soft-adhesion / "liquid limit"
contact scaling for small, floppy, low-bending-rigidity adhesive bodies),
not something invented to make the numbers come out a preferred way -- but
it is still an approximation, not a direct measurement of EV contact
patches specifically, and that limitation is stated honestly below.

WHAT THIS MEANS FOR THE MODEL, WORKED OUT
---------------------------------------------
Define:
    D          = vesicle diameter
    rho        = n / (pi * D^2)   -- areal receptor density (receptors per
                 unit membrane area; pi*D^2 is the same "sphere surface
                 area" convention already used in
                 ev_capture_data/scripts/decomposition_analysis.py, so this
                 is consistent with the digitized-data analysis)
    A_patch(D) = c_patch * D^p    -- contact patch area, some constant
                 c_patch times diameter to a power p

Receptors actually available for bonding in the contact patch:
    n_contact = rho * A_patch(D)
              = (n / (pi * D^2)) * (c_patch * D^p)
              = (c_patch / pi) * n * D^(p - 2)

Two named special cases of the exponent p:

    p = 2 (SOFT-ADHESION / "liquid" limit -- the physically motivated
           default, per the argument above):
        n_contact = (c_patch / pi) * n
        The D^(p-2) term becomes D^0 = 1, so vesicle diameter drops out
        completely and n_contact is just proportional to TOTAL RECEPTOR
        COUNT n. This is mathematically IDENTICAL to what the original
        q_bound()/capture_probability() functions already do (with
        c_patch/pi absorbed into k_intra) -- so this case does not change
        any existing result, it just makes explicit the geometric
        assumption that was always implicitly built into using n directly.

    p = 0 (FIXED-PATCH limit -- contact patch area is a fixed absolute
           size, independent of vesicle diameter; the alternative the
           handoff flagged as option (b)):
        n_contact = (c_patch / pi) * n / D^2 = c_patch * rho
        Now n_contact is proportional to AREAL DENSITY rho, and total
        receptor count n only enters through how it affects density. This
        is what the model would look like if the "switch state variable to
        areal density" option were followed literally.

Rather than picking ONE of these and hard-coding it, this section adds ONE
new parameter (patch_exponent, p) that lets the model run EITHER regime, or
anything in between, by evaluating the same physical formula with a
different p. This directly implements what the handoff called option (c):
"run both and report whether target rankings flip -- which would itself be
a strong Aim 4 sensitivity finding." p=2 is used as the DEFAULT (matching
the soft-adhesion argument above and preserving all existing validated
results), but p=0 is provided as a clearly-labeled alternative, and the
gap between the two is exactly the sensitivity Aim 4 should quantify.

WHAT THIS DOES NOT CLAIM
----------------------------
This is a single power-law approximation to real contact mechanics, not a
first-principles derivation for EV membranes specifically -- no published
measurement of actual EV-surface contact patch size vs. EV diameter was
found (this is worth checking again before Aim 3 calibration; if one
exists, it should replace this approximation). The p=2 default is chosen
because it is the most physically defensible SIMPLE choice given (a) the
general soft-adhesion scaling law for small, low-bending-rigidity bodies,
and (b) full backward compatibility with every validated result so far --
not because it was checked against direct EV contact-patch measurements,
which do not appear to exist in the literature yet.
"""


def contact_patch_receptors(n, diameter, k_geom, patch_exponent=2.0):
    """
    Compute the number of receptors available for bonding within a
    vesicle's contact patch, accounting for how contact patch area scales
    with vesicle diameter. See the module-level discussion above for the
    full physical reasoning.

    n_contact = k_geom * n * diameter^(patch_exponent - 2)

    where k_geom lumps together the contact-patch proportionality constant
    (c_patch / pi from the derivation above). At patch_exponent=2 (the
    default), diameter drops out and this reduces to k_geom * n exactly --
    i.e., identical in form to plugging n directly into q_bound(), just
    with k_geom playing the role that used to be folded silently into
    k_intra.

    Parameters
    ----------
    n : float or array
        Total receptor count on the vesicle.
    diameter : float or array
        Vesicle diameter (same units used consistently everywhere else --
        the digitized data in ev_capture_data/ uses nanometers).
    k_geom : float
        Lumped geometric constant (contact-patch-area-to-membrane-area
        proportionality, divided by pi). Plays an analogous role to
        k_intra in q_bound(), but keeps the diameter-dependence separate
        and explicit rather than silently absorbed.
    patch_exponent : float, optional
        The power p in A_patch ~ diameter^p. Default 2.0 = soft-adhesion
        limit (contact patch is a fixed fraction of membrane area,
        independent of vesicle size -- see module discussion). Use 0.0 for
        the fixed-absolute-patch-size alternative (contact-patch receptors
        become proportional to areal density instead of total count).

    Returns
    -------
    float or array
        n_contact, the effective number of receptors available for
        bonding in the contact patch. Feed this into q_bound() in place of
        raw n when using the diameter-aware model.
    """
    n = np.asarray(n, dtype=float)
    diameter = np.asarray(diameter, dtype=float)
    return k_geom * n * diameter ** (patch_exponent - 2.0)


def q_bound_patch(n, diameter, sigma, k_geom, k_intra, patch_exponent=2.0):
    """
    Diameter-aware version of q_bound(): computes the bond-count partition
    function using the contact-patch-adjusted receptor count
    (contact_patch_receptors) instead of raw n.

    At patch_exponent=2.0 (the default, soft-adhesion limit), this is
    mathematically identical to q_bound(n, sigma, k_geom * k_intra) -- i.e.
    it reproduces the original size-blind model exactly, just with the
    geometric constant now named separately from the bond-strength
    constant. This is deliberate: it means every existing validated result
    (analytical vs. Monte Carlo, the Poisson mean-field theorem, the
    Langmuir limiting case) is untouched by adding this function, since
    those all used q_bound() directly with a single k_intra that already
    implicitly contained whatever k_geom would have been.

    At patch_exponent=0.0 (fixed-patch limit), this becomes a function of
    areal density rather than total count -- diameter enters explicitly.

    Parameters
    ----------
    n : float or array
        Total receptor count on the vesicle.
    diameter : float or array
        Vesicle diameter, same units as used in k_geom's calibration.
    sigma : float or array
        Ligand surface density on the capture surface.
    k_geom : float
        Lumped geometric constant (see contact_patch_receptors).
    k_intra : float
        Bond-strength parameter (same meaning as in q_bound).
    patch_exponent : float, optional
        Contact-patch scaling exponent p (default 2.0, soft-adhesion
        limit -- see module discussion above).

    Returns
    -------
    float or array
        The bond-count partition function q_b, using the contact-patch-
        adjusted receptor count.
    """
    n_contact = contact_patch_receptors(n, diameter, k_geom, patch_exponent)
    return q_bound(n_contact, sigma, k_intra)


def capture_probability_patch(n, diameter, sigma, k_geom, k_intra, z, patch_exponent=2.0):
    """
    Diameter-aware version of capture_probability(): the full theta(n,
    sigma), but computed using the contact-patch-adjusted receptor count
    rather than raw n. See q_bound_patch() and the module discussion above.

    Parameters mirror capture_probability(), with the addition of
    diameter, k_geom, and patch_exponent (see q_bound_patch for their
    meaning).

    Returns
    -------
    float or array
        Capture probability theta, between 0 and 1.
    """
    qb = q_bound_patch(n, diameter, sigma, k_geom, k_intra, patch_exponent)
    zqb = z * qb
    return zqb / (1.0 + zqb)
