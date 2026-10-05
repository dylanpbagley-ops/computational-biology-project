"""
k_intra_z_plausibility_bounds.py

WHY THIS FILE EXISTS
---------------------
This is NOT a calibration script and NOT a fitted model. It is the
documented FALLBACK that the v6 handoff named explicitly: "if no clean
source is found, consider whether K_INTRA/z/k_geom can instead be
constrained by physical plausibility bounds (e.g., known antibody-antigen
Kd ranges) rather than fit to a specific number -- worth discussing as a
fallback, not silently defaulting to it." This session searched hard for a
genuine single-marker ABSOLUTE calibration source (see the handoff's
"Immediate next steps" item 1 and the search log below) and did not find
one clean enough to fit. The user was asked directly (via AskUserQuestion,
matching the "stress-test before executing" working style) and chose this
fallback over forcing an imperfect fit or continuing to search indefinitely.

WHAT WAS SEARCHED, AND WHY EACH CANDIDATE DIDN'T WORK AS A CALIBRATION
------------------------------------------------------------------------
Three real candidate sources were found and evaluated. None cleanly
provides "absolute single-marker capture efficiency at a known ligand
surface density" -- which is literally what K_INTRA/z/k_geom need.

1. Mizenko et al. 2021's OWN raw Figure 2e serial-dilution data (re-read
   directly this session, not just the summary used for the (rejected)
   cross-marker co-occurrence calibration attempt from the v5->v6
   transition -- see binding_model.py's other docstring section for that
   separate, already-rejected attempt). This IS a genuine single-marker
   (anti-CD81 capture / anti-CD81 detect) dose-response: EVs at known
   NTA-measured bulk concentrations (4x10^7 to 4x10^9 particles/mL) were
   incubated on a fixed antibody spot, giving captured-particle counts per
   spot. No cross-marker ratio problem this time. BUT: the ExoView chip's
   antibody spotting density (sigma in our model's units) is a proprietary
   commercial-platform detail that is not published anywhere findable
   (checked: paper's own methods, manufacturer specs, product literature).
   Without sigma, there is no way to separate "how strong is the bond"
   (K_INTRA) from "how much ligand is there" (sigma) -- only their
   product is observable, and even that requires assuming a receptor
   copy-number distribution for these particular EVs, which Mizenko's
   paper doesn't report for CD81 specifically outside the tetraspanin
   panel context.

2. Lopez Baltazar et al. 2024 (ACS Sensors, DOI 10.1021/acssensors.4c00558)
   and the underlying 2024 Cornell PhD dissertation (Lopez Baltazar,
   "Direct Detection of Extracellular Vesicle Transmembrane Protein -
   Antibody Interactions Using Surface Plasmon Resonance Biosensors",
   DOI 10.7298/bssz-8415, open access via Cornell eCommons,
   https://hdl.handle.net/1813/116519). This is the CLOSEST physical match
   found: a real surface-immobilized (SPR chip) anti-CD81 antibody assay,
   with a directly measured equilibrium dissociation constant Kd = 1.92 nM
   for EV binding, at "40% surface coverage" of polyclonal antibody on a
   mixed self-assembled monolayer, with a reported linear dynamic range of
   1.9x10^8-1.9x10^9 EVs/mL.

   UPDATE (v7->v8 handoff transition): the dissertation's full text WAS
   retrieved this session (via Claude-in-Chrome, which was connected this
   time -- see the handoff for the retrieval method: Google Docs' PDF
   viewer proxying the eCommons bitstream URL, since Chrome's own PDF
   handling triggered a download rather than an inline view). Appendix
   A.1 ("Antibody Surface Coverage Estimation," p.136) gives the actual
   antibody surface density in physical units that v7 could not recover:
   the SPR chip's MAXIMUM close-packed antibody density is 1.11x10^12
   Abs/cm^2 (derived from real antibody dimensions, not assumed), so the
   "40% surface coverage" condition that produced Kd=1.92 nM corresponds
   to a real, physical sigma = 0.40 x 1.11x10^12 = 4.44x10^11 Abs/cm^2
   (equivalent to 104 ng/cm^2, also stated directly in the text). This
   closes the SPECIFIC gap v7 named ("only '40% coverage of an unspecified
   maximum'"). See SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2 below.

   HOWEVER -- and this matters, checked carefully rather than assumed --
   having a real sigma alongside a real Kd does NOT, on its own, unlock a
   clean K_INTRA fit, for a reason distinct from (but structurally similar
   to) why Kakavoulia's FIDA data was rejected. The SPR Kd=1.92nM is a 3D
   BULK EV CONCENTRATION (the EV particle concentration in the flowing
   sample at which surface capture reaches half-maximum), not a 2D areal
   binding constant. It was measured across the natural, uncontrolled
   RAW 264.7 EV population -- i.e., it is an ensemble average over
   whatever CD81 copy-number distribution those EVs actually have, and it
   also folds in mass-transport/diffusion to the surface, exactly the kind
   of "everything other than the multivalent bond-counting itself" that
   this project's own z parameter is defined to absorb (see
   binding_model.py's capture_probability docstring). So Kd, properly
   understood, calibrates z's CONCENTRATION axis (see point 3 below, now
   updated with real linear-range and control-experiment concentrations
   from the same dissertation), not K_INTRA directly. Converting a 3D
   bulk Kd into a 2D single-bond K_INTRA would still require an
   encounter-volume/contact-geometry model for EVs that does not exist --
   the same missing piece already documented in point 1 of this section
   for K_INTRA generally, and in binding_model.py's contact-patch
   discussion for k_geom.

   WHAT THIS SESSION'S RETRIEVAL DOES CONCRETELY IMPROVE: sigma is no
   longer an unknown "commercial black box" number for this source (it
   was for Mizenko's ExoView chip, and remains so) -- it is now a real,
   citable, physically-derived antibody areal density. This lets the
   dimensionless-exponent bound (point 2 below) be divided through by a
   REAL sigma instead of left as a bare per-receptor sigma*K_INTRA
   product, producing an actual K_INTRA range in physical units
   (cm^2/antibody) for the first time -- see
   k_intra_bounds_from_real_sigma() below. This is a genuine tightening
   of the bound, not a fit: it is still ~3-4 orders of magnitude wide
   because the n range (1-100) and the target-saturation range
   (0.01-1) are both still assumed, not measured, for K_INTRA
   specifically. If a future session finds a genuinely single-marker,
   surface-immobilized, CONTROLLED-copy-number EV binding curve (i.e.
   one where n is known, not ensemble-averaged), that would be the
   source that finally separates K_INTRA from n and sigma cleanly --
   this is not yet that source, but it is closer than anything found
   in the v6->v7 session.

3. Kakavoulia et al. 2026 (J Extracell Vesicles, DOI 10.1002/jev2.70269,
   PMC13088886, published this same month) -- found via a fresh literature
   search, using Flow-Induced Dispersion Analysis (FIDA) to measure
   anti-CD63 antibody binding to three EV sources. This is the RICHEST
   single-marker dataset found: real EC50 values (1.9x10^8 particles/mL
   for engineered CD63-eGFP EVs; 8.4x10^8 and 2.4x10^9 particles/mL for
   two natural ciMSC-EV preparations) AND, unusually, real per-vesicle
   receptor copy-number estimates for the SAME markers in the SAME paper
   (86+/-6 to 104+/-31 sites/EV, and 67+/-9 to 72+/-12 sites/EV, for the
   two ciMSC preps). This is the only source found this session pairing a
   real dose-response with a real copy-number estimate for one marker.
   HOWEVER: this is explicitly a SOLUTION-PHASE technique (free antibody
   + free EVs in a microfluidic channel, no surface immobilization at
   all) -- not a surface-capture measurement, which is what K_INTRA/sigma
   actually describe. The paper's own Discussion states this directly:
   "Unlike immobilization-based platforms (SPR, BLI, QCM-D), FIDA performs
   kinetic measurements entirely in-solution... Prior comparisons of FIDA
   with SPR/BLI for anti-HLA antibodies in patient serum revealed modestly
   lower affinities for FIDA, attributable to its physiological assay
   environment." Using this data to calibrate a SURFACE-capture model's
   sigma-dependent parameter would repeat, in a softer form, the exact
   mistake already caught and corrected for Mizenko's cross-marker
   co-occurrence numbers: mapping data from the wrong physical regime onto
   a parameter it doesn't actually describe. The user was shown this
   reasoning directly and chose NOT to force this calibration (see
   AskUserQuestion result this session) -- use it below only as a
   qualitative cross-check on the CONCENTRATION axis (see "What z's
   plausible range looks like" below), not as a K_INTRA fit.

WHAT THIS FILE ACTUALLY DOES: BOUNDS, NOT FITS
-------------------------------------------------
Rather than pin down K_INTRA as a single number (which would require a
unit conversion from a 3D solution Kd to a 2D surface encounter
probability that nobody has measured for EVs specifically -- the exact
kind of missing conversion factor that makes candidate #2 above
incomplete even with its real Kd number), this file bounds the
DIMENSIONLESS PRODUCT n*sigma*K_INTRA -- the actual quantity that controls
q_bound's behavior (see binding_model.py's q_bound: q_b = exp(n*sigma*K_INTRA) - 1).

The bounding argument, worked out explicitly so a future session can
check or tighten it:

1. Real EV-specific single-bond affinities are nanomolar (Kd = 1.92 nM,
   Lopez Baltazar, candidate #2 above; general antibody-antigen range
   0.1-10 nM per multiple sources found this session; RGD peptide by
   contrast is ~1000x weaker at 12-89 nM per the handoff's existing
   Anal. Chem. DOI 10.1021/acs.analchem.7b00554 citation). This is a real,
   EV-specific anchor for "how strong is one bond" -- but converting it
   into K_INTRA's actual units (1/(receptor count x ligand surface
   density)) requires a contact-patch volume and encounter-time model that
   does not exist for EVs (same gap already documented in binding_model.py
   for k_geom/patch_exponent -- see that file's module docstring).

2. A weaker but genuinely useful QUALITATIVE constraint comes from
   Mizenko's own finding, repurposed honestly (not as a numeric fit, only
   as a qualitative bound on the REGIME the model should sit in): real
   single-marker capture demonstrably misses a large fraction of
   marker-bearing EVs (this is the entire premise of Mizenko's 10x/5x
   co-occurrence finding -- single-tetraspanin capture badly undercounts
   EpCAM+/Her2+ EVs relative to non-specific capture). This means real
   assays are NOT deep in the saturated regime (theta ~ 1) for typical
   single-digit-to-low-hundreds receptor copy numbers. So n*sigma*K_INTRA
   should be O(0.01-1), not >>1, for n in the range 1-100 (the range
   documented elsewhere in this project for real EV surface markers -- see
   the handoff's literature table). This bound is honestly WIDE (four
   orders of magnitude once divided through by plausible n) -- it rules
   out "always saturated" and "always negligible" but does not pin down
   a specific operating point within that range. Presented honestly as
   such below, not dressed up as tighter than it is.

3. What CAN be bounded more usefully is the CONCENTRATION axis that z (or
   more precisely, z combined with incubation time/transport, which z is
   supposed to absorb) needs to operate over. Two independent papers,
   using two different platforms, converge on the same order of magnitude:
   Mizenko's own serial dilution operated over 4x10^7 to 4x10^9
   particles/mL (finding their practical single-particle detection ceiling
   around 2.5x10^8-5x10^8 particles/mL); Kakavoulia's EC50s (candidate #3
   above, used ONLY for this cross-check, not for a direct fit -- see the
   caveat above) sit at 1.9x10^8 to 2.4x10^9 particles/mL. This is a real,
   convergent, cross-platform finding: real EV capture/binding assays
   operate with their half-saturation point roughly at 10^8-10^9
   particles/mL. This does not directly hand us a numeric z (z is
   dimensionless and would need an assumed incubation time / diffusion
   model to convert a concentration into it), but it is a genuine,
   citable anchor for "what does a realistic operating concentration look
   like," useful for choosing DEMO figure parameters that are at least in
   the right realistic ballpark, and worth preserving as a finding in its
   own right.

   UPDATE (v7->v8 handoff transition): the Lopez Baltazar dissertation's
   retrieved full text adds a THIRD independent platform to this
   convergence, and it is now used as a genuine (not merely cross-check)
   member of the concentration-axis finding, since -- unlike Kakavoulia's
   FIDA data -- this source's antibody side really is surface-immobilized,
   the correct physical regime. The dissertation's Chapter 2 reports a
   linear dynamic range of 1.9x10^8-1.9x10^9 EVs/mL (matching the
   already-known headline number) and its Appendix A.3 control
   experiments additionally used EV concentrations of 4.9x10^9, 5.9x10^9,
   and 1.9x10^10 EVs/mL -- extending the previously-known convergent range
   upward by roughly one order of magnitude. All three sources (Mizenko,
   Kakavoulia, Lopez Baltazar) now span, together, roughly 4x10^7 to
   2x10^10 particles/mL, with every source's OWN reported "working" or
   half-saturation range still clustering within 10^8-10^10 -- the
   convergence claim is strengthened, not just repeated, by this addition.

WHAT THIS FILE DOES NOT CLAIM
----------------------------------
This is not a calibration. K_INTRA, z, and k_geom remain formally
uncalibrated in binding_model.py -- this file does not change that
module's actual default constants, and no existing test, demo, or figure
should be described as "calibrated" because of this file. What this file
adds is a documented, checkable RANGE that any future calibration attempt
should land inside of, and a record of what was searched and why it
didn't resolve to a single number.

UPDATE (v7->v8 handoff transition): the Lopez Baltazar dissertation's
full text was retrieved this session (see point 2 above) and DOES
tighten this file's bounds in one specific, honest way: K_INTRA's
plausibility range can now be expressed in actual physical units
(cm^2/antibody) using a REAL antibody surface density
(SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2 = 4.44x10^11 Abs/cm^2), rather than
left as a bare dimensionless per-receptor product with no real sigma to
divide through by. See k_intra_bounds_from_real_sigma() below. This is
still NOT a fit -- the same ~3-4-order-of-magnitude width remains,
because the n range (1-100) and target-saturation range (0.01-1) that
define the bound are still assumed physical ranges, not values read off
a real single-marker titration curve with known n. The retrieval also
did NOT resolve the deeper problem: Lopez Baltazar's Kd is a 3D bulk EV
concentration (ensemble-averaged over an uncontrolled EV population),
not a 2D areal binding constant, so it calibrates z's concentration axis
(strengthened further, see point 3 above) rather than K_INTRA itself. A
genuinely single-marker, surface-immobilized, CONTROLLED-copy-number
binding curve (known n, not an ensemble average) remains the missing
piece for an actual K_INTRA fit -- not yet found.
"""

import numpy as np


# ---------------------------------------------------------------------------
# Real numbers found this session, kept together for easy reuse/citation.
# ---------------------------------------------------------------------------

# Candidate #2: Lopez Baltazar et al. 2024, ACS Sensors, DOI 10.1021/acssensors.4c00558
# and Lopez Baltazar PhD dissertation, Cornell, DOI 10.7298/bssz-8415.
# Real, EV-specific, surface-immobilized (SPR) anti-CD81 dissociation constant.
LOPEZ_BALTAZAR_KD_M = 1.92e-9  # molar; anti-CD81, "40% surface coverage" polyclonal
LOPEZ_BALTAZAR_LINEAR_RANGE_PARTICLES_PER_ML = (1.9e8, 1.9e9)

# NEW (v7->v8 handoff transition): real antibody surface densities, retrieved
# from the dissertation's Appendix A.1 ("Antibody Surface Coverage
# Estimation," p.136) via Claude-in-Chrome this session. sigma_max is the
# maximum close-packed antibody density the SPR chip's surface chemistry
# can physically hold, derived from real antibody dimensions (not assumed);
# sigma_at_kd is the actual density at the "40% surface coverage" condition
# that produced the Kd=1.92nM headline number above -- this is the specific
# number v7 could not find ("40% coverage of an unspecified maximum").
SIGMA_MAX_LOPEZ_BALTAZAR_ABS_PER_CM2 = 1.11e12  # Abs/cm^2, 100% close-packed coverage
LOPEZ_BALTAZAR_COVERAGE_FRACTION_AT_KD = 0.40
SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2 = (
    LOPEZ_BALTAZAR_COVERAGE_FRACTION_AT_KD * SIGMA_MAX_LOPEZ_BALTAZAR_ABS_PER_CM2
)  # = 4.44e11 Abs/cm^2; also reported directly in the text as 104 ng/cm^2

# Additional real EV concentrations from the dissertation's control
# experiments (Appendix A, Figures A.2-A.3), extending the concentration-
# axis cross-check (point 3) beyond the previously-known linear range.
LOPEZ_BALTAZAR_CONTROL_CONCENTRATIONS_PARTICLES_PER_ML = (4.9e9, 5.9e9, 1.9e10)

# General antibody-antigen affinity range (multiple sources, not EV-specific)
GENERAL_ANTIBODY_KD_RANGE_M = (0.1e-9, 10e-9)

# Existing handoff citation for comparison: RGD peptide is ~1000x weaker than
# antibodies (Anal. Chem. DOI 10.1021/acs.analchem.7b00554)
RGD_KD_RANGE_M = (12.2e-9, 89e-9)  # linear RGD ~89nM, heptapeptide GRGDSPK ~12.2nM

# Candidate #3 (Kakavoulia et al. 2026, J Extracell Vesicles, DOI 10.1002/jev2.70269,
# PMC13088886) -- used ONLY as a cross-check on the concentration axis, per the
# reasoning above. NOT used to fit K_INTRA or z directly.
KAKAVOULIA_EC50_ENGINEERED_PARTICLES_PER_ML = 1.9e8   # CD63-eGFP EVs
KAKAVOULIA_EC50_NATURAL_PARTICLES_PER_ML = (8.4e8, 2.4e9)  # ciMSC EVs, two preps
KAKAVOULIA_CD63_SITES_PER_EV = (67, 104)  # approximate range across both ciMSC preps

# Candidate #1 (Mizenko et al. 2021, J Nanobiotechnology, DOI 10.1186/s12951-021-00987-1)
# own serial dilution range, re-read directly this session (Figure 2e).
MIZENKO_DILUTION_RANGE_PARTICLES_PER_ML = (4e7, 4e9)
MIZENKO_PRACTICAL_SINGLE_PARTICLE_CEILING_PARTICLES_PER_ML = (2.5e8, 5e8)


def dimensionless_exponent_bounds(n_values=(1, 10, 50, 100),
                                    target_exponents=(0.01, 0.1, 1.0)):
    """
    Compute the plausibility bound on sigma*K_INTRA (per receptor), given a
    target value for the dimensionless exponent n*sigma*K_INTRA that
    q_bound's exp(n*sigma*K_INTRA) actually depends on.

    This is NOT a fit -- it is a table showing how wide the honest bound is:
    for a given assumed receptor copy number n, what sigma*K_INTRA would be
    needed to put the model at a given point on the saturation curve
    (0.01 = far from saturation, 1.0 = order-1, entering the saturating
    regime). See the module docstring, point 2, for the qualitative
    argument (from Mizenko's own undersaturation finding) that real single-
    marker capture sits toward the LOW end of this range for typical n.

    Parameters
    ----------
    n_values : tuple of int
        Receptor copy numbers to evaluate at (default spans the range
        documented elsewhere in this project as typical for EV surface
        markers, roughly 1-100).
    target_exponents : tuple of float
        Target values of n*sigma*K_INTRA to solve for sigma*K_INTRA at.

    Returns
    -------
    dict
        Keys are (n, target_exponent) tuples, values are the implied
        sigma*K_INTRA per receptor.
    """
    results = {}
    for n in n_values:
        for target in target_exponents:
            results[(n, target)] = target / n
    return results


def concentration_axis_cross_check():
    """
    Report the real, convergent finding on the concentration axis: three
    independent papers/platforms (Mizenko's serial dilution, Kakavoulia's
    FIDA EC50s, and -- new as of the v7->v8 handoff transition -- Lopez
    Baltazar's SPR linear range and control-experiment concentrations) all
    operate with half-saturation/practical-ceiling/working concentrations
    clustering in the same order-of-magnitude band, roughly 1e8-1e10
    particles/mL overall.

    This does NOT hand us a numeric z (see module docstring point 3 for
    why), but is a genuine, citable, cross-platform convergence worth
    preserving and using to sanity-check any demo figure's assumed
    "realistic" concentration regime. Unlike Kakavoulia's data (solution-
    phase FIDA, included only as a cross-check), Lopez Baltazar's SPR
    platform is genuinely surface-immobilized -- the correct physical
    regime -- so its numbers are used here as a full member of the
    convergence claim, not just a cross-check.

    Returns
    -------
    dict
        Summary of the ranges found, for easy printing/citation.
    """
    return {
        "mizenko_dilution_range_particles_per_mL": MIZENKO_DILUTION_RANGE_PARTICLES_PER_ML,
        "mizenko_practical_ceiling_particles_per_mL": MIZENKO_PRACTICAL_SINGLE_PARTICLE_CEILING_PARTICLES_PER_ML,
        "kakavoulia_ec50_engineered_particles_per_mL": KAKAVOULIA_EC50_ENGINEERED_PARTICLES_PER_ML,
        "kakavoulia_ec50_natural_range_particles_per_mL": KAKAVOULIA_EC50_NATURAL_PARTICLES_PER_ML,
        "lopez_baltazar_linear_range_particles_per_mL": LOPEZ_BALTAZAR_LINEAR_RANGE_PARTICLES_PER_ML,
        "lopez_baltazar_control_concentrations_particles_per_mL": LOPEZ_BALTAZAR_CONTROL_CONCENTRATIONS_PARTICLES_PER_ML,
        "convergent_order_of_magnitude": "1e8-1e10 particles/mL",
        "caveat": (
            "This is a real cross-platform convergence on WHAT CONCENTRATION "
            "REGIME real EV capture/binding assays operate in, not a fit for "
            "z itself -- z is dimensionless and converting a concentration "
            "into it needs an assumed incubation time / transport model this "
            "project does not have. Kakavoulia's numbers are solution-phase "
            "(FIDA); Lopez Baltazar's are surface-immobilized (SPR, the "
            "correct regime) but ensemble-averaged over an uncontrolled EV "
            "population, so neither directly calibrates K_INTRA -- see "
            "module docstring for the full reasoning."
        ),
    }


def k_intra_bounds_from_real_sigma(n_values=(1, 10, 50, 100),
                                     target_exponents=(0.01, 0.1, 1.0),
                                     sigma_abs_per_cm2=SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2):
    """
    Convert the dimensionless-exponent plausibility bound (see
    dimensionless_exponent_bounds) into an actual K_INTRA range in physical
    units (cm^2 per antibody, per receptor), using a REAL antibody surface
    density rather than leaving the bound as a bare per-receptor
    sigma*K_INTRA product with no real sigma to divide through by.

    NEW as of the v7->v8 handoff transition, made possible by retrieving
    the Lopez Baltazar dissertation's Appendix A.1, which gives the actual
    antibody surface density (not just "40% of an unspecified maximum")
    for the first time this project has had access to.

    THIS IS STILL NOT A FIT. It answers "if the model's dimensionless
    exponent n*sigma*K_INTRA needs to be O(target) for n receptors, and
    sigma is fixed at the real value real single-marker EV-capture assays
    actually use, what would K_INTRA need to be?" -- not "here is K_INTRA."
    The n range and target-saturation range are still assumed physical
    ranges (see module docstring point 2 for the qualitative argument
    behind them), not read off a real controlled-n titration curve. The
    resulting K_INTRA range is honestly just as wide (~3-4 orders of
    magnitude) as the dimensionless bound it comes from -- dividing by a
    constant real sigma does not narrow the bound's WIDTH, only its units
    and its connection to a real, physical antibody density instead of an
    arbitrary one.

    Parameters
    ----------
    n_values : tuple of int
        Receptor copy numbers to evaluate at (same default as
        dimensionless_exponent_bounds).
    target_exponents : tuple of float
        Target values of n*sigma*K_INTRA (same default as
        dimensionless_exponent_bounds).
    sigma_abs_per_cm2 : float
        Antibody surface density to divide through by. Defaults to the
        real Lopez Baltazar value (4.44e11 Abs/cm^2, the "40% coverage"
        condition that produced Kd=1.92nM) -- pass a different value to
        check sensitivity to this choice (e.g. the 100% max packing value,
        SIGMA_MAX_LOPEZ_BALTAZAR_ABS_PER_CM2, for an upper-density bound).

    Returns
    -------
    dict
        Keys are (n, target_exponent) tuples, values are the implied
        K_INTRA in cm^2 per antibody (per receptor) -- i.e., units of
        1/(receptor count x Abs/cm^2), consistent with binding_model.py's
        q_bound(n, sigma, k_intra) where sigma is also in Abs/cm^2.
    """
    dimensionless_bounds = dimensionless_exponent_bounds(n_values, target_exponents)
    return {
        key: sigma_k_per_receptor / sigma_abs_per_cm2
        for key, sigma_k_per_receptor in dimensionless_bounds.items()
    }


def print_bounds_summary():
    """Human-readable summary of everything this file establishes."""
    print("=" * 70)
    print("K_INTRA / z / k_geom -- PLAUSIBILITY BOUNDS (not a calibration)")
    print("=" * 70)
    print()
    print(f"Real EV-specific single-bond Kd (Lopez Baltazar, anti-CD81, SPR): "
          f"{LOPEZ_BALTAZAR_KD_M*1e9:.2f} nM")
    print(f"General antibody-antigen Kd range: "
          f"{GENERAL_ANTIBODY_KD_RANGE_M[0]*1e9:.1f}-{GENERAL_ANTIBODY_KD_RANGE_M[1]*1e9:.1f} nM")
    print(f"RGD peptide Kd range (for comparison, ~1000x weaker): "
          f"{RGD_KD_RANGE_M[0]*1e9:.1f}-{RGD_KD_RANGE_M[1]*1e9:.1f} nM")
    print()
    print("Dimensionless exponent (n*sigma*K_INTRA) plausibility table:")
    bounds = dimensionless_exponent_bounds()
    for (n, target), sigma_k in sorted(bounds.items()):
        print(f"  n={n:4d}, target exponent={target:5.2f}  ->  "
              f"sigma*K_INTRA = {sigma_k:.4f} per receptor")
    print()
    print(f"NEW (v7->v8): real antibody surface density found "
          f"(Lopez Baltazar dissertation Appendix A.1):")
    print(f"  sigma_max (100% close-packed) = {SIGMA_MAX_LOPEZ_BALTAZAR_ABS_PER_CM2:.2e} Abs/cm^2")
    print(f"  sigma at Kd=1.92nM (40% coverage) = {SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2:.2e} Abs/cm^2")
    print()
    print("K_INTRA plausibility table using this REAL sigma (still not a fit -- see docstring):")
    k_intra_bounds = k_intra_bounds_from_real_sigma()
    for (n, target), k_intra in sorted(k_intra_bounds.items()):
        print(f"  n={n:4d}, target exponent={target:5.2f}  ->  "
              f"K_INTRA = {k_intra:.3e} cm^2/Ab per receptor")
    print()
    print("Concentration-axis cross-check (for z / demo realism, not a fit):")
    cc = concentration_axis_cross_check()
    print(f"  Mizenko dilution range: {cc['mizenko_dilution_range_particles_per_mL']} particles/mL")
    print(f"  Kakavoulia EC50s: engineered={cc['kakavoulia_ec50_engineered_particles_per_mL']:.1e}, "
          f"natural={cc['kakavoulia_ec50_natural_range_particles_per_mL']} particles/mL")
    print(f"  Lopez Baltazar linear range: {cc['lopez_baltazar_linear_range_particles_per_mL']} particles/mL")
    print(f"  Lopez Baltazar control concentrations: {cc['lopez_baltazar_control_concentrations_particles_per_mL']} particles/mL")
    print(f"  Convergent order of magnitude: {cc['convergent_order_of_magnitude']}")


if __name__ == "__main__":
    print_bounds_summary()
