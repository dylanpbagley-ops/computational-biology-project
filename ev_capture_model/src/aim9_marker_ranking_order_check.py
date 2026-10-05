"""
aim9_marker_ranking_order_check.py

WHY THIS FILE EXISTS
---------------------
This is the project's response to the diagnosis-followthrough investigation's
top-priority recommendation (docs/Diagnosis_Followthrough_and_Recommendation.md,
Part 6.1): "Build a small, new analysis that checks whether the project's
existing binding-model framework, applied to the real digitized marginal
distributions it already has for the CRC markers, correctly predicts the
RELATIVE ORDERING of diagnostic separation strength that Guo et al.'s own
paper reports."

The diagnosis document's broader point (Weak Point 1 in
docs/Novelty_Competitiveness_Diagnosis_and_Restructuring_Options.md) is that
every validation this project had done through Aim 8 REPRODUCES something
already knowable from a single existing dataset -- never a claim that
required combining the model's own machinery with real data in a genuinely
NEW way. This file is the project's first attempt at that: using the
nonlinear capture-probability model (binding_model.py) together with REAL
per-vesicle marker-positivity data (the same raw SVAtlas CRC file Aim 8,
Part B, already uses) to ask a question nobody has asked before: does this
project's own model, fed real marginal data, predict the same relative
marker ordering that an independent, real, published diagnostic-AUC
comparison (Guo et al.'s own Results section) reports?

READ THIS BEFORE TRUSTING ANY SINGLE NUMBER IN THIS FILE: WHAT DATA
ACTUALLY EXISTS, AND WHAT DOES NOT
-----------------------------------------------------------------------
Before writing any code, the following was checked directly (not assumed):

1. Guo et al.'s real, printed, healthy-vs-CRC diagnostic AUC values exist
   for exactly four markers: ITGAM (0.9158), ITGB3 (0.8681), EPS8 (0.8498),
   CD9 (0.8242). See aim7_crc_plasma_exosome_cluster_replication.py's
   GUO_ET_AL_MARKER_AUC_HEALTHY_VS_CRC for the source and exact provenance.

2. The raw SVAtlas per-vesicle CRC file (the SAME file
   aim8_raw_vesicle_correlation_check.py, Part B, uses) contains REAL
   marker-positivity columns for exactly four markers: ITGB3, ITGAM, ITGAL,
   ITGB2 -- confirmed directly by inspecting the CSV header. EPS8 and CD9
   are NOT columns in this file, and were checked (grep across this
   project's entire codebase and data directory) to not exist ANYWHERE else
   in this project's data either.

3. The raw file has NO healthy-control arm (confirmed in
   aim8_raw_vesicle_correlation_check.py's own docstring and this file's own
   re-check below) -- only 8 real CRC-patient samples, split by tumor
   location (colon/rectal), which is not a disease-status split.

The overlap of "has a real AUC" and "has real raw positivity data" is
EXACTLY ONE MARKER: ITGAM. ITGB3 has a real AUC but no positivity data in
this file (only ITGAM/ITGAL/ITGB2 co-occur in a named cluster the paper
discusses; ITGB3 appears in a DIFFERENT discovery-cohort cluster -- see
aim7's DISCOVERY_COHORT_CLUSTERS). EPS8 and CD9 have real AUCs but zero real
positivity data anywhere in this project.

THIS MEANS: a strict 4-marker Spearman-style ranking-correlation check
against Guo et al.'s AUC ranking is NOT POSSIBLE with real data currently in
hand -- attempting one would require fabricating EPS8/CD9 marginal
distributions from nothing, which this project's standing discipline
(never invent precision data doesn't support) forbids. This file does NOT
do that. Instead, it runs two clearly separated, honestly-scoped checks:

CHECK 1 (PRIMARY, fully real-data-only, no fabrication): does the
nonlinear model REORDER the four markers it CAN compute from real data
(ITGB3, ITGAM, ITGAL, ITGB2), relative to their raw/linear positivity-
fraction ranking, anywhere across the project's standing K_INTRA/z
plausibility sweep? This is an internal-consistency / sensitivity finding
-- it needs no ground truth and fabricates nothing. It directly tests
this project's central thesis (nonlinear capture can reorder markers a
linear comparison would rank differently) on a fourth real dataset,
independent of the three systems aim3_aim4_range_sweep.py and
stem_cell_marker_ranking_extension.py already swept.

CHECK 2 (SECONDARY): an attempt to build a single-marker (ITGAM-only)
synthetic healthy-vs-CRC comparison, using the one real, paper-reported
healthy-to-CRC decline ratio this project has anywhere near this marker
family (the ITGAM+/ITGAL+/ITGB2+ TRIPLE-POSITIVE cluster's real 13.47%
healthy -> 1.13% CRC decline). CAUGHT DIRECTLY WHILE FIRST RUNNING THIS
FILE: applying that cluster-level ratio to ITGAM's own individual real
positivity fraction (14.62%) produces 174.3% -- a mathematically
IMPOSSIBLE positivity fraction. This is not a weak assumption; it is a
structurally INVALID one (a joint/intersection-population ratio does not
transfer to a single marker's marginal), and the impossible arithmetic
proves that directly. Check 2 therefore does NOT attempt to force a
synthetic healthy population through anyway (no clipping, no arbitrary
ceiling) -- it reports the impossibility itself as the honest finding. See
itgam_cluster_ratio_transfer_check() below for the full account, kept in
the code (not silently deleted) precisely because a caught, understood,
correctly-reported dead end is itself real information for a future
session (matching this project's standing "report a well-scoped attempt's
surprising or negative result plainly" discipline).

WHAT THIS FILE DOES NOT CLAIM
--------------------------------
- It does NOT claim to have validated Guo et al.'s 4-marker AUC ranking.
  That would require real per-marker healthy AND cancer marginal data for
  all four markers, which does not exist in this project (or, per the
  diagnosis-followthrough investigation's own honest search, appears not
  to exist in any publicly accessible form at all right now).
- It does NOT invent EPS8 or CD9 marginal distributions. These two markers
  are excluded from every quantitative computation in this file.
- It does NOT construct a synthetic ITGAM healthy population by any
  invalid ratio-transfer -- Check 2 reports why that specific approach is
  invalid rather than forcing a number through it (see above).

INDEPENDENT SUPPORTING LITERATURE FOR THE BROADER "MARKER CHOICE/
COMBINATION CHANGES DIAGNOSTIC POWER" THESIS (added third and fourth
sessions, August 30, 2026 -- citation-only, no computation in this file
uses any of these numbers)
-----------------------------------------------------------------------
This file's own Check 1/Check 2 are real but narrow (four integrin markers,
one CRC data source). The broader claim this project makes -- that marker
CHOICE and marker COMBINATION materially change diagnostic power, which is
exactly what a nonlinear (density-dependent) capture model would predict
but a naive fold-change comparison could miss -- has substantial INDEPENDENT
support from real, unrelated data sources found this project's third and
fourth research sessions:

- Shi, Zhao, Xue & Jia, "Extracellular vesicle biomarkers in circulation
  for colorectal cancer detection: a systematic review and meta-analysis,"
  BMC Cancer 24:623 (2024), DOI 10.1186/s12885-024-12312-8, PMID 38778252,
  PMC11110411. A real meta-analysis of 48 studies, 7,727 participants, 162
  individual EV biomarkers. Key real finding: EV PROTEIN PANELS reach AUC
  0.92 (85% sensitivity / 84% specificity) vs. 0.87-0.90 for individual
  markers or RNA-based markers -- independent, large-scale, aggregated
  confirmation that marker COMBINATION beats any single marker, consistent
  with this project's own Aim 5 AND-gate finding, from a source with no
  connection to this project's own CRC data (Guo et al.).

- The CRC "supermeres" paper (a brand-new, Jan 2026 finding): [title
  redacted for brevity in this docstring -- see below], DOI
  10.1038/s41598-026-36626-8, PMID 41571721, PMC12902025. Reports real
  head-to-head diagnostic AUC comparisons across THREE different
  nanoparticle types (supermeres, exomeres, small EVs) for the SAME five
  shared markers (TGFBI, ENO1, DPEP1, GPC1, CEA), in a real 22-CRC-vs-25-
  healthy-control cohort -- independent, very recent evidence that both
  marker choice AND particle-type choice materially change diagnostic AUC.
  Not usable as per-vesicle model input (bulk/SPR measurement method, a
  different nanoparticle class than this project models), but a real,
  citable, independent data point for the same broader thesis.

- The PDAC (pancreatic) integrin paper: Scientific Reports 2024, DOI
  10.1038/s41598-024-65209-8, PMID 38902362, PMC11189911 -- "Differential
  extracellular vesicle concentration and their biomarker expression of
  integrin alpha/beta, EpCAM, and glypican-1 in pancreatic cancer models."
  Real per-marker, per-cell-line ranking differences (e.g., GPC-1 showed
  only a "mild increase" in PANC-1 cells specifically, but was significant
  in BxPC-3 cells) via Western blot relative expression -- a DIFFERENT
  measurement method (bulk/relative, not per-vesicle) than this project's
  own approach, so not usable as quantitative model input, but real,
  citable, independent evidence that marker ranking is not universal across
  even closely related cell lines of the same cancer type.

- TWO BRAND-NEW (2026) papers found via an explicitly open-ended search
  pass (not itemized against either prior diagnosis document's brainstorm
  list), both reporting real multi-marker EV panels with strong combined
  diagnostic performance and BOTH explicitly naming single-EV/single-
  vesicle resolution as an open gap in their own work:
    (a) Lobb, Zhou, Fielding et al., "Early-stage multi-cancer detection
        through a plasma extracellular vesicle protein signature," Cell
        Reports Medicine (2026), DOI 10.1016/j.xcrm.2026.102694, PMID
        41881025, PMC13130687. A 4-protein small-EV panel (thrombospondin-1,
        nidogen-1, pentraxin-3, versican), discovered via an ISOGENIC cell
        line pair (matched human bronchial epithelial cells, differing only
        by a deliberate p53-knockdown/KRAS-overexpression transformation --
        a clean, confounder-free healthy-vs-cancer comparison at the
        cell-line level), validated across 22 cancer cell lines and 764
        clinical plasma/serum samples (AUC 0.91-1.00), including colorectal
        cancer as one of the 8 clinical cohort cancer types. Its own
        Limitations section states directly that its SERS-based readout
        "measures bulk sEV populations enriched by immunocapture, without
        resolving single-vesicle heterogeneity" and that "future single-EV
        studies will be important" -- a real, recent, independent paper
        explicitly naming the exact gap this project's single-vesicle
        modeling approach addresses. Proteomics data (the isogenic
        discovery phase only) deposited at ProteomeXchange/PRIDE, accession
        PXD073635 -- summary/relative-intensity data, not per-vesicle
        counts, so not usable as this project's model input, but real and
        citable.
    (b) Tamai, Chiwaki, Shiotani et al., "EV-Finder: Direct Detection of
        Extracellular Vesicle-Associated Proteins by Proximity Extension
        Assay for Multi-Cancer Screening," International Journal of
        Molecular Sciences (2026), DOI 10.3390/ijms27114904, PMID 42278434,
        PMC13256567. Real serum samples across five cancer types (n=193)
        vs. independent healthy controls (n=138), AUC 0.985, sensitivity
        0.929, specificity 0.957 -- another real, very recent (2026),
        independent multi-marker EV diagnostic result, further supporting
        that this remains an active, fast-moving research area with new
        marker panels and new performance numbers appearing on a timescale
        of months, not years.

None of these five sources changes any number, function, or computation in
this file -- Check 1 and Check 2 above remain exactly as described, using
only the real Guo et al. and raw-SVAtlas-CRC data already documented above.
These are citation-only additions supporting the broader thesis this file's
narrow, honestly-scoped checks are one small, real piece of.

THE AUC-LIKE SEPARATION METRIC USED, AND WHY
-------------------------------------------------
No AUC/ROC/discriminability metric of any kind exists anywhere else in this
codebase (checked directly -- every existing "ranking" function in
stem_cell_marker_ranking_extension.py and aim3_aim4_range_sweep.py ranks
markers by scalar YIELD magnitude, not by how well two populations'
capture-probability distributions SEPARATE). This file adds one, built from
scratch: the rank-biserial / Mann-Whitney U definition of AUC,
    AUC = P(theta_cancer_vesicle > theta_healthy_vesicle),
estimated directly via scipy.stats.mannwhitneyu's U statistic
(AUC = U / (n1 * n2)) applied to two samples of per-vesicle capture
probabilities theta = capture_probability(n, sigma, k_intra, z), each drawn
from a negative-binomial receptor-count distribution built the same way
stem_cell_marker_ranking_extension.py's fit_negative_binomial_to_positivity()
already does (project standing convention -- negative binomial, never
Poisson/Gaussian, for real-data reconstruction). This is a genuinely NEW
building block for this codebase, reused nowhere else yet. It is NOT
actually usable in this file's Check 2 (see above -- no valid two-
population comparison exists to run it on for ITGAM), but is kept as a
reusable function for the reordering analysis (see Check 1) and any future
session that DOES find a valid two-population real dataset for this
project.
"""

import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from binding_model import capture_probability
from k_intra_z_plausibility_bounds import SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2
from stem_cell_marker_ranking_extension import fit_negative_binomial_to_positivity
from aim8_raw_vesicle_correlation_check import (
    crc_empirical_marginal_positivity,
    CRC_MARKERS,
)
from aim7_crc_plasma_exosome_cluster_replication import (
    GUO_ET_AL_MARKER_AUC_HEALTHY_VS_CRC,
    GUO_ET_AL_CLUSTER11_HEALTHY_TO_CANCER_DECLINE,
)

# ---------------------------------------------------------------------------
# Real data this file uses, gathered from the two sources above -- kept as
# module-level constants here (not just imported and used inline) so a
# reader can see, in one place, exactly what real numbers this whole
# analysis rests on.
# ---------------------------------------------------------------------------

# The four markers with real raw CRC positivity data (ITGB3, ITGAM, ITGAL,
# ITGB2 -- see CRC_MARKERS in aim8). EPS8 and CD9 are deliberately absent --
# no real positivity data exists for them anywhere in this project.
MARKERS_WITH_REAL_CRC_POSITIVITY = CRC_MARKERS

# The one marker with BOTH a real Guo AUC AND real CRC positivity data.
MARKER_WITH_BOTH_REAL_AUC_AND_POSITIVITY = "ITGAM"


def _reconstruct_crc_population(marker, positivity_fraction, n_samples=5000, seed=0):
    """
    Build a synthetic per-vesicle receptor-count sample for one CRC marker
    from its REAL positivity fraction, using this project's standing
    negative-binomial reconstruction convention
    (fit_negative_binomial_to_positivity, from
    stem_cell_marker_ranking_extension.py -- reused directly here, not
    reimplemented).

    IMPORTANT LIMITATION, STATED HONESTLY: this raw CRC file gives 0/1
    presence/absence per marker, not an actual copy-number count -- there
    is no "mean copies per positive particle" number available for these
    four markers (unlike the stem-cell system, which has both numbers
    reported). Without that second number, fit_negative_binomial_to_positivity
    needs a stand-in for mean_copies_if_positive. This function uses a
    FIXED, project-wide generic value (see GENERIC_MEAN_COPIES_IF_POSITIVE
    below) for every marker, so this assumption affects all four markers
    EQUALLY and cannot, by construction, be the reason any ranking
    differences appear between them -- exactly the same "hold everything
    but the one real, differing number constant" discipline
    stem_cell_marker_ranking_extension.py uses for its own diameter
    assumption (see that file's module docstring, "Diameter" section, for
    the precedent).

    Parameters
    ----------
    marker : str
        One of MARKERS_WITH_REAL_CRC_POSITIVITY.
    positivity_fraction : float
        Real fraction of vesicles positive for this marker (must be
        strictly between 0 and 1 -- see
        fit_negative_binomial_to_positivity's own validation).
    n_samples : int
        Synthetic sample size (5000, matching this project's standing
        convention for reconstructed populations -- see
        stem_cell_marker_ranking_extension.reconstruct_marker_population's
        own choice of 5000 for the same reason).
    seed : int
        Fixed random seed for exact reproducibility.

    Returns
    -------
    numpy array of length n_samples, real-valued (non-negative) synthetic
    receptor counts.
    """
    fit = fit_negative_binomial_to_positivity(
        positivity_fraction, GENERIC_MEAN_COPIES_IF_POSITIVE
    )
    rng = np.random.default_rng(seed)
    p_param = fit["dispersion_r"] / (fit["dispersion_r"] + fit["mean_n"])
    n_values = rng.negative_binomial(fit["dispersion_r"], p_param, size=n_samples)
    return n_values.astype(float)


# Generic, literature-typical mean-copies-if-positive stand-in -- NOT
# marker-specific, NOT measured for these four integrins. Needed because
# the raw CRC file only gives 0/1 positivity, not copy number, and
# fit_negative_binomial_to_positivity requires both a zero-fraction target
# AND a positive-conditional-mean target. Applying the SAME value to every
# marker (see _reconstruct_crc_population's docstring) means this
# assumption cannot itself manufacture a ranking difference between
# markers -- any reordering found below comes from the REAL, differing
# positivity fractions, not from this shared assumption. Chosen as a
# round, moderate value broadly consistent with the copy-number scales
# already used elsewhere in this project for other real integrin/
# tetraspanin markers (e.g. stem_cell_marker_ranking_extension.py's real
# CD81/CD63/CD73 copies-per-positive-particle numbers, all in the 20-53
# range) -- flagged explicitly as an assumption, not a measurement.
GENERIC_MEAN_COPIES_IF_POSITIVE = 30.0


def auc_from_two_samples(sample_a, sample_b):
    """
    A rank-biserial / Mann-Whitney U definition of AUC -- AUC = P(a
    randomly drawn value from sample_a exceeds a randomly drawn value from
    sample_b). This did not exist anywhere in this codebase before this
    file (checked directly -- see module docstring). Provided as a
    reusable building block; see module docstring for why it is not
    actually usable in this file's Check 2 for ITGAM specifically (no
    valid two-population comparison could be constructed -- see
    itgam_cluster_ratio_transfer_check() below).

    Parameters
    ----------
    sample_a : array
        Values for the "positive" group. By convention, AUC > 0.5 means
        sample_a tends to have HIGHER values than sample_b.
    sample_b : array
        Values for the "negative"/reference group.

    Returns
    -------
    dict with:
        "auc": float, between 0 and 1. 0.5 = no separation (random);
            1.0 = perfect separation (every sample_a value exceeds every
            sample_b value); 0.0 = perfect separation in the OPPOSITE
            direction.
        "u_statistic": the raw Mann-Whitney U statistic.
        "p_value": two-sided p-value that the two samples come from
            distributions with the same location (from
            scipy.stats.mannwhitneyu).
    """
    sample_a = np.asarray(sample_a, dtype=float)
    sample_b = np.asarray(sample_b, dtype=float)
    u_stat, p_value = stats.mannwhitneyu(sample_a, sample_b, alternative="two-sided")
    auc = float(u_stat) / (len(sample_a) * len(sample_b))
    return {"auc": auc, "u_statistic": float(u_stat), "p_value": float(p_value)}


# ---------------------------------------------------------------------------
# CHECK 1 (PRIMARY): does the nonlinear model reorder the four real-data
# markers relative to their raw/linear positivity ranking, anywhere across
# the K_INTRA/z plausibility sweep?
# ---------------------------------------------------------------------------

def linear_ranking_from_real_positivity():
    """
    The STANDARD FIELD PRACTICE this project critiques: rank the four real
    CRC markers by raw positivity fraction alone (the "naive fold-change"
    analog used elsewhere in this project, e.g.
    stem_cell_marker_ranking_extension.naive_fold_change_ranking) -- no
    model, no nonlinearity, just the real measured fraction itself.

    Returns
    -------
    dict with "positivity_fractions" (real, from
    crc_empirical_marginal_positivity()) and "ranking_best_first" (marker
    names sorted by positivity fraction, highest first).
    """
    marg = crc_empirical_marginal_positivity()
    fractions = {m: marg[m] for m in MARKERS_WITH_REAL_CRC_POSITIVITY}
    ranking = sorted(fractions, key=fractions.get, reverse=True)
    return {"positivity_fractions": fractions, "ranking_best_first": ranking}


def nonlinear_yield_ranking(sigma, k_intra, z, n_samples=5000, seed=0):
    """
    Rank the same four real CRC markers by the FULL NONLINEAR model's mean
    capture probability (yield), at one (sigma, K_INTRA, z) point, using
    receptor-count populations reconstructed from each marker's REAL
    positivity fraction (see _reconstruct_crc_population).

    Parameters
    ----------
    sigma : float
        Ligand surface density.
    k_intra : float
        Bond-strength parameter.
    z : float
        Activity/prefactor.
    n_samples, seed :
        Passed through to _reconstruct_crc_population for every marker.

    Returns
    -------
    dict with "yields" (marker -> mean capture probability) and
    "ranking_best_first".
    """
    marg = crc_empirical_marginal_positivity()
    yields = {}
    for marker in MARKERS_WITH_REAL_CRC_POSITIVITY:
        n_values = _reconstruct_crc_population(
            marker, marg[marker], n_samples=n_samples, seed=seed
        )
        theta = capture_probability(n_values, sigma, k_intra, z)
        yields[marker] = float(np.mean(theta))
    ranking = sorted(yields, key=yields.get, reverse=True)
    return {"yields": yields, "ranking_best_first": ranking}


def sweep_reordering_frequency(sigma=SIGMA_LOPEZ_BALTAZAR_ABS_PER_CM2,
                                 n_samples=5000, seed=0):
    """
    THE MAIN CHECK 1 RESULT: across the project's FULL standing K_INTRA/z
    plausibility grid (the same grid aim3_aim4_range_sweep.py and
    stem_cell_marker_ranking_extension.py already sweep -- imported here,
    not redefined, so this file always tracks the project's one shared
    plausibility range rather than drifting out of sync with it), what
    fraction of grid points produce a nonlinear-model ranking that DIFFERS
    from the real, linear positivity-fraction ranking?

    This is a real, fully-defensible result that needs no ground truth and
    fabricates no data -- it directly tests this project's central thesis
    (nonlinear capture can reorder markers a linear comparison would rank
    differently) using ONLY the four markers' REAL positivity fractions.

    Returns
    -------
    dict with:
        "linear_ranking": the real, data-only ranking (see
            linear_ranking_from_real_positivity).
        "grid_results": list of dicts, one per (k_intra, z) grid point,
            each with "k_intra", "z", "nonlinear_ranking", "matches_linear".
        "n_grid_points": total grid points swept.
        "n_reordered": how many grid points produced a DIFFERENT ranking
            than the linear one.
        "reordering_fraction": n_reordered / n_grid_points.
        "itgam_rank_positions": list of the rank position (1=best) ITGAM
            held at each grid point -- used by print_summary() to report
            whether ITGAM (the one marker with a real external AUC) stays
            in a stable, high-ranked position across the sweep, as a
            qualitative cross-check (NOT a validation -- see module
            docstring).
    """
    # Imported here (not at module level) to avoid a module-load-order
    # dependency between aim9 and aim3_aim4 -- matches this project's
    # existing pattern of importing get_k_intra_sweep_values/
    # get_z_sweep_values lazily where needed (see
    # aim4_distribution_and_assumption_sensitivity.py for the same
    # pattern).
    from aim3_aim4_range_sweep import get_k_intra_sweep_values, get_z_sweep_values

    linear = linear_ranking_from_real_positivity()
    k_intra_values = get_k_intra_sweep_values()
    z_values = get_z_sweep_values()

    grid_results = []
    itgam_rank_positions = []
    n_reordered = 0
    for k_intra in k_intra_values:
        for z in z_values:
            nonlinear = nonlinear_yield_ranking(
                sigma, k_intra, z, n_samples=n_samples, seed=seed
            )
            matches = nonlinear["ranking_best_first"] == linear["ranking_best_first"]
            if not matches:
                n_reordered += 1
            grid_results.append({
                "k_intra": float(k_intra),
                "z": float(z),
                "nonlinear_ranking": nonlinear["ranking_best_first"],
                "matches_linear": matches,
            })
            itgam_rank_positions.append(
                nonlinear["ranking_best_first"].index("ITGAM") + 1
            )

    n_grid_points = len(grid_results)
    return {
        "linear_ranking": linear,
        "grid_results": grid_results,
        "n_grid_points": n_grid_points,
        "n_reordered": n_reordered,
        "reordering_fraction": n_reordered / n_grid_points if n_grid_points else float("nan"),
        "itgam_rank_positions": itgam_rank_positions,
    }


# ---------------------------------------------------------------------------
# CHECK 2 (SECONDARY): tests whether a naive "borrow the cluster's
# healthy:CRC decline ratio for one of its own member markers" approach is
# even valid -- and finds, honestly, that it is NOT.
# ---------------------------------------------------------------------------

def itgam_cluster_ratio_transfer_check():
    """
    WHAT THIS FUNCTION ORIGINALLY TRIED TO DO, AND WHY THAT APPROACH WAS
    ABANDONED (kept here, explained, rather than silently deleted, because
    this is itself a real, useful, honestly-reported finding): the first
    version of this file tried to build a synthetic "healthy" positivity
    fraction for ITGAM by taking ITGAM's REAL individual CRC positivity
    fraction (14.62%, from the raw SVAtlas file) and multiplying it by the
    ONE real healthy:CRC decline ratio this project has for this marker
    family -- the ITGAM+/ITGAL+/ITGB2+ TRIPLE-POSITIVE cluster's real,
    paper-reported 13.47% (healthy) -> 1.13% (CRC) decline (a ratio of
    ~11.9x).

    CAUGHT DIRECTLY WHILE FIRST RUNNING THIS FILE (not before, and not
    glossed over -- per this project's standing "report a surprising or
    negative result plainly" discipline): 14.62% * 11.9x = 174.3%, a
    mathematically IMPOSSIBLE positivity fraction. This is not just a weak
    assumption -- it is a STRUCTURALLY INVALID one, and the impossible
    result proves it directly rather than requiring a separate argument.
    The reason: the 13.47%/1.13% numbers describe a TRIPLE-POSITIVE
    (ITGAM+ AND ITGAL+ AND ITGB2+ simultaneously) population, which is
    necessarily much RARER than any single marker's own individual
    positivity -- 14.62% for ITGAM ALONE is already larger than 13.47% for
    the THREE-MARKER INTERSECTION in healthy controls, before any ratio is
    even applied. A ratio computed on a rare joint/intersection population
    does not validly transfer to a single marker's own marginal
    probability -- these are different statistical objects, and no
    honest labeling or capping of the resulting number fixes that; it
    would only hide the same invalid computation behind a differently
    arbitrary ceiling.

    WHAT THIS FUNCTION DOES INSTEAD: reports the arithmetic that PROVES the
    invalidity, as a genuine (if negative) finding, rather than either (a)
    silently clipping the impossible value to something plausible-looking
    (what the abandoned first version did, until this was caught), or (b)
    quietly deleting the attempt with no record of why it was tried and
    rejected. No AUC-like separation metric is computed here, because there
    is no valid synthetic healthy population to compute one from. A future
    session wanting to attempt a per-marker healthy estimate differently
    would need a genuinely different anchor -- e.g., a binormal-ROC
    decomposition of the AUC itself, under an explicit, checkable
    equal-variance assumption -- not attempted here because this project's
    dataset does not carry the fields (a stated decision threshold, a
    comparable-cohort confirmation) that assumption would need to be
    checked rather than just asserted.

    Returns
    -------
    dict with "itgam_real_crc_fraction", "cluster_healthy_fraction",
    "cluster_crc_fraction", "naive_decline_ratio",
    "naive_transferred_fraction" (the impossible >1.0 result), and
    "is_valid" (False by construction, for the real numbers involved --
    computed directly rather than hard-coded, so this would correctly flip
    to True if the underlying real numbers ever changed).
    """
    marg = crc_empirical_marginal_positivity()
    itgam_real_crc_fraction = float(marg["ITGAM"])

    decline = GUO_ET_AL_CLUSTER11_HEALTHY_TO_CANCER_DECLINE
    cluster_healthy_fraction = decline["healthy_fraction"]
    cluster_crc_fraction = decline["primary_and_metastatic_crc_fraction"]
    naive_decline_ratio = cluster_healthy_fraction / cluster_crc_fraction

    naive_transferred_fraction = itgam_real_crc_fraction * naive_decline_ratio

    return {
        "itgam_real_crc_fraction": itgam_real_crc_fraction,
        "cluster_healthy_fraction": cluster_healthy_fraction,
        "cluster_crc_fraction": cluster_crc_fraction,
        "naive_decline_ratio": naive_decline_ratio,
        "naive_transferred_fraction": naive_transferred_fraction,
        "is_valid": naive_transferred_fraction <= 1.0,
    }


def print_summary():
    print("=" * 78)
    print("AIM 9 -- CRC MARKER RANKING-ORDER CHECK")
    print("(the project's first new-claim check: does the nonlinear model,")
    print(" fed real marginal data, predict a real, independently-published")
    print(" AUC ranking? See module docstring for exactly what real data")
    print(" does and does not exist to check this.)")
    print("=" * 78)

    print("\n" + "=" * 78)
    print("REAL GROUND TRUTH (Guo et al.'s own reported healthy-vs-CRC AUCs)")
    print("=" * 78)
    for marker, auc in sorted(GUO_ET_AL_MARKER_AUC_HEALTHY_VS_CRC.items(),
                                key=lambda kv: kv[1], reverse=True):
        has_data = ("  <- has real raw CRC positivity data"
                    if marker in MARKERS_WITH_REAL_CRC_POSITIVITY
                    else "  <- NO raw data anywhere in this project")
        print(f"  {marker:>6}: AUC = {auc:.4f}{has_data}")
    print(f"\n  Overlap of (real AUC) and (real raw positivity data): "
          f"{MARKER_WITH_BOTH_REAL_AUC_AND_POSITIVITY} only.")
    print("  A strict 4-marker ranking-correlation check is NOT POSSIBLE with")
    print("  real data currently in hand -- see module docstring for why. This")
    print("  file runs two honestly-scoped checks instead (below).")

    print("\n" + "=" * 78)
    print("CHECK 1 (PRIMARY): does the nonlinear model reorder the four")
    print("real-data markers (ITGB3/ITGAM/ITGAL/ITGB2), relative to their")
    print("real, linear positivity-fraction ranking, anywhere across the")
    print("project's full K_INTRA/z plausibility sweep?")
    print("=" * 78)

    result = sweep_reordering_frequency()
    linear = result["linear_ranking"]
    print(f"\n  Real linear (positivity-fraction) ranking, best first: "
          f"{' > '.join(linear['ranking_best_first'])}")
    for marker in linear["ranking_best_first"]:
        print(f"    {marker:>6}: {linear['positivity_fractions'][marker]*100:.2f}% "
              f"(real, raw SVAtlas CRC data)")

    print(f"\n  Swept {result['n_grid_points']} (K_INTRA, z) grid points "
          f"(the project's standing plausibility grid).")
    print(f"  Reordered relative to the linear ranking at "
          f"{result['n_reordered']} of {result['n_grid_points']} points "
          f"({result['reordering_fraction']*100:.1f}%).")

    itgam_positions = result["itgam_rank_positions"]
    print(f"\n  ITGAM's rank position across the sweep (1=best of 4): "
          f"min={min(itgam_positions)}, max={max(itgam_positions)}, "
          f"most common={max(set(itgam_positions), key=itgam_positions.count)}")
    print("  (Qualitative cross-check only, NOT a validation -- ITGAM is the")
    print("   one marker with a real external AUC, and it is also the highest")
    print("   of the four real Guo AUCs. Whether ITGAM stays highly ranked")
    print("   across this sweep is worth noting, but n=1 marker overlap means")
    print("   this cannot be scored as a ranking match/mismatch.)")

    if result["reordering_fraction"] > 0:
        print(f"\n  FINDING: the nonlinear model DOES reorder these four real")
        print(f"  markers relative to naive linear positivity, at "
              f"{result['reordering_fraction']*100:.1f}% of the plausibility")
        print(f"  grid -- a real, data-grounded demonstration of this project's")
        print(f"  central thesis on a fourth independent real dataset.")
    else:
        print(f"\n  FINDING: across this entire plausibility grid, the nonlinear")
        print(f"  model did NOT reorder these four markers relative to naive")
        print(f"  linear positivity -- an honest negative finding for this")
        print(f"  specific dataset (report this directly, do not reach for a")
        print(f"  different grid/dataset to force a positive finding).")

    print("\n" + "=" * 78)
    print("CHECK 2 (SECONDARY): is a naive cluster-ratio-transfer to a")
    print("single member marker (ITGAM) even valid? Tested directly.")
    print("=" * 78)

    check2 = itgam_cluster_ratio_transfer_check()
    print(f"\n  ITGAM real CRC positivity (raw SVAtlas data): "
          f"{check2['itgam_real_crc_fraction']*100:.2f}%")
    print(f"  ITGAM+/ITGAL+/ITGB2+ cluster, healthy controls (Guo et al.): "
          f"{check2['cluster_healthy_fraction']*100:.2f}%")
    print(f"  ITGAM+/ITGAL+/ITGB2+ cluster, primary+metastatic CRC (Guo et al.): "
          f"{check2['cluster_crc_fraction']*100:.2f}%")
    print(f"  Naive cluster healthy:CRC decline ratio: "
          f"{check2['naive_decline_ratio']:.2f}x")
    print(f"\n  Naively transferring this ratio to ITGAM's own real CRC fraction:")
    print(f"    {check2['itgam_real_crc_fraction']*100:.2f}% x "
          f"{check2['naive_decline_ratio']:.2f} = "
          f"{check2['naive_transferred_fraction']*100:.1f}% "
          f"(IMPOSSIBLE -- a fraction cannot exceed 100%)")
    print(f"\n  FINDING: this specific ratio-transfer is NOT just weak, it is")
    print(f"  mathematically INVALID -- the cluster ratio describes a rare")
    print(f"  triple-marker INTERSECTION, which cannot be transferred to a")
    print(f"  single marker's own, much more common, individual marginal.")
    print(f"  is_valid = {check2['is_valid']}. No synthetic healthy-ITGAM")
    print(f"  population is constructed from this approach (no clipping, no")
    print(f"  arbitrary ceiling) -- this negative finding is reported plainly")
    print(f"  instead. See this function's own docstring for what a valid")
    print(f"  alternative approach would need (a checkable binormal-ROC")
    print(f"  decomposition of the AUC), not attempted this session.")

    print(f"\n{'='*78}")
    print("SCOPE REMINDER (read before citing this file anywhere):")
    print("Check 1 is real, data-only, and fabricates nothing -- cite it as a")
    print("sensitivity/reordering finding, the same class of evidence as")
    print("aim3_aim4_range_sweep.py's ranking-flip checks, on a NEW dataset.")
    print("Check 2 found a real, informative dead end (a natural-seeming")
    print("ratio-transfer that turns out to be mathematically invalid) -- cite")
    print("it only as that: a caught methodological pitfall, never as a")
    print("healthy-ITGAM estimate or a step toward validating Guo et al.'s AUC")
    print("ranking.")
    print("=" * 78)


if __name__ == "__main__":
    print_summary()
