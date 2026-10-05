"""
aim8_raw_vesicle_correlation_check.py

WHY THIS FILE EXISTS
---------------------
This file collects TWO real raw-per-vesicle-data checks, both made possible
by the user finding SVAtlas's raw-data download page
(https://www.svatlas.org/download.php) -- a page distinct from the JSON
API this project's Aim 6/7 pivot originally used, which only ever returned
already-summarized percentages. This page has ACTUAL per-vesicle data
files: one real row per individual physical vesicle, not a population
summary statistic.

PART A (MASEV tetraspanins, PANC-1/CAPAN-2) converts Aim 6's swept
correlation parameter into a directly MEASURED number from real single-
vesicle data.

PART B (CRC integrins, NEW) does the analogous thing for Aim 7: computes
real per-vesicle joint statistics for the ITGB3/ITGAM/ITGAL/ITGB2 markers
Aim 7's replication check used, directly from a real 24,000-vesicle
dataset -- rather than relying only on the published paper's already-
summarized cluster percentages the way Aim 7 originally did.

Both parts follow the exact same honesty discipline: a real raw dataset
being FROM the same paper as an earlier Aim does NOT mean it's the SAME
sample/pipeline/gating as that Aim's own summary-statistic numbers --
this file checks and reports the real numeric relationship between the two
data sources for both Parts, rather than assuming they must match.

===========================================================================
PART A: MASEV TETRASPANIN CORRELATION (unchanged from this file's first
version; see original module-level discussion for full detail)
===========================================================================

Aim 6 (aim6_masev_tetraspanin_validation.py) fits this project's Gaussian-
copula joint-reconstruction machinery to REAL marginal positivity
fractions, then SWEEPS a shared pairwise correlation parameter because no
independent measurement of that correlation existed -- only the resulting
joint category PERCENTAGES (Figure 4's summary statistics) were available.

The raw file used here (downloaded directly from
https://raw.svatlas.org/Disease/PRJ20230304.zip, project PRJ20230304,
Spitzberg et al. 2023) is the source data behind that SAME PAPER's own
Figure 6A -- described in the paper's own text as "representative maps of
12,000 single EV analyzed for 12 biomarkers... No clustering was applied
in these examples to illustrate the heterogeneity" -- confirmed by direct
inspection to contain 3000 real vesicle rows for EACH of PANC-1 and
CAPAN-2, with individual 0/1 marker-positivity columns for CD9, CD63,
CD81, and nine other real membrane markers.

IMPORTANT: this is a DIFFERENT figure than the one Aim 6 validates against
(Figure 4, N=5000, presumably gated/thresholded) -- Figure 6A is
explicitly ungated/unclustered raw display data. This module's own
verification_against_figure4_marginals() function checks and reports the
real, non-trivial difference between the two (up to 13.1 percentage
points for PANC-1's CD9 marginal) -- do not read this raw file as "proof"
of Aim 6's exact numbers, read it as an independent real measurement of
the SAME qualitative question from a related, not identical, real sample.

The paper's own precomputed "combined" columns (CD63-CD81, etc.) are
NOT used here (checked directly: they do not equal the logical AND of
the individual marker columns, and the paper's text does not explicitly
define what they represent -- see this module's original docstring
discussion). All correlation numbers below come only from the twelve
unambiguous individual-marker columns.

===========================================================================
PART B: CRC INTEGRIN CLUSTER CORRELATION (NEW)
===========================================================================

Aim 7 (aim7_crc_plasma_exosome_cluster_replication.py) uses PUBLISHED,
ALREADY-SUMMARIZED cluster percentages from Guo et al. 2023's Supplementary
Information PDF (Table S7 and Supplementary Fig. 12) to check whether a
named cancer-associated integrin cluster (ITGAM+/ITGAL+/ITGB2+) replicates
across two independent real patient cohorts. It never touched a real
individual vesicle's own marker calls directly -- everything was already
pre-clustered and pre-summarized by the paper's own analysis pipeline.

The raw file used in this Part (downloaded directly from
https://raw.svatlas.org/Disease/PRJ20230130.zip, project PRJ20230130,
the SAME Guo et al. 2023 paper Aim 7 already cites, PMID 37040507, DOI
10.34133/research.0041) is a genuinely different kind of data: a real,
long-format proximity-barcoding-assay table with 234,164 individual
marker-detection events across 24,000 REAL individual vesicles from 8
REAL patient plasma samples, covering roughly 180 real markers (not just
a hand-picked few) -- including, critically, ITGB3, ITGAM, ITGAL, and
ITGB2, the exact four markers Aim 7's replication check discusses.

WHAT THIS RAW FILE DOES AND DOES NOT LET US CHECK
------------------------------------------------------
Checked directly (not assumed) before writing any code here: the 8
patient samples in this specific downloadable file are labeled ONLY by
tumor location ("CO" = colon, "RE" = rectal) -- NOT by healthy-vs-cancer
status. The paper's own abstract (retrieved via PubMed) confirms the
FULL study design included healthy control (HC) donors and TWO
independent validation cohorts, but this particular raw file appears to
contain only the CRC-patient portion of that design (or possibly just
the discovery cohort) -- there is no way, from this file alone, to
recover a real healthy-vs-cancer comparison. This is a genuine, confirmed
limitation of what SVAtlas chose to package for download, not a bug in
this module or a gap in the search that produced it.

Given that limitation, what this file DOES let us do, honestly:
1. Compute the REAL marginal positivity fraction for each of the four
   integrin markers, directly from 24,000 real individual vesicles --
   not from the paper's own already-clustered percentages.
2. Compute the REAL pairwise (and full three-way) correlation between
   ITGAM, ITGAL, and ITGB2 -- the three markers Aim 7's replicating
   cluster11 requires ALL THREE of, directly measured rather than
   inferred from a clustering algorithm's output.
3. Compute the REAL fraction of vesicles that are simultaneously
   ITGAM+/ITGAL+/ITGB2+ (the raw-data analog of Aim 7's discovery-
   cohort cluster11), and compare it honestly to Aim 7's own published
   number for that cluster (4.58% of the whole discovery cohort) --
   again expecting a real, reportable difference rather than an exact
   match, for the same reasons as Part A (different clustering/gating
   pipeline, not the same processed output).
4. Check whether tumor location (colon vs. rectal, the one real grouping
   variable this file DOES have) shows any real difference in these
   markers -- a genuinely new, real comparison this file supports that
   the published summary tables did not report at all.
"""

import csv
import os

import numpy as np
from scipy import stats

_DATA_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..", "ev_capture_data", "raw_svatlas"
)

# ---------------------------------------------------------------------------
# PART A: MASEV tetraspanins (PANC-1 / CAPAN-2) -- unchanged from this
# file's first version.
# ---------------------------------------------------------------------------

RAW_DATA_FILES = {
    "PANC-1": os.path.join(_DATA_DIR, "panc1_fig6a_raw.csv"),
    "CAPAN-2": os.path.join(_DATA_DIR, "capan2_fig6a_raw.csv"),
}

TETRASPANIN_COLUMNS = {"CD9": "CD9", "CD63": "CD63", "CD81": "CD81sc"}

AIM6_COMMITTED_CORRELATION = 0.2


def load_raw_vesicle_data(cell_line):
    """
    Load the real per-vesicle marker-positivity table for one cell line
    from the downloaded SVAtlas raw-data CSV (converted directly from the
    original "Source data.xlsx" spreadsheet's "Figure 6a" sheet -- see
    module docstring, Part A, for provenance).

    Parameters
    ----------
    cell_line : str
        "PANC-1" or "CAPAN-2".

    Returns
    -------
    dict with:
        "CD9", "CD63", "CD81": numpy int arrays (0/1), one entry per real
            vesicle, length = real number of vesicles in the file.
        "n_vesicles": int, the real count.
    """
    if cell_line not in RAW_DATA_FILES:
        raise ValueError(
            f"cell_line must be one of {list(RAW_DATA_FILES.keys())}, got "
            f"{cell_line!r} -- this raw file only covers PANC-1 and "
            f"CAPAN-2 (see module docstring; ASPC1/A549 are not in this "
            f"particular SVAtlas raw file)."
        )

    path = RAW_DATA_FILES[cell_line]
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Raw data file not found at {path}. This file is derived from "
            f"a user-downloaded SVAtlas raw-data zip (PRJ20230304.zip from "
            f"https://raw.svatlas.org/Disease/PRJ20230304.zip) and is NOT "
            f"regenerated automatically -- see module docstring."
        )

    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    result = {"n_vesicles": len(rows)}
    for marker, col_name in TETRASPANIN_COLUMNS.items():
        result[marker] = np.array([int(row[col_name]) for row in rows], dtype=int)

    return result


def empirical_marginal_positivity(cell_line):
    """
    Compute the REAL empirical marginal positivity fraction for each
    tetraspanin, directly from the raw per-vesicle data -- a genuine
    measurement (a simple mean of a real 0/1 column), not a value derived
    from Figure 4's published percentages.

    Parameters
    ----------
    cell_line : str
        "PANC-1" or "CAPAN-2".

    Returns
    -------
    dict with "CD9", "CD63", "CD81" -> real empirical positivity fraction
    (0-1), and "n_vesicles" -> the real sample size this was computed
    from.
    """
    data = load_raw_vesicle_data(cell_line)
    return {
        "CD9": float(data["CD9"].mean()),
        "CD63": float(data["CD63"].mean()),
        "CD81": float(data["CD81"].mean()),
        "n_vesicles": data["n_vesicles"],
    }


def empirical_pairwise_correlation(cell_line):
    """
    THE MAIN MEASUREMENT PART A PROVIDES: the REAL empirical pairwise
    correlation between each pair of tetraspanins' binary positivity
    calls, computed directly from real single-vesicle data.

    Parameters
    ----------
    cell_line : str
        "PANC-1" or "CAPAN-2".

    Returns
    -------
    dict keyed by marker-pair tuple (e.g. ("CD9", "CD63")) -> dict with
    "correlation" (float, the Pearson/phi coefficient) and "p_value"
    (float, two-sided test that correlation == 0).
    """
    data = load_raw_vesicle_data(cell_line)
    markers = ("CD9", "CD63", "CD81")
    results = {}
    for i, marker_a in enumerate(markers):
        for marker_b in markers[i + 1:]:
            r, p = stats.pearsonr(data[marker_a], data[marker_b])
            results[(marker_a, marker_b)] = {"correlation": float(r), "p_value": float(p)}
    return results


def average_measured_correlation_across_pairs(cell_line):
    """
    Average the three real pairwise correlations (CD9-CD63, CD9-CD81,
    CD63-CD81) into a single number, for direct comparison against Aim
    6's SINGLE shared committed correlation (0.2).

    Parameters
    ----------
    cell_line : str
        "PANC-1" or "CAPAN-2".

    Returns
    -------
    float, the mean of the three real pairwise correlations.
    """
    pairwise = empirical_pairwise_correlation(cell_line)
    return float(np.mean([v["correlation"] for v in pairwise.values()]))


def compare_measured_correlation_to_aim6_committed(cell_line):
    """
    Real measured average pairwise correlation (from raw per-vesicle
    data) vs. Aim 6's committed correlation (0.2, swept and chosen to
    best fit Figure 4's summary percentages).

    Parameters
    ----------
    cell_line : str
        "PANC-1" or "CAPAN-2".

    Returns
    -------
    dict with "measured_average_correlation", "aim6_committed_correlation",
    "absolute_difference", "pairwise_detail".
    """
    measured = average_measured_correlation_across_pairs(cell_line)
    pairwise = empirical_pairwise_correlation(cell_line)
    return {
        "measured_average_correlation": measured,
        "aim6_committed_correlation": AIM6_COMMITTED_CORRELATION,
        "absolute_difference": abs(measured - AIM6_COMMITTED_CORRELATION),
        "pairwise_detail": pairwise,
    }


def verification_against_figure4_marginals():
    """
    Compares this raw file's real empirical marginals against the
    marginals Aim 6 derived from Figure 4's published percentages.

    Returns
    -------
    dict keyed by cell_line -> dict keyed by marker -> dict with
    "figure6a_raw_data", "figure4_summary_stats", "absolute_difference".
    """
    import sys
    sys.path.insert(0, os.path.dirname(__file__))
    from aim6_masev_tetraspanin_validation import marginals_from_real_data

    comparison = {}
    for cell_line in RAW_DATA_FILES:
        fig6a = empirical_marginal_positivity(cell_line)
        fig4 = marginals_from_real_data(cell_line)
        comparison[cell_line] = {}
        for marker in ("CD9", "CD63", "CD81"):
            comparison[cell_line][marker] = {
                "figure6a_raw_data": fig6a[marker],
                "figure4_summary_stats": fig4[marker],
                "absolute_difference": abs(fig6a[marker] - fig4[marker]),
            }
    return comparison


# ---------------------------------------------------------------------------
# PART B: CRC integrins (ITGB3/ITGAM/ITGAL/ITGB2) -- NEW.
# ---------------------------------------------------------------------------

CRC_RAW_DATA_FILE = os.path.join(_DATA_DIR, "crc", "crc_raw_vesicle_markers.csv")

CRC_MARKERS = ("ITGB3", "ITGAM", "ITGAL", "ITGB2")

# Aim 7's own published discovery-cohort cluster11 fraction (ITGAM+/ITGAL+/
# ITGB2+, whole cohort) -- hard-coded here for direct comparison, not
# imported, so this comparison is explicit (see aim7's DISCOVERY_COHORT_
# CLUSTERS["cluster11"]["fraction"] for where this number originates).
AIM7_DISCOVERY_CLUSTER11_FRACTION = 0.0458


def load_crc_raw_vesicle_data():
    """
    Load the real per-vesicle integrin-marker positivity table from the
    downloaded SVAtlas raw-data CSV (converted from the original
    "colorectal cancer.csv" long-format proximity-barcoding-assay table
    -- see module docstring, Part B, for provenance and for the honest
    limitation that this file has NO healthy-control samples, only 8 CRC
    patient samples labeled by tumor location).

    Returns
    -------
    dict with:
        "ITGB3", "ITGAM", "ITGAL", "ITGB2": numpy int arrays (0/1), one
            entry per real vesicle.
        "sample_code": numpy array of str, the real tumor-location label
            ("CO_01".."CO_05", "RE_01".."RE_03") for each vesicle -- the
            one real grouping variable this file supports (NOT healthy-
            vs-cancer).
        "n_vesicles": int, the real count (24000).
    """
    if not os.path.exists(CRC_RAW_DATA_FILE):
        raise FileNotFoundError(
            f"Raw CRC data file not found at {CRC_RAW_DATA_FILE}. This "
            f"file is derived from a user-downloaded SVAtlas raw-data zip "
            f"(PRJ20230130.zip from "
            f"https://raw.svatlas.org/Disease/PRJ20230130.zip) and is NOT "
            f"regenerated automatically -- see module docstring."
        )

    with open(CRC_RAW_DATA_FILE, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    result = {"n_vesicles": len(rows)}
    for marker in CRC_MARKERS:
        result[marker] = np.array([int(row[marker]) for row in rows], dtype=int)
    result["sample_code"] = np.array([row["sample_code"] for row in rows])
    return result


def crc_empirical_marginal_positivity():
    """
    Real empirical marginal positivity fraction for each of the four CRC
    integrin markers, directly from 24,000 real individual vesicles.

    Returns
    -------
    dict with "ITGB3", "ITGAM", "ITGAL", "ITGB2" -> real positivity
    fraction (0-1), and "n_vesicles".
    """
    data = load_crc_raw_vesicle_data()
    result = {marker: float(data[marker].mean()) for marker in CRC_MARKERS}
    result["n_vesicles"] = data["n_vesicles"]
    return result


def crc_empirical_pairwise_correlation():
    """
    Real empirical pairwise correlation between every pair of the four
    CRC integrin markers, directly measured from real single-vesicle
    data.

    Returns
    -------
    dict keyed by marker-pair tuple -> dict with "correlation", "p_value".
    """
    data = load_crc_raw_vesicle_data()
    results = {}
    for i, marker_a in enumerate(CRC_MARKERS):
        for marker_b in CRC_MARKERS[i + 1:]:
            r, p = stats.pearsonr(data[marker_a], data[marker_b])
            results[(marker_a, marker_b)] = {"correlation": float(r), "p_value": float(p)}
    return results


def crc_cluster11_like_fraction():
    """
    THE HEADLINE CRC MEASUREMENT: the REAL fraction of individual
    vesicles that are simultaneously ITGAM+/ITGAL+/ITGB2+ -- the raw-
    data analog of Aim 7's discovery-cohort cluster11 (a named cancer-
    associated, monocyte/macrophage-derived cluster), computed directly
    from real per-vesicle marker calls rather than inherited from the
    published paper's own clustering algorithm.

    Called "cluster11-LIKE" (not "cluster11") throughout this module,
    deliberately -- this is a simple triple-AND threshold on three real
    marker columns, NOT a reproduction of whatever clustering algorithm
    (e.g. FlowSOM, as seen in the paper's own supplementary figures)
    actually produced the published "cluster11" grouping. The real
    measured fraction is compared honestly to Aim 7's published number
    below, expecting a real difference (different methods), not an exact
    match.

    Returns
    -------
    dict with:
        "measured_fraction": float, real ITGAM+/ITGAL+/ITGB2+ fraction
            among all 24,000 real vesicles.
        "n_positive": int, real count of vesicles meeting this criterion.
        "n_total": int, total real vesicle count.
        "aim7_discovery_cluster11_fraction": echoed from
            AIM7_DISCOVERY_CLUSTER11_FRACTION.
        "absolute_difference": |measured - Aim 7's published value|.
    """
    data = load_crc_raw_vesicle_data()
    mask = (data["ITGAM"] == 1) & (data["ITGAL"] == 1) & (data["ITGB2"] == 1)
    measured = float(mask.mean())
    return {
        "measured_fraction": measured,
        "n_positive": int(mask.sum()),
        "n_total": data["n_vesicles"],
        "aim7_discovery_cluster11_fraction": AIM7_DISCOVERY_CLUSTER11_FRACTION,
        "absolute_difference": abs(measured - AIM7_DISCOVERY_CLUSTER11_FRACTION),
    }


def crc_marker_positivity_by_tumor_location():
    """
    A genuinely NEW real comparison this raw file supports that Aim 7's
    published summary tables never reported: does marker positivity
    differ by tumor location (colon "CO" vs. rectal "RE")? This is NOT a
    healthy-vs-cancer comparison (this file has no healthy samples -- see
    module docstring) -- it is a real, if narrower, comparison within the
    CRC-patient population this file actually contains.

    Returns
    -------
    dict keyed by marker -> dict with "CO_fraction", "RE_fraction",
    "difference" (CO - RE), and a chi-square test's p_value for
    independence between tumor location and marker positivity.
    """
    data = load_crc_raw_vesicle_data()
    is_co = np.array([code.startswith("CO") for code in data["sample_code"]])
    is_re = np.array([code.startswith("RE") for code in data["sample_code"]])

    results = {}
    for marker in CRC_MARKERS:
        co_frac = float(data[marker][is_co].mean())
        re_frac = float(data[marker][is_re].mean())

        # 2x2 contingency table: location (CO/RE) x marker (+/-)
        co_pos = int(data[marker][is_co].sum())
        co_neg = int(is_co.sum() - co_pos)
        re_pos = int(data[marker][is_re].sum())
        re_neg = int(is_re.sum() - re_pos)
        table = np.array([[co_pos, co_neg], [re_pos, re_neg]])
        chi2, p_value, _, _ = stats.chi2_contingency(table)

        results[marker] = {
            "CO_fraction": co_frac,
            "RE_fraction": re_frac,
            "difference": co_frac - re_frac,
            "p_value": float(p_value),
        }
    return results


def print_summary():
    print("=" * 78)
    print("AIM 8 -- REAL RAW-VESICLE CHECKS (PART A: MASEV tetraspanins,")
    print("PART B: CRC integrins)")
    print("=" * 78)

    print("\n" + "=" * 78)
    print("PART A: MASEV TETRASPANIN CORRELATION")
    print("=" * 78)

    print("\n--- Honest cross-check: does this raw file match Figure 4's ---")
    print("--- summary-statistic marginals Aim 6 was built from? ---")
    comparison = verification_against_figure4_marginals()
    for cell_line, markers in comparison.items():
        print(f"\n  {cell_line}:")
        print(f"  {'marker':>8} | {'Fig 6a raw':>10} | {'Fig 4 summary':>13} | {'|diff|':>8}")
        print("  " + "-" * 48)
        for marker, vals in markers.items():
            print(f"  {marker:>8} | {vals['figure6a_raw_data']*100:>9.1f}% | "
                  f"{vals['figure4_summary_stats']*100:>12.1f}% | "
                  f"{vals['absolute_difference']*100:>7.1f}pp")
    print("\n  These do NOT match exactly -- confirmed and expected (different")
    print("  figure, different N, Figure 6a is explicitly unclustered/ungated).")

    print(f"\n  HEADLINE CHECK: real measured tetraspanin correlation vs.")
    print(f"  Aim 6's committed value ({AIM6_COMMITTED_CORRELATION})")
    for cell_line in RAW_DATA_FILES:
        data = load_raw_vesicle_data(cell_line)
        print(f"\n  --- {cell_line} (real N={data['n_vesicles']} vesicles) ---")
        result = compare_measured_correlation_to_aim6_committed(cell_line)
        for (m1, m2), detail in result["pairwise_detail"].items():
            sig = "p<0.001" if detail["p_value"] < 0.001 else f"p={detail['p_value']:.4f}"
            print(f"    {m1}-{m2}: r={detail['correlation']:+.3f} ({sig})")
        print(f"    Average: {result['measured_average_correlation']:+.3f} "
              f"(Aim 6 committed: {result['aim6_committed_correlation']:+.3f}, "
              f"|diff|={result['absolute_difference']:.3f})")

    print("\n" + "=" * 78)
    print("PART B: CRC INTEGRIN CLUSTER CHECK (NEW)")
    print("=" * 78)

    marg = crc_empirical_marginal_positivity()
    print(f"\n--- Real marginal positivity (N={marg['n_vesicles']} vesicles, "
          f"8 real CRC patient samples, Guo et al. 2023) ---")
    for marker in CRC_MARKERS:
        print(f"  {marker}: {marg[marker]*100:.2f}%")

    print(f"\n--- Real pairwise correlations among the four integrin markers ---")
    pairwise = crc_empirical_pairwise_correlation()
    for (m1, m2), detail in pairwise.items():
        sig = "p<0.001" if detail["p_value"] < 0.001 else f"p={detail['p_value']:.4f}"
        print(f"  {m1}-{m2}: r={detail['correlation']:+.3f} ({sig})")

    print(f"\n--- HEADLINE CHECK: real ITGAM+/ITGAL+/ITGB2+ fraction vs. ---")
    print(f"--- Aim 7's published discovery-cohort cluster11 fraction ---")
    cluster_check = crc_cluster11_like_fraction()
    print(f"  Real measured fraction: {cluster_check['measured_fraction']*100:.3f}% "
          f"(n={cluster_check['n_positive']} of {cluster_check['n_total']} real vesicles)")
    print(f"  Aim 7's published cluster11 fraction: "
          f"{cluster_check['aim7_discovery_cluster11_fraction']*100:.2f}%")
    print(f"  Absolute difference: {cluster_check['absolute_difference']*100:.2f}pp")
    print(f"  (Expected to differ -- a simple triple-AND threshold on raw")
    print(f"   markers is NOT the same computation as the paper's own")
    print(f"   clustering algorithm; see module docstring.)")

    print(f"\n--- NEW real comparison: marker positivity by tumor location ---")
    print(f"--- (colon vs. rectal -- NOT healthy-vs-cancer, this file has ---")
    print(f"--- no healthy samples) ---")
    by_location = crc_marker_positivity_by_tumor_location()
    print(f"  {'marker':>8} | {'colon %':>8} | {'rectal %':>9} | {'diff':>7} | p-value")
    print("  " + "-" * 55)
    for marker, detail in by_location.items():
        sig = "p<0.001" if detail["p_value"] < 0.001 else f"p={detail['p_value']:.4f}"
        print(f"  {marker:>8} | {detail['CO_fraction']*100:>7.2f}% | "
              f"{detail['RE_fraction']*100:>8.2f}% | "
              f"{detail['difference']*100:>+6.2f}pp | {sig}")

    print(f"\n{'='*78}")
    print("SCOPE REMINDER (read before citing either part anywhere):")
    print("Part A and Part B are both REAL, independently measured checks from")
    print("DIFFERENT real datasets than the summary statistics Aim 6/Aim 7")
    print("originally validated against -- both from the same published papers,")
    print("but not the same underlying sample/pipeline. Report both as")
    print("'an independently measured real number, compared honestly to what")
    print("the earlier Aim's summary-statistic analysis found' -- not as proof")
    print("that the earlier Aim's exact numbers were correct. Part B's tumor-")
    print("location comparison is a genuinely NEW real finding (not previously")
    print("reported anywhere in this project or the source paper's own summary")
    print("tables) but is explicitly NOT a healthy-vs-cancer comparison --")
    print("this raw file contains no healthy control samples.")
    print("=" * 78)


if __name__ == "__main__":
    print_summary()
