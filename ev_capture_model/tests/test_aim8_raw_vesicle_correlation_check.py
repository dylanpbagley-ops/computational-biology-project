"""
test_aim8_raw_vesicle_correlation_check.py

Regression tests for src/aim8_raw_vesicle_correlation_check.py -- the
project's first check using genuinely raw (not summary-statistic) real
per-vesicle data. See that file's own module docstring for the full
account of what this validates, how it differs from Aim 6, and why this
raw file is NOT the same dataset as Figure 4's summary percentages.

WHAT THESE TESTS LOCK IN, AND WHY EACH ONE MATTERS
-------------------------------------------------------
1. The raw data file loads correctly and has the expected real sample
   size (3000 vesicles per cell line) -- a basic sanity check that the
   CSV conversion from the original spreadsheet didn't drop or duplicate
   rows.
2. load_raw_vesicle_data() rejects an unknown cell line with a clear
   error (this raw file only covers PANC-1/CAPAN-2, not ASPC1/A549).
3. empirical_marginal_positivity() returns values in [0, 1] and matches
   an independent hand-recomputation from the raw arrays directly (not
   just re-deriving the same function's own logic).
4. empirical_pairwise_correlation() returns valid correlations in
   [-1, 1] with real p-values, and the three pairs' correlations are
   POSITIVE for both cell lines -- the headline qualitative finding
   (tetraspanins really do co-occur more than chance on individual real
   vesicles, matching the qualitative direction Aim 6's model assumed,
   even though the exact magnitude differs from Aim 6's committed value).
5. average_measured_correlation_across_pairs() is the correct mean of the
   three pairwise values (checked via independent computation).
6. verification_against_figure4_marginals() correctly imports and
   compares against Aim 6's own real Figure 4 data, and correctly reports
   these as DIFFERING (not matching) -- since if this test ever found
   them matching to within a tight tolerance, that would mean something
   about the data provenance had silently changed and needs
   investigation, not celebration.

PART B additionally locks in (see aim8_raw_vesicle_correlation_check.py's
module docstring, Part B, for full provenance of the CRC raw data):
7. The CRC raw data file loads with the expected real sample size
   (24,000 vesicles across 8 real patient samples, 3,000 each).
8. crc_empirical_marginal_positivity() values are bounded and match an
   independent hand recomputation.
9. crc_empirical_pairwise_correlation() returns valid, bounded
   correlations with real p-values.
10. crc_cluster11_like_fraction() is computed correctly (hand-checked)
    and honestly reports a real difference from Aim 7's published
    cluster11 fraction (this is expected -- see docstring -- a simple
    triple-AND threshold is not the same computation as the paper's own
    clustering algorithm).
11. crc_marker_positivity_by_tumor_location() returns valid fractions and
    p-values, and only covers the two real labels this file actually has
    (CO/RE) -- not a healthy-vs-cancer comparison.
"""

import sys
import os

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from aim8_raw_vesicle_correlation_check import (
    RAW_DATA_FILES,
    TETRASPANIN_COLUMNS,
    AIM6_COMMITTED_CORRELATION,
    load_raw_vesicle_data,
    empirical_marginal_positivity,
    empirical_pairwise_correlation,
    average_measured_correlation_across_pairs,
    compare_measured_correlation_to_aim6_committed,
    verification_against_figure4_marginals,
    CRC_MARKERS,
    AIM7_DISCOVERY_CLUSTER11_FRACTION,
    load_crc_raw_vesicle_data,
    crc_empirical_marginal_positivity,
    crc_empirical_pairwise_correlation,
    crc_cluster11_like_fraction,
    crc_marker_positivity_by_tumor_location,
)


# ---------------------------------------------------------------------
# Test 1: raw data loads with the expected real sample size.
# ---------------------------------------------------------------------

def test_raw_data_has_expected_sample_size():
    for cell_line in RAW_DATA_FILES:
        data = load_raw_vesicle_data(cell_line)
        assert data["n_vesicles"] == 3000, (
            f"{cell_line}: expected 3000 real vesicles (confirmed directly "
            f"when this raw file was first inspected), got "
            f"{data['n_vesicles']}. If this number has genuinely changed "
            f"(e.g. the source CSV was regenerated from an updated "
            f"download), update this test deliberately -- do not just "
            f"loosen the check."
        )
        for marker in ("CD9", "CD63", "CD81"):
            assert len(data[marker]) == data["n_vesicles"], (
                f"{cell_line}/{marker}: array length {len(data[marker])} "
                f"does not match n_vesicles {data['n_vesicles']}."
            )
            assert set(np.unique(data[marker])).issubset({0, 1}), (
                f"{cell_line}/{marker}: expected only 0/1 values, got "
                f"{set(np.unique(data[marker]))}."
            )


# ---------------------------------------------------------------------
# Test 2: unknown cell line is rejected clearly.
# ---------------------------------------------------------------------

def test_load_raw_vesicle_data_rejects_unknown_cell_line():
    try:
        load_raw_vesicle_data("ASPC1")
        assert False, (
            "Expected a ValueError for 'ASPC1' -- this raw file only "
            "covers PANC-1/CAPAN-2 (see module docstring)."
        )
    except ValueError:
        pass


# ---------------------------------------------------------------------
# Test 3: empirical marginals are correct and bounded.
# ---------------------------------------------------------------------

def test_empirical_marginal_positivity_matches_hand_computation():
    for cell_line in RAW_DATA_FILES:
        data = load_raw_vesicle_data(cell_line)
        result = empirical_marginal_positivity(cell_line)
        for marker in ("CD9", "CD63", "CD81"):
            hand_computed = float(data[marker].mean())
            assert abs(result[marker] - hand_computed) < 1e-9, (
                f"{cell_line}/{marker}: empirical_marginal_positivity() "
                f"returned {result[marker]}, hand computation from the "
                f"same raw array gives {hand_computed}."
            )
            assert 0.0 <= result[marker] <= 1.0, (
                f"{cell_line}/{marker}: positivity fraction "
                f"{result[marker]} out of bounds [0, 1]."
            )


# ---------------------------------------------------------------------
# Test 4: pairwise correlations are valid and positive (headline
# qualitative finding).
# ---------------------------------------------------------------------

def test_pairwise_correlations_are_bounded():
    for cell_line in RAW_DATA_FILES:
        pairwise = empirical_pairwise_correlation(cell_line)
        for pair, detail in pairwise.items():
            assert -1.0 <= detail["correlation"] <= 1.0, (
                f"{cell_line} {pair}: correlation {detail['correlation']} "
                f"out of bounds [-1, 1]."
            )
            assert 0.0 <= detail["p_value"] <= 1.0, (
                f"{cell_line} {pair}: p_value {detail['p_value']} out of "
                f"bounds [0, 1]."
            )


def test_headline_finding_all_pairs_positively_correlated():
    """
    THE HEADLINE QUALITATIVE FINDING: every one of the three tetraspanin
    pairs, in both real cell lines, shows a REAL, statistically
    significant POSITIVE correlation (not just a positive point estimate
    that could be noise) -- confirming, from genuinely raw single-vesicle
    data, the qualitative assumption Aim 6's model made (tetraspanins
    co-occur more than chance) even though the measured MAGNITUDE (~0.08-
    0.15) differs from Aim 6's committed value (0.2). A significance
    threshold of p < 0.01 is used -- deliberately looser than the actual
    p < 0.005 found for every pair, so this is a real regression guard,
    not a brittle exact-value check.
    """
    for cell_line in RAW_DATA_FILES:
        pairwise = empirical_pairwise_correlation(cell_line)
        for pair, detail in pairwise.items():
            assert detail["correlation"] > 0, (
                f"{cell_line} {pair}: expected a positive correlation, "
                f"got {detail['correlation']:.4f} -- this is the headline "
                f"finding; if this now fails, that is a real and "
                f"important result to report honestly."
            )
            assert detail["p_value"] < 0.01, (
                f"{cell_line} {pair}: correlation {detail['correlation']:.4f} "
                f"is not statistically significant at p<0.01 "
                f"(p={detail['p_value']:.4f})."
            )


# ---------------------------------------------------------------------
# Test 5: average correlation is the correct mean.
# ---------------------------------------------------------------------

def test_average_correlation_is_correct_mean():
    for cell_line in RAW_DATA_FILES:
        pairwise = empirical_pairwise_correlation(cell_line)
        expected_average = float(np.mean([v["correlation"] for v in pairwise.values()]))
        actual_average = average_measured_correlation_across_pairs(cell_line)
        assert abs(actual_average - expected_average) < 1e-9, (
            f"{cell_line}: average_measured_correlation_across_pairs() "
            f"returned {actual_average}, independent hand computation "
            f"gives {expected_average}."
        )


def test_comparison_to_aim6_committed_correlation_is_consistent():
    for cell_line in RAW_DATA_FILES:
        result = compare_measured_correlation_to_aim6_committed(cell_line)
        assert result["aim6_committed_correlation"] == AIM6_COMMITTED_CORRELATION
        expected_diff = abs(
            result["measured_average_correlation"] - AIM6_COMMITTED_CORRELATION
        )
        assert abs(result["absolute_difference"] - expected_diff) < 1e-9, (
            f"{cell_line}: reported absolute_difference does not match "
            f"a direct recomputation."
        )


# ---------------------------------------------------------------------
# Test 6: cross-check against Aim 6's Figure 4 data correctly reports a
# DIFFERENCE (these are two different real datasets, not identical).
# ---------------------------------------------------------------------

def test_figure6a_and_figure4_marginals_are_meaningfully_different():
    """
    This raw file (Figure 6a's source data) and Aim 6's Figure 4 data are
    DIFFERENT real datasets from the same paper (see module docstring for
    the full account -- different N, different gating). This test checks
    that the cross-check function correctly reports a REAL difference
    (not a bug making them look identical, and not a bug making the
    comparison itself fail) -- at least one marker/cell-line combination
    should show a difference of at least 2 percentage points, which is
    what was actually found (up to 13.1pp for PANC-1 CD9). If this test
    ever finds the two datasets suddenly matching almost exactly, that is
    NOT a good sign -- it would mean the raw CSV or Aim 6's data changed
    in a way that needs investigating, not celebrating as "validation
    improved."
    """
    comparison = verification_against_figure4_marginals()
    max_diff = 0.0
    for cell_line, markers in comparison.items():
        for marker, vals in markers.items():
            assert 0.0 <= vals["figure6a_raw_data"] <= 1.0
            assert 0.0 <= vals["figure4_summary_stats"] <= 1.0
            max_diff = max(max_diff, vals["absolute_difference"])
    assert max_diff > 0.02, (
        f"Expected at least one marker/cell-line combination to differ by "
        f"more than 2 percentage points between Figure 6a's raw data and "
        f"Figure 4's summary statistics (this is expected -- they are "
        f"different real datasets), but the largest difference found was "
        f"only {max_diff*100:.2f}pp."
    )


# ---------------------------------------------------------------------
# PART B: CRC integrin tests (NEW)
# ---------------------------------------------------------------------

def test_crc_raw_data_has_expected_sample_size():
    data = load_crc_raw_vesicle_data()
    assert data["n_vesicles"] == 24000, (
        f"Expected 24000 real vesicles (confirmed directly when this raw "
        f"file was first processed: 8 samples x 3000 vesicles each), got "
        f"{data['n_vesicles']}. If this has genuinely changed, update this "
        f"test deliberately -- do not just loosen the check."
    )
    for marker in CRC_MARKERS:
        assert len(data[marker]) == data["n_vesicles"], (
            f"{marker}: array length {len(data[marker])} does not match "
            f"n_vesicles {data['n_vesicles']}."
        )
        assert set(np.unique(data[marker])).issubset({0, 1}), (
            f"{marker}: expected only 0/1 values, got "
            f"{set(np.unique(data[marker]))}."
        )
    assert len(data["sample_code"]) == data["n_vesicles"]
    unique_samples = set(data["sample_code"])
    assert len(unique_samples) == 8, (
        f"Expected 8 real patient samples (CO_01-05, RE_01-03), got "
        f"{len(unique_samples)}: {sorted(unique_samples)}"
    )


def test_crc_empirical_marginal_positivity_matches_hand_computation():
    data = load_crc_raw_vesicle_data()
    result = crc_empirical_marginal_positivity()
    for marker in CRC_MARKERS:
        hand_computed = float(data[marker].mean())
        assert abs(result[marker] - hand_computed) < 1e-9, (
            f"{marker}: crc_empirical_marginal_positivity() returned "
            f"{result[marker]}, hand computation gives {hand_computed}."
        )
        assert 0.0 <= result[marker] <= 1.0, (
            f"{marker}: positivity fraction {result[marker]} out of "
            f"bounds [0, 1]."
        )


def test_crc_pairwise_correlations_are_bounded():
    pairwise = crc_empirical_pairwise_correlation()
    # 4 markers -> 6 unique pairs
    assert len(pairwise) == 6, f"Expected 6 marker pairs, got {len(pairwise)}"
    for pair, detail in pairwise.items():
        assert -1.0 <= detail["correlation"] <= 1.0, (
            f"{pair}: correlation {detail['correlation']} out of bounds "
            f"[-1, 1]."
        )
        assert 0.0 <= detail["p_value"] <= 1.0, (
            f"{pair}: p_value {detail['p_value']} out of bounds [0, 1]."
        )


def test_crc_cluster11_like_fraction_is_correct_and_honestly_different():
    """
    Hand-checks crc_cluster11_like_fraction()'s triple-AND computation
    directly against the raw arrays, and confirms it honestly reports a
    real difference from Aim 7's published cluster11 fraction (4.58%) --
    expected, since this is a simple threshold on raw markers, not a
    reproduction of the paper's own clustering algorithm. If this test
    ever finds the two numbers matching almost exactly, that's not a good
    sign -- it would mean something about the computation changed to
    silently mimic the published number rather than measuring
    independently, and needs investigating.
    """
    data = load_crc_raw_vesicle_data()
    hand_mask = (data["ITGAM"] == 1) & (data["ITGAL"] == 1) & (data["ITGB2"] == 1)
    hand_fraction = float(hand_mask.mean())

    result = crc_cluster11_like_fraction()
    assert abs(result["measured_fraction"] - hand_fraction) < 1e-9, (
        f"crc_cluster11_like_fraction() returned "
        f"{result['measured_fraction']}, hand computation gives "
        f"{hand_fraction}."
    )
    assert result["n_positive"] == int(hand_mask.sum())
    assert result["n_total"] == data["n_vesicles"]
    assert result["aim7_discovery_cluster11_fraction"] == AIM7_DISCOVERY_CLUSTER11_FRACTION
    assert 0.0 <= result["measured_fraction"] <= 1.0

    expected_diff = abs(result["measured_fraction"] - AIM7_DISCOVERY_CLUSTER11_FRACTION)
    assert abs(result["absolute_difference"] - expected_diff) < 1e-9

    # Honest expectation: real measured fraction should meaningfully
    # differ from the published summary-table percentage (different
    # methods -- see docstring). Confirmed in exploratory analysis to be
    # about 0.5% measured vs. 4.58% published.
    assert result["absolute_difference"] > 0.01, (
        f"Expected a real, meaningful difference (>1pp) between the "
        f"measured triple-AND fraction and Aim 7's published cluster11 "
        f"fraction (different methods -- see docstring), but found only "
        f"{result['absolute_difference']*100:.3f}pp difference. This "
        f"could mean the underlying data or comparison changed in a way "
        f"that needs investigating."
    )


def test_crc_marker_positivity_by_tumor_location_is_valid():
    data = load_crc_raw_vesicle_data()
    is_co = np.array([code.startswith("CO") for code in data["sample_code"]])
    is_re = np.array([code.startswith("RE") for code in data["sample_code"]])
    assert is_co.sum() + is_re.sum() == data["n_vesicles"], (
        "Expected every real vesicle to be labeled either CO (colon) or "
        "RE (rectal) -- these are the only two real tumor-location labels "
        "this file contains."
    )

    by_location = crc_marker_positivity_by_tumor_location()
    assert set(by_location.keys()) == set(CRC_MARKERS)
    for marker, detail in by_location.items():
        assert 0.0 <= detail["CO_fraction"] <= 1.0
        assert 0.0 <= detail["RE_fraction"] <= 1.0
        assert 0.0 <= detail["p_value"] <= 1.0
        expected_diff = detail["CO_fraction"] - detail["RE_fraction"]
        assert abs(detail["difference"] - expected_diff) < 1e-9


if __name__ == "__main__":
    # Allow running without pytest, matching this project's other test
    # files' convention.
    test_functions = [
        test_raw_data_has_expected_sample_size,
        test_load_raw_vesicle_data_rejects_unknown_cell_line,
        test_empirical_marginal_positivity_matches_hand_computation,
        test_pairwise_correlations_are_bounded,
        test_headline_finding_all_pairs_positively_correlated,
        test_average_correlation_is_correct_mean,
        test_comparison_to_aim6_committed_correlation_is_consistent,
        test_figure6a_and_figure4_marginals_are_meaningfully_different,
        test_crc_raw_data_has_expected_sample_size,
        test_crc_empirical_marginal_positivity_matches_hand_computation,
        test_crc_pairwise_correlations_are_bounded,
        test_crc_cluster11_like_fraction_is_correct_and_honestly_different,
        test_crc_marker_positivity_by_tumor_location_is_valid,
    ]
    for test_fn in test_functions:
        try:
            test_fn()
            print(f"PASS: {test_fn.__name__}")
        except AssertionError as e:
            print(f"FAIL: {test_fn.__name__}\n  {e}")
