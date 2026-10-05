"""
test_aim7_crc_plasma_exosome_cluster_replication.py

Regression tests for src/aim7_crc_plasma_exosome_cluster_replication.py --
the project's second real external-data check (see that file's own module
docstring for why this is a REPLICATION check, not a predictive
validation like aim6, and what that distinction means).

WHAT THESE TESTS LOCK IN, AND WHY EACH ONE MATTERS
-------------------------------------------------------
1. The real data itself is internally consistent: the discovery cohort's
   12 clusters sum to 100%, and the validation cohort's 15 clusters sum
   to ~100% for BOTH groups (CRC and HC) -- a transcription sanity check,
   the same role aim6's own sum-to-one test plays.
2. clusters_with_marker() correctly finds clusters by marker name,
   checked against a hand-verified expectation from Table S7 (e.g. ITGB3
   appears in both cluster4 and cluster8).
3. cancer_marker_cluster_enrichment() computes the enrichment ratio
   correctly (a direct division, checked against hand computation) and
   the THE HEADLINE FINDING: the validation cohort's cluster11 (the same
   cluster number the discovery cohort found to be ITGAM+/ITGAL+/ITGB2+
   monocyte/macrophage-derived) is enriched in CRC patients relative to
   healthy controls by a real margin, i.e. this replicates directionally
   across the two cohorts.
4. within_cluster_group_composition() returns the real numbers correctly
   and rejects an unknown cluster name.
5. The within-cluster composition (Part 3) and the per-patient-fraction
   framing (Part 2) are checked to be CONSISTENT WITH EACH OTHER in
   DIRECTION (not numerically identical -- they are different framings,
   see module docstring -- but both should agree on whether cluster11 is
   CRC-enriched and cluster3 is HC-enriched, since a real biological
   effect should show the same direction under either framing).
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from aim7_crc_plasma_exosome_cluster_replication import (
    DISCOVERY_COHORT_CLUSTERS,
    VALIDATION_COHORT_CLUSTERS,
    VALIDATION_COHORT_CLUSTERS_COMBINED,
    CLUSTER_GROUP_COMPOSITION,
    discovery_cohort_total_fraction,
    validation_cohort_total_fraction,
    clusters_with_marker,
    cancer_marker_cluster_enrichment,
    within_cluster_group_composition,
)


# ---------------------------------------------------------------------
# Test 1: real data sums to ~1.0.
# ---------------------------------------------------------------------

def test_discovery_cohort_sums_to_one():
    total = discovery_cohort_total_fraction()
    assert abs(total - 1.0) < 1e-9, (
        f"Discovery cohort's 12 clusters sum to {total:.6f}, expected "
        f"exactly 1.0 (Table S7's own percentages sum to exactly 100.0%, "
        f"unlike some of aim6's real data, which has a documented rounding "
        f"artifact -- if this now fails, check the transcribed numbers "
        f"for a typo)."
    )


def test_validation_cohort_sums_to_approximately_one_for_both_groups():
    """
    Tolerance is 2e-3 (0.2 percentage points), matching aim6's own
    precedent for real published one-decimal-place percentages that sum
    to slightly off 100% due to rounding (CRC sums to exactly 100.00%,
    HC to 99.99% here -- both well within this tolerance, checked
    directly against the transcribed legend numbers).
    """
    for group in ("CRC", "HC"):
        total = validation_cohort_total_fraction(group)
        assert abs(total - 1.0) < 2e-3, (
            f"Validation cohort group {group!r}: 15 clusters sum to "
            f"{total*100:.2f}%, expected ~100% within rounding."
        )


def test_validation_cohort_combined_sums_to_approximately_one():
    total = sum(VALIDATION_COHORT_CLUSTERS_COMBINED.values())
    assert abs(total - 1.0) < 2e-3, (
        f"Validation cohort combined-group clusters sum to {total*100:.2f}%, "
        f"expected ~100%."
    )


# ---------------------------------------------------------------------
# Test 2: clusters_with_marker() correctness.
# ---------------------------------------------------------------------

def test_clusters_with_marker_finds_itgb3_in_two_clusters():
    # Hand-checked directly against Table S7: ITGB3 is listed as a
    # defining marker for both cluster4 (ITGB3, ITGA6) and cluster8
    # (ITGB3 alone).
    result = clusters_with_marker("ITGB3")
    assert set(result) == {"cluster4", "cluster8"}, (
        f"Expected ITGB3 to appear in exactly cluster4 and cluster8 per "
        f"Table S7, got {result}."
    )


def test_clusters_with_marker_finds_itga6_in_two_clusters():
    # ITGA6 appears in cluster4 (ITGB3, ITGA6) AND cluster12 (ITGA6 alone)
    # per Table S7.
    result = clusters_with_marker("ITGA6")
    assert set(result) == {"cluster4", "cluster12"}, (
        f"Expected ITGA6 to appear in exactly cluster4 and cluster12 per "
        f"Table S7, got {result}."
    )


def test_clusters_with_marker_returns_empty_for_unused_marker():
    result = clusters_with_marker("CD9")
    assert result == (), (
        f"CD9 is not listed as a defining marker for any discovery-cohort "
        f"cluster in Table S7 -- expected an empty tuple, got {result}."
    )


# ---------------------------------------------------------------------
# Test 3: cancer_marker_cluster_enrichment() -- the headline finding.
# ---------------------------------------------------------------------

def test_enrichment_ratio_is_direct_division():
    result = cancer_marker_cluster_enrichment("cluster11")
    expected_ratio = (
        VALIDATION_COHORT_CLUSTERS["CRC"]["cluster11"]
        / VALIDATION_COHORT_CLUSTERS["HC"]["cluster11"]
    )
    assert abs(result["enrichment_ratio"] - expected_ratio) < 1e-9, (
        f"enrichment_ratio {result['enrichment_ratio']} does not match "
        f"direct hand computation {expected_ratio}."
    )


def test_headline_finding_cluster11_enriched_in_crc():
    """
    THE HEADLINE REPLICATION FINDING: cluster11 in the INDEPENDENT
    validation cohort is enriched in CRC patients relative to healthy
    controls by a real, non-trivial margin. A threshold of 1.5x is used
    (deliberately looser than the actual ~2.08x found) so this is a real
    regression guard against the finding disappearing, not a brittle
    exact-value check.
    """
    result = cancer_marker_cluster_enrichment("cluster11")
    assert result["enrichment_ratio"] > 1.5, (
        f"Expected validation-cohort cluster11 to show enrichment ratio "
        f"> 1.5x (CRC fraction / HC fraction), got "
        f"{result['enrichment_ratio']:.2f}x -- this is the headline "
        f"replication finding; if this now fails against the real "
        f"transcribed numbers, that is a real and important result to "
        f"report honestly, not to paper over."
    )


def test_all_clusters_have_positive_enrichment_ratios():
    """
    Basic sanity check: every enrichment ratio should be a positive
    finite number (both real fractions being divided are themselves
    positive percentages, so this should hold for all 15 clusters by
    construction -- a check against, e.g., an accidental zero or
    negative value from a transcription error).
    """
    for cluster_name in VALIDATION_COHORT_CLUSTERS["CRC"]:
        result = cancer_marker_cluster_enrichment(cluster_name)
        assert result["enrichment_ratio"] > 0, (
            f"{cluster_name}: enrichment ratio {result['enrichment_ratio']} "
            f"is not positive."
        )


# ---------------------------------------------------------------------
# Test 4: within_cluster_group_composition() correctness and error
# handling.
# ---------------------------------------------------------------------

def test_within_cluster_composition_sums_to_approximately_one():
    for cluster_name in CLUSTER_GROUP_COMPOSITION:
        comp = within_cluster_group_composition(cluster_name)
        total = comp["CRC"] + comp["HC"]
        assert abs(total - 1.0) < 1e-6, (
            f"{cluster_name}: CRC + HC composition = {total}, expected 1.0."
        )


def test_within_cluster_composition_rejects_unknown_cluster():
    try:
        within_cluster_group_composition("cluster99")
        assert False, "Expected ValueError for an unknown cluster name."
    except ValueError:
        pass


# ---------------------------------------------------------------------
# Test 5: the two framings (Part 2's per-patient-fraction enrichment vs.
# Part 3's within-cluster composition) agree in DIRECTION for the two
# clusters where both are available (cluster3 and cluster11) -- a real
# cross-check between two independently-reported numbers from the same
# figure, not just re-deriving one from the other (they are NOT
# mathematically equivalent quantities -- see module docstring -- so
# agreement in direction is itself a meaningful consistency check, not a
# tautology).
# ---------------------------------------------------------------------

def test_enrichment_direction_agrees_with_within_cluster_composition():
    for cluster_name in CLUSTER_GROUP_COMPOSITION:
        enrichment = cancer_marker_cluster_enrichment(cluster_name)
        composition = within_cluster_group_composition(cluster_name)

        enrichment_says_crc_enriched = enrichment["enrichment_ratio"] > 1.0
        composition_says_crc_majority = composition["CRC"] > composition["HC"]

        assert enrichment_says_crc_enriched == composition_says_crc_majority, (
            f"{cluster_name}: Part 2 framing (enrichment_ratio="
            f"{enrichment['enrichment_ratio']:.2f}x, "
            f"{'CRC-enriched' if enrichment_says_crc_enriched else 'HC-enriched'}) "
            f"disagrees in DIRECTION with Part 3 framing (within-cluster "
            f"composition {composition['CRC']*100:.1f}% CRC / "
            f"{composition['HC']*100:.1f}% HC, "
            f"{'CRC-majority' if composition_says_crc_majority else 'HC-majority'}) "
            f"-- these are two independently-reported real numbers from the "
            f"same supplementary figure; a direction mismatch would mean a "
            f"transcription error in one of them."
        )


if __name__ == "__main__":
    # Allow running without pytest, matching this project's other test
    # files' convention.
    test_functions = [
        test_discovery_cohort_sums_to_one,
        test_validation_cohort_sums_to_approximately_one_for_both_groups,
        test_validation_cohort_combined_sums_to_approximately_one,
        test_clusters_with_marker_finds_itgb3_in_two_clusters,
        test_clusters_with_marker_finds_itga6_in_two_clusters,
        test_clusters_with_marker_returns_empty_for_unused_marker,
        test_enrichment_ratio_is_direct_division,
        test_headline_finding_cluster11_enriched_in_crc,
        test_all_clusters_have_positive_enrichment_ratios,
        test_within_cluster_composition_sums_to_approximately_one,
        test_within_cluster_composition_rejects_unknown_cluster,
        test_enrichment_direction_agrees_with_within_cluster_composition,
    ]
    for test_fn in test_functions:
        try:
            test_fn()
            print(f"PASS: {test_fn.__name__}")
        except AssertionError as e:
            print(f"FAIL: {test_fn.__name__}\n  {e}")
