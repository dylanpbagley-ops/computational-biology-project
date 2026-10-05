"""
test_remaining_literature_table_scoping.py

NEW (v14->v15). Regression test for
src/remaining_literature_table_scoping.py -- "Immediate next steps" item
5's checked finding that none of the three remaining literature-table
systems (SK-OV-3, HeLa/HepG2/H1650/SH-SY5Y, MDA-MB-231) currently has the
specific pair of summary statistics (a single-marker positivity fraction
AND a mean copies-per-positive-particle number) the range-sweep
reconstruction template needs.

This is a lightweight structural test, not a numerical one -- there is no
computation to independently re-derive here (this file is a documented,
checked LITERATURE finding, not a model output), but the STRUCTURE of the
finding (which systems are usable, which aren't, and why) should not
silently drift if someone edits the scoping file later without updating
this test to match.

HOW TO RUN
    pytest test_remaining_literature_table_scoping.py -v
or:
    python test_remaining_literature_table_scoping.py
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from remaining_literature_table_scoping import REMAINING_SYSTEMS_SCOPING_RESULT


def test_all_three_remaining_systems_are_scoped():
    """
    All three systems named in "Immediate next steps" item 5 (SK-OV-3,
    HeLa/HepG2/H1650/SH-SY5Y, MDA-MB-231) must have a scoping entry -- if
    a future session adds a fourth system or renames one of these, this
    test should be updated deliberately, not silently pass with a missing
    entry.
    """
    expected_systems = {"SK-OV-3", "HeLa_HepG2_H1650_SH-SY5Y", "MDA-MB-231"}
    assert set(REMAINING_SYSTEMS_SCOPING_RESULT.keys()) == expected_systems, (
        f"Expected exactly the three systems named in 'Immediate next "
        f"steps' item 5 -- got {set(REMAINING_SYSTEMS_SCOPING_RESULT.keys())}."
    )


def test_none_of_the_three_are_currently_usable_for_reconstruction():
    """
    Locks in this session's checked finding: as of v14->v15, none of the
    three remaining systems has both required statistics. If this ever
    becomes True for one of them (e.g. a future session finds a companion
    paper with the missing copy-number data), this test SHOULD fail and
    prompt an update -- that would be good news worth updating this test
    and the handoff to reflect, not a regression to silently patch around.
    """
    for system, info in REMAINING_SYSTEMS_SCOPING_RESULT.items():
        assert info["usable_for_range_sweep_reconstruction"] is False, (
            f"{system} is now marked usable for range-sweep reconstruction "
            f"-- if this is a genuine new finding (e.g. a new source with "
            f"both required statistics was found), update this test "
            f"deliberately and build the actual reconstruction; do not "
            f"just adjust the assertion to make it pass without also "
            f"doing that work."
        )


def test_every_entry_states_both_required_statistics_explicitly():
    """
    Every scoping entry must explicitly record whether it has each of the
    two required statistics (positivity fraction, mean copies per
    positive particle) -- this is the specific, checkable pair the
    reconstruction method needs (see module docstring), and every entry
    should be explicit about both rather than only stating an overall
    conclusion.
    """
    for system, info in REMAINING_SYSTEMS_SCOPING_RESULT.items():
        assert "has_positivity_fraction_for_single_marker" in info
        assert "has_mean_copies_per_positive_particle" in info
        assert "reason" in info and len(info["reason"]) > 20, (
            f"{system}'s scoping entry should have a substantive, "
            f"specific reason recorded, not a placeholder."
        )


def test_hela_hepg2_h1650_shsy5y_has_positivity_but_not_copy_number():
    """
    Specifically locks in the most nuanced of the three findings: the
    Breitwieser et al. 2022 (PMC9369185) paper DOES report a usable
    positivity-fraction analogue (fraction of particles per antibody
    spot) for all four cell lines, but has NO copy-number data anywhere.
    This is different from the other two systems (which fail BOTH
    requirements) and worth keeping distinct in case a future session
    only checks the overall "usable" boolean and misses that this system
    is one missing statistic away from being usable, not two.
    """
    info = REMAINING_SYSTEMS_SCOPING_RESULT["HeLa_HepG2_H1650_SH-SY5Y"]
    assert info["has_positivity_fraction_for_single_marker"] is True, (
        "Expected this system to have a usable positivity-fraction "
        "analogue (checked directly against the full retrieved paper "
        "text this session) -- if this is now False, re-check against "
        "the actual paper before changing it."
    )
    assert info["has_mean_copies_per_positive_particle"] is False, (
        "Expected this system to still lack any copy-number figure -- if "
        "a future session finds one (e.g. in a companion paper), this "
        "should change to True AND the range-sweep reconstruction should "
        "actually be built, not just this flag flipped."
    )


if __name__ == "__main__":
    tests = [
        test_all_three_remaining_systems_are_scoped,
        test_none_of_the_three_are_currently_usable_for_reconstruction,
        test_every_entry_states_both_required_statistics_explicitly,
        test_hela_hepg2_h1650_shsy5y_has_positivity_but_not_copy_number,
    ]

    n_passed = 0
    n_failed = 0
    for test_func in tests:
        print(f"\nRunning {test_func.__name__} ...")
        try:
            test_func()
            print(f"PASSED: {test_func.__name__}")
            n_passed += 1
        except AssertionError as e:
            print(f"FAILED: {test_func.__name__}\n   {e}")
            n_failed += 1

    print(f"\n{'='*60}\n{n_passed} passed, {n_failed} failed out of {len(tests)} tests\n{'='*60}")
    if n_failed > 0:
        sys.exit(1)
