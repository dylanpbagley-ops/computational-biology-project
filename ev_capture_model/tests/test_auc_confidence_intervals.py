"""Tests for the bootstrap confidence intervals on headline separations.

These lock in three things: that the bootstrap's point estimates still reproduce the
separations documented in the results document, that the intervals behave sanely, and --
specifically -- that reference B's haematopoietic flag is evaluated against reference B's
own tissue-level homes rather than borrowed from reference A. That borrowing was a real
bug caught during implementation: it silently inflated B_two_factor from 0.751 to 0.866,
because reference A's cell-type home fires the flag on far more targets than B's
whole-tissue home does. The test below would have caught it immediately.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'aim13_scoping'))

import auc_confidence_intervals as CI  # noqa: E402


@pytest.fixture(scope='module')
def results():
    return CI.compute_all()


def test_point_estimates_match_documented_separations(results):
    """The bootstrap must not change the point estimates -- they are the same statistic."""
    assert results['A_one_factor_top_share']['point'] == pytest.approx(0.8804, abs=1e-3)
    assert results['A_two_factor']['point'] == pytest.approx(0.9139, abs=1e-3)
    assert results['B_two_factor']['point'] == pytest.approx(0.7512, abs=1e-3)


def test_reference_B_uses_its_own_tissue_home_not_reference_A(results):
    """Regression test for the borrowed-flag bug.

    If reference B's expendability flag were taken from reference A's cell-type homes,
    B_two_factor would come out near 0.866 instead of 0.751. Pin the correct value and
    assert we are nowhere near the wrong one.
    """
    b = results['B_two_factor']['point']
    assert b == pytest.approx(0.7512, abs=1e-3)
    assert abs(b - 0.866) > 0.05, 'reference B appears to be using reference A flags'


def test_intervals_bracket_their_point_estimates(results):
    for k, v in results.items():
        if k.startswith('_'):
            continue
        assert v['ci_lo'] <= v['point'] <= v['ci_hi'], f'{k}: point outside its own CI'


def test_intervals_stay_inside_the_unit_interval(results):
    for k, v in results.items():
        if k.startswith('_'):
            continue
        assert 0.0 <= v['ci_lo'] <= 1.0
        assert 0.0 <= v['ci_hi'] <= 1.0


def test_intervals_are_wide_at_this_sample_size(results):
    """Honesty check: at 19v11 these intervals should NOT look precise. If a future change
    makes them narrow, something has gone wrong (e.g. resampling cell types instead of
    targets, which would fabricate precision the design does not support)."""
    for k in ('A_one_factor_top_share', 'A_two_factor', 'B_two_factor'):
        width = results[k]['ci_hi'] - results[k]['ci_lo']
        assert width > 0.15, f'{k}: CI width {width:.3f} is implausibly narrow for n=19v11'


def test_reference_B_lower_bound_reaches_toward_chance(results):
    """The independent-reference replication is significant but its effect size is only
    weakly constrained -- the lower bound sits near 0.5. This is a real limitation and the
    write-up must not present 0.751 as a firm estimate."""
    assert results['B_two_factor']['ci_lo'] < 0.60


def test_bootstrap_is_deterministic_under_fixed_seed():
    a = [0.9, 0.8, 0.7, 0.6, 0.55]
    b = [0.4, 0.3, 0.2, 0.1]
    first = CI.bootstrap_ci(a, b, n_boot=500, seed=7)
    second = CI.bootstrap_ci(a, b, n_boot=500, seed=7)
    assert first == second


def test_perfectly_separated_groups_give_interval_at_one():
    a = [10.0] * 6
    b = [1.0] * 6
    pt, lo, hi, _ = CI.bootstrap_ci(a, b, n_boot=500, seed=0)
    assert pt == pytest.approx(1.0)
    assert lo == pytest.approx(1.0)
    assert hi == pytest.approx(1.0)


def test_identical_groups_centre_on_chance():
    vals = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    pt, lo, hi, _ = CI.bootstrap_ci(vals, list(vals), n_boot=2000, seed=0)
    assert pt == pytest.approx(0.5, abs=1e-9)
    assert lo < 0.5 < hi
