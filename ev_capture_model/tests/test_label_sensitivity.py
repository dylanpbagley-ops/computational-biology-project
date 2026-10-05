"""Tests for the label-sensitivity analysis.

Locks in the claim the write-up will make: the headline separation survives every
adversarial relabeling tested, and remains significant throughout.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'aim13_scoping'))

import label_sensitivity as LS  # noqa: E402


@pytest.fixture(scope='module')
def results():
    return LS.compute_all()


def test_frozen_scenario_reproduces_documented_numbers(results):
    """Scenario 'as_frozen' must reproduce section 23 exactly -- if it does not, the
    sensitivity analysis is not testing the same thing the results document reports."""
    f = results['as_frozen']
    assert f['one_factor']['sep'] == pytest.approx(0.8804, abs=1e-3)
    assert f['two_factor']['sep'] == pytest.approx(0.9139, abs=1e-3)
    assert f['one_factor']['n_safe'] == 19
    assert f['one_factor']['n_toxic'] == 11


def test_every_scenario_stays_in_the_reported_band(results):
    """The write-up will cite a 0.835-0.926 range. Pin it."""
    s = results['_summary']
    assert s['sep_min'] == pytest.approx(0.835, abs=0.005)
    assert s['sep_max'] == pytest.approx(0.926, abs=0.005)


def test_every_scenario_remains_significant(results):
    """The whole point of the analysis: no relabeling tested breaks significance."""
    assert results['_summary']['p_max'] < 0.01
    for name in LS.SCENARIOS:
        for k in ('one_factor', 'two_factor'):
            assert results[name][k]['p'] < 0.01, f'{name}/{k} lost significance'


def test_two_factor_beats_one_factor_in_every_scenario(results):
    """Expendability should add signal regardless of which labels are perturbed."""
    for name in LS.SCENARIOS:
        one = results[name]['one_factor']['sep']
        two = results[name]['two_factor']['sep']
        assert two > one, f'{name}: two-factor ({two:.3f}) did not beat one-factor ({one:.3f})'


def test_relabelings_actually_move_the_group_sizes(results):
    """Guard against a scenario silently doing nothing (e.g. a typo'd gene symbol that
    matches no target, which would make the analysis look robust for the wrong reason)."""
    frozen = results['as_frozen']['one_factor']
    for name, changes in LS.SCENARIOS.items():
        if not changes:
            continue
        r = results[name]['one_factor']
        moved = (r['n_safe'] != frozen['n_safe']) or (r['n_toxic'] != frozen['n_toxic'])
        assert moved, f'{name}: relabeling changed no group size -- gene symbol typo?'


def test_all_perturbed_genes_exist_in_the_label_file(results):
    """A changed gene that is not in the panel would be a silent no-op."""
    import target_labels as TL
    known = {r[0] for r in TL.TARGETS}
    for name, changes in LS.SCENARIOS.items():
        for gene in changes:
            assert gene in known, f'{name}: {gene} is not in target_labels.TARGETS'


def test_summary_records_the_taxonomy_caveat(results):
    """The caveat is load-bearing: this analysis perturbs labels within the taxonomy, it
    does not validate the taxonomy. Keep it attached to the numbers."""
    assert 'taxonomy' in results['_summary']['caveat'].lower()
