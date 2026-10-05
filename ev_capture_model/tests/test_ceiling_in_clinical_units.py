"""Tests for the ceiling expressed in clinical units.

The load-bearing test here is `test_roc_curve_reproduces_published_aucs`: if the ROC curve
built by this module is correct, the area under it must match the AUCs in
myeloma_ceiling_results.json, which were computed by entirely separate code (a rank-sum
formula, not a curve). That cross-check is what licenses reading any other quantity off
the same curve.
"""
import json
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, '..', 'src', 'aim13_scoping'))

import ceiling_in_clinical_units as C  # noqa: E402

PUBLISHED = os.path.join(HERE, '..', 'src', 'aim13_scoping', 'myeloma_ceiling_results.json')


@pytest.fixture(scope='module')
def gated():
    return C.load_populations()


def test_roc_curve_reproduces_published_aucs(gated):
    """Integrating this module's ROC must recover the independently-computed AUCs."""
    names, ix, samples = gated
    pub = json.load(open(PUBLISHED))
    checked = 0
    for n in ('CD269', 'CD38', 'CD319', 'CD307'):
        for lab in ('T', 'MY', 'NK', 'NB'):
            vals = [C.roc_auc_from_curve(s['tumor'][:, ix[n]], s['comps'][lab][:, ix[n]],
                                         n_grid=501)
                    for s in samples if lab in s['comps']]
            if not vals:
                continue
            assert float(np.median(vals)) == pytest.approx(pub[n][lab], abs=2e-3), \
                f'{n}/{lab}: ROC area disagrees with published AUC'
            checked += 1
    assert checked >= 12, 'too few AUC cross-checks actually ran'


def test_max_kill_is_monotone_in_harm_budget():
    """Allowing more harm can never permit less killing."""
    rng = np.random.default_rng(0)
    tumor = rng.normal(1.0, 1.0, 4000)
    normal = rng.normal(0.0, 1.0, 4000)
    k1 = C.max_kill_at_harm(tumor, normal, 0.01)
    k5 = C.max_kill_at_harm(tumor, normal, 0.05)
    k20 = C.max_kill_at_harm(tumor, normal, 0.20)
    assert k1 <= k5 <= k20


def test_identical_distributions_give_kill_equal_to_budget():
    """With no separation, you can only kill as much tumour as you kill normal."""
    rng = np.random.default_rng(1)
    x = rng.normal(0, 1, 20000)
    y = rng.normal(0, 1, 20000)
    assert C.max_kill_at_harm(x, y, 0.05) == pytest.approx(0.05, abs=0.01)


def test_perfect_separation_allows_total_kill():
    assert C.max_kill_at_harm(np.full(500, 10.0), np.full(500, 0.0), 0.01) == pytest.approx(1.0)


def test_gate_markers_are_decircularised(gated):
    """CD138 and CD269 define the tumour gate, so each must be scored against a gate
    rebuilt without it. Regression guard: the circular version inflates both markedly."""
    r = C.compute_all()
    for m in C.GATE_MARKERS:
        assert r[m]['0.01']['gate_decircularised'] is True
    # circular CD269 at 1% harm was 1.6%; de-circularised it must be materially lower
    assert r['CD269']['0.01']['ceiling'] < 0.014


def test_nonexpendable_ceiling_is_never_below_uniform_ceiling():
    """Dropping a constraint cannot reduce what is achievable."""
    r = C.compute_all()
    for n, v in r.items():
        if n.startswith('_'):
            continue
        for b, d in v.items():
            if d.get('ceiling_nonexpendable_only') is None:
                continue
            assert d['ceiling_nonexpendable_only'] >= d['ceiling'] - 1e-12, \
                f'{n}/{b}: relaxing the expendable constraint lowered the ceiling'


def test_all_values_are_fractions():
    r = C.compute_all()
    for n, v in r.items():
        if n.startswith('_'):
            continue
        for b, d in v.items():
            assert 0.0 <= d['ceiling'] <= 1.0
            for k, x in d['per_compartment'].items():
                assert 0.0 <= x <= 1.0


def test_meta_records_the_calibration_invariance_claim():
    r = C.compute_all()
    note = r['_meta']['units_note'].lower()
    assert 'fraction' in note and 'calibrat' in note
