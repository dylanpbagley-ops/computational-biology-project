"""Pins the second continuous-expendability attempt (precommitted schemes E1-E4)."""
import json, os, sys
import pytest
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', 'src', 'aim13_scoping')
sys.path.insert(0, SRC)
import continuous_weights_v2 as C  # noqa: E402

R = json.load(open(os.path.join(SRC, 'continuous_weights_v2.json')))


def test_baselines():
    b = R['baselines']
    assert b['one_factor'] == pytest.approx(0.880, abs=0.001)
    assert b['binary_two_factor'] == pytest.approx(0.914, abs=0.001)
    assert b['old_continuous_home'] == pytest.approx(0.809, abs=0.001)


@pytest.mark.parametrize('k,v', [('E1', 0.818), ('E2', 0.699), ('E3', 0.885), ('E4', 0.880)])
def test_scheme_separations_pinned(k, v):
    assert R['schemes'][k]['sep'] == pytest.approx(v, abs=0.001)


def test_no_scheme_meets_adoption_rule_so_binary_stays_primary():
    # If this fails, a continuous scheme now beats binary AND its shuffled-weight null: revisit the paper wording.
    assert R['decision']['binary_stays_primary'] is True
    for v in R['schemes'].values():
        assert v['sep'] < R['baselines']['binary_two_factor']


def test_soft_expendability_adds_nothing_beyond_confinement():
    # E3/E4 sit at the one-factor level and do not beat their shuffled-weight nulls.
    for k in ('E3', 'E4'):
        assert R['schemes'][k]['sep'] - R['baselines']['one_factor'] < 0.01
        assert R['schemes'][k]['p_vs_shuffled_weights'] > 0.05


def test_live_recompute_matches_json():
    d = json.load(open(os.path.join(SRC, 'wholebody_v2.json')))
    S = {g: v for g, v in C.shares().items() if g in d['reference_A_cell_type']}
    for k in ('E1', 'E3'):
        s, _, _ = C.sep_of(C.scores(S, k), d['labels'])
        assert s == pytest.approx(R['schemes'][k]['sep'], abs=1e-9)
