"""Pins the third continuous-expendability attempt (measured cell lifespans, precommitted schemes X1-X6)."""
import json, os, sys
import pytest
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', 'src', 'aim13_scoping')
sys.path.insert(0, SRC)
import external_turnover_weights as X  # noqa: E402
import continuous_weights_v2 as C2  # noqa: E402

R = json.load(open(os.path.join(SRC, 'external_turnover_weights_test.json')))


def test_baselines():
    assert R['baselines']['one_factor'] == pytest.approx(0.880, abs=0.001)
    assert R['baselines']['binary_two_factor'] == pytest.approx(0.914, abs=0.001)


@pytest.mark.parametrize('k,v', [('X1', 0.541), ('X2', 0.498), ('X3', 0.833), ('X4', 0.809), ('X5', 0.766), ('X6', 0.718)])
def test_scheme_separations_pinned(k, v):
    assert R['schemes'][k]['sep'] == pytest.approx(v, abs=0.001)


def test_no_scheme_meets_adoption_rule_so_binary_stays_primary():
    assert R['decision']['binary_stays_primary'] is True
    for v in R['schemes'].values():
        assert v['sep'] < R['baselines']['binary_two_factor']


def test_turnover_adds_nothing_beyond_confinement():
    # Every confinement + turnover scheme is BELOW confinement alone, and none beats its shuffled null.
    for k in ('X3', 'X4', 'X5', 'X6'):
        assert R['schemes'][k]['sep'] < R['baselines']['one_factor']
        assert R['schemes'][k]['p_vs_shuffled_weights'] > 0.05


def test_turnover_alone_is_at_chance():
    for k in ('X1', 'X2'):
        assert R['schemes'][k]['p_label_permutation'] > 0.05


def test_mapping_coverage_pinned():
    assert R['n_types'] == 154 and R['n_types_mapped'] == 77
    assert R['targets_with_mapped_home'] == 19 and R['targets_total'] == 30


def test_lifespan_weights():
    assert X.raw_weight('neurons' if False else 'brain excitatory neurons') == pytest.approx(1.0)
    assert X.raw_weight('colonocytes') == pytest.approx(0.0975, abs=0.001)   # log10(3.4)/log10(3e5)
    assert X.raw_weight('syncytiotrophoblasts') is None                      # unmapped -> neutral weight


def test_live_recompute_matches_json():
    d = json.load(open(os.path.join(SRC, 'wholebody_v2.json')))
    A = d['reference_A_cell_type']
    S = {g: v for g, v in C2.shares().items() if g in A}
    types = sorted({c for t, _ in S.values() for c in t})
    import external_turnover_weights_test as T
    wf = X.weight_fun(types)
    for k in ('X3', 'X5'):
        sc = T.score_all(S, A, k, wf, X.flag_fun)
        s, _, _ = C2.sep_of(sc, d['labels'])
        assert s == pytest.approx(R['schemes'][k]['sep'], abs=1e-9)
