"""Pins round 4 (blood-flag cleanup and the search for a more effective rule)."""
import json, os, sys
import pytest
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', 'src', 'aim13_scoping')
sys.path.insert(0, SRC)
import lineage_map as LM  # noqa: E402
import target_labels as TL  # noqa: E402

A = json.load(open(os.path.join(SRC, 'flag_cleanup.json')))
B = json.load(open(os.path.join(SRC, 'rule_improvement.json')))
C = json.load(open(os.path.join(SRC, 'rule_improvement_c.json')))


def test_m2_reproduces_home_rule_on_every_panel_target():
    for g, t in A['targets'].items():
        assert t['M2_flag'] == t['home_flag'], g
    assert A['schemes']['M2']['sep'] == pytest.approx(0.914, abs=0.001)


def test_m1_pinned_and_loses_cd19_cd22():
    assert A['schemes']['M1']['sep'] == pytest.approx(0.904, abs=0.001)
    assert not A['targets']['CD19']['M1_flag'] and not A['targets']['CD22']['M1_flag']


def test_cleanup_adopted_and_no_toxic_target_flagged():
    assert A['decision']['cleanup_adopted'] is True
    for k in ('home_rule', 'M1', 'M2'):
        assert A['flag_counts'][k]['toxic_flagged'] == 0


def test_no_scheme_promoted_in_either_search_round():
    assert B['decision']['promoted'] == [] and C['decision']['promoted'] == []
    for R in (B, C):
        best = max(v['sep'] for v in R['schemes'].values())
        assert best == pytest.approx(0.914, abs=0.001)   # nothing beats the current rule


def test_tau_is_rank_identical_to_top_share():
    assert B['schemes']['C3']['sep'] == B['schemes']['C0']['sep']
    assert B['schemes']['C3+flag']['sep'] == B['schemes']['C0+flag']['sep']


@pytest.mark.parametrize('k,v', [('C1', 0.545), ('C2', 0.794), ('C4', 0.818), ('C5', 0.876), ('C6', 0.761), ('C7', 0.622), ('C8', 0.856)])
def test_round_4b_separations(k, v):
    assert B['schemes'][k]['sep'] == pytest.approx(v, abs=0.001)


@pytest.mark.parametrize('k,v', [('V1', 0.775), ('V2', 0.703), ('V3', 0.699), ('V4', 0.861), ('D1', 0.775), ('D2', 0.854), ('D3', 0.823)])
def test_round_4c_separations(k, v):
    assert C['schemes'][k]['sep'] == pytest.approx(v, abs=0.001)


def test_lineage_map_partition():
    assert LM.lineage('b-cells') == 'blood' and LM.lineage('colonocytes') == 'epithelial'
    assert LM.lineage('oocytes') == 'germline' and LM.lineage('hematopoietic stem cells') == 'blood_other'
    assert TL.haematopoietic('hematopoietic stem cells') is False   # frozen list unchanged
