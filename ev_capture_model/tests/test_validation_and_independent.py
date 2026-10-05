"""Pins round 5: validation targets, label-flip sensitivity, naming test, threshold-harm rules."""
import json, os, sys
import pytest
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', 'src', 'aim13_scoping')
sys.path.insert(0, SRC)
import target_labels as TL  # noqa: E402
import target_labels_validation as V  # noqa: E402

VT = json.load(open(os.path.join(SRC, 'validation_test.json')))
VA = json.load(open(os.path.join(SRC, 'validation_addendum.json')))
IT = json.load(open(os.path.join(SRC, 'independent_tests.json')))
TH = json.load(open(os.path.join(SRC, 'threshold_harm.json')))


def test_validation_labels_do_not_overlap_the_panel_and_tiers_are_as_frozen():
    panel = {t[0] for t in TL.TARGETS}
    assert not ({t[0] for t in V.VALIDATION} & panel)
    assert V.PRIMARY_SAFE == ['MAGEA4', 'SSTR2', 'CCR4', 'CTAG1B', 'IL2RA'] and V.PRIMARY_TOXIC == ['GUCY2C', 'PSCA']


def test_validation_result_pinned():
    p = VT['primary']
    assert p['sep_rule'] == pytest.approx(0.875, abs=0.001) and p['p_exact_permutation'] == pytest.approx(0.133, abs=0.001)
    assert all(p['toxic_below_median_safe'].values())
    assert VA['sep_with_absent_gene_as_safest'] == pytest.approx(0.900, abs=0.001)
    assert VT['missing_genes'] == ['CTAG1B', 'MAGEA12']


@pytest.mark.parametrize('k,v', [('F0 none', 0.914), ('F1 NECTIN4->toxic', 0.926), ('F2 F3,FOLR1->safe', 0.884),
                                 ('F3 F1+F2', 0.890), ('F4 TNFRSF17->toxic', 0.838)])
def test_label_flip_sensitivity(k, v):
    assert VT['flip_sensitivity'][k]['sep'] == pytest.approx(v, abs=0.001)
    assert VT['flip_sensitivity'][k]['p'] < 0.001


def test_naming_summary():
    s = IT['naming']['summary']
    assert (s['N1']['median_best_rank'], s['N3']['median_best_rank']) == (7.0, 4.0)
    assert s['N3']['better_than_N1'] == 8 and s['N3']['worse_than_N1'] == 2
    assert s['N3']['sign_test_p'] == pytest.approx(0.109, abs=0.001)
    assert s['N4']['sign_test_p'] > 0.5
    assert IT['naming']['per_target']['CD47']['N2'] == 137.0     # the informative failure


def test_vital_share_candidates_not_adopted():
    assert IT['rules']['R1']['sep'] == pytest.approx(0.435, abs=0.001) and IT['rules']['R2']['sep'] == pytest.approx(0.507, abs=0.001)
    assert not IT['decision']['R1']['adopt'] and not IT['decision']['R2']['adopt']


@pytest.mark.parametrize('k,v', [('TH1(1)', 0.804), ('TH2(1)', 0.785), ('TH1(10)', 0.856), ('TH2(10)', 0.837), ('TH1(50)', 0.859), ('TH2(50)', 0.833)])
def test_threshold_harm_pinned_and_not_promoted(k, v):
    assert TH['schemes'][k]['sep'] == pytest.approx(v, abs=0.001)
    assert TH['schemes'][k]['sep'] < TH['schemes']['frozen rule']['sep']
    assert TH['decision']['promoted'] == []
