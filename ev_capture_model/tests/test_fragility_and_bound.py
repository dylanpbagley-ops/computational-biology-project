"""Pins section 40: label fragility (exhaustive label swaps) and the top-share breadth lemma."""
import hashlib, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', 'src', 'aim13_scoping')
sys.path.insert(0, SRC)
import label_fragility as LF  # noqa: E402


def J(n):
    return json.load(open(os.path.join(SRC, n)))


FR, BL = J('label_fragility.json'), J('breadth_bound_lemma.json')


def test_frozen_rule_reproduced_from_unrounded_scores():
    genes, score, lab = LF.load()
    auc, p = LF.rank_stats(score, lab == 1)
    assert abs(auc - 0.9138755980861244) < 1e-12 and p < 1e-4
    assert len(genes) == 30 and int(lab.sum()) == 19


def test_exact_mann_whitney_table_is_a_distribution():
    t = LF.exact_table(5, 4)
    assert abs(t[0] - 1.0) < 1e-12 and abs(t[-1] - 1 / 126) < 1e-12   # C(9,4) = 126 arrangements


def test_fragility_index_and_k_tables():
    assert FR['fragility_index'] == 4
    for k in ('1', '2', '3'):
        assert FR['F1'][k]['frac_p_lt_0.05'] == 1.0
    assert 0.98 < FR['F1']['4']['frac_p_lt_0.05'] < 0.99
    assert [FR['F1'][str(k)]['n_sets'] for k in range(1, 7)] == [30, 435, 4060, 27405, 142506, 593775]


def test_every_significance_losing_four_swap_moves_an_approved_blood_target():
    x = FR['k4_sets_losing_significance']
    assert x['n'] == 349 and x['n_containing_a_blood_safe_target'] == 349 and len(x['blood_safe_targets']) == 8


def test_weak_label_subsets():
    s = FR['F2_summary']
    assert len(FR['F2_all_subsets']) == 256
    assert abs(s['all_eight']['sep'] - 0.744) < 1e-9 and 0.04 < s['all_eight']['p'] < 0.05
    assert s['by_n_flipped_max_p']['4'] < 0.05 < s['by_n_flipped_max_p']['5']


def test_breadth_lemma_has_no_violations_and_is_tight():
    assert BL['violations_panel'] == 0 and BL['violations_random'] == 0
    assert BL['tightness_all_within_1'] and BL['n_checks_panel'] == 132
    assert all(v < -0.9 for v in BL['spearman_top_share_vs_E'].values())


def test_precommit_hashes():
    for line in open(os.path.join(SRC, 'label_fragility_precommitment.md.hashes')):
        h, name = line.split()
        assert hashlib.sha256(open(os.path.join(SRC, name), 'rb').read()).hexdigest() == h, name
