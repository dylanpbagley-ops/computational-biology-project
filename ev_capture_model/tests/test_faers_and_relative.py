"""Pins round 6 (section 39): weights diagnostic, measured FAERS severity weights, relative-exposure rules, naming test."""
import hashlib, json, os, sys
import pytest
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', 'src', 'aim13_scoping')
sys.path.insert(0, SRC)
import faers_severity as F  # noqa: E402


def J(name):
    return json.load(open(os.path.join(SRC, name)))


WD, FC, FR, RE, RN = J('weights_diagnostic.json'), J('faers_counts.json'), J('faers_rules.json'), J('relative_exposure.json'), J('relative_naming.json')


def test_system_and_cell_type_mapping_is_complete_and_disjoint():
    assert len(F.SYSTEM_PTS) == 28
    assert len(F.CELL_SYSTEM) == 147 and len(F.UNMAPPED_TYPES) == 7
    assert not set(F.CELL_SYSTEM) & set(F.UNMAPPED_TYPES)
    assert len(F.CELL_SYSTEM) + len(F.UNMAPPED_TYPES) == 154


@pytest.mark.parametrize('tissue,cell,expected', [
    ('Liver', 'hepatocytes', 'liver'), ('Kidney', 'cells in glomeruli', 'kidney_glomerulus'),
    ('Skin', 'melanocytes', 'melanocyte'), ('Adipose tissue', 'adipocytes', None),
    ('Retina', 'photoreceptor cells', 'retina_nerve'), ('Eye', 'corneal epithelial cells', 'eye_surface')])
def test_ihc_system_examples(tissue, cell, expected):
    assert F.ihc_system(tissue, cell) == expected


def test_stored_counts_are_consistent_and_recompute_the_weights():
    s = FC['systems']
    assert set(s) == set(F.SYSTEM_PTS)
    assert (s['liver']['reports'], s['liver']['deaths'], s['liver']['severe']) == (94244, 26518, 35188)
    for k, v in s.items():
        assert 0 < v['deaths'] <= v['severe'] <= v['reports'], k
    w = F.measured_weights()
    assert max(w.values()) == pytest.approx(1.0) and min(w.values()) > 0
    assert max(w, key=w.get) == 'lung' and min(w, key=w.get) == 'eye_surface'
    for k in w:
        assert FR['weights']['death'][k] == pytest.approx(w[k], abs=1e-12), k
        assert FR['weights']['severe'][k] == pytest.approx(s[k]['severe'] / s[k]['reports'] / max(v['severe'] / v['reports'] for v in s.values()), abs=1e-12), k
    assert w['heart_muscle'] == pytest.approx(0.756, abs=0.001) and w['liver'] == pytest.approx(0.737, abs=0.001)


def test_frozen_inputs_still_match_their_precommitted_hashes():
    n_checked = 0
    for stem in ('faers_precommitment.md', 'relative_exposure_precommitment.md', 'weights_diagnostic_precommitment.md'):
        path = os.path.join(SRC, stem + '.hashes')
        assert os.path.exists(path), stem
        for line in open(path):
            if not line.strip():
                continue
            digest, name = line.split(None, 1)
            f = os.path.join(SRC, name.strip())
            assert os.path.exists(f), (stem, name.strip())
            assert hashlib.sha256(open(f, 'rb').read()).hexdigest() == digest, (stem, name.strip())
            n_checked += 1
    assert n_checked >= 15


def test_weights_diagnostic_information_sits_in_the_non_blood_group():
    assert WD['frozen_sep'] == pytest.approx(0.914, abs=0.001)
    assert WD['blood_flag_only_sep'] == pytest.approx(0.711, abs=0.001)
    assert WD['frozen_within_nonblood_sep'] == pytest.approx(0.851, abs=0.001)
    assert (WD['n_nonblood_safe'], WD['n_nonblood_toxic']) == (11, 11)
    assert WD['blood_flagged'] == {'safe': 8, 'toxic': 0}


def test_family_weights_cannot_help_much_and_do_not_generalise():
    for k in ('F', 'G'):
        assert WD[k]['oracle_sep'] - WD['frozen_sep'] < 0.03, k         # hindsight upper bound is tiny
        assert WD[k]['nested_cv_sep'] < WD['frozen_sep'] - 0.10, k      # honest cross-validation is far lower
        assert WD[k]['null_oracle_mean'] > 0.6, k                       # the optimiser fits noise well
    assert WD['oracle_F_misordered_pairs'] == 16
    assert sum(WD['oracle_F_misordered_toxic_counts'].values()) == 16 == sum(WD['oracle_F_misordered_safe_counts'].values())


def test_measured_weights_are_uncorrelated_with_hand_weights():
    from collections import defaultdict
    import numpy as np
    from scipy.stats import spearmanr
    import expendability_weights as EW
    hs = defaultdict(list)
    for cell, s in F.CELL_SYSTEM.items():
        hs[s].append(EW.weight(cell))
    ok = [s for s in F.SYSTEM_PTS if hs[s]]
    assert len(ok) == 27                                          # thyroid has no cell type among the 154
    for k in ('death', 'severe'):
        rho = spearmanr([FR['weights'][k][s] for s in ok], [np.mean(hs[s]) for s in ok])[0]
        assert abs(rho) < 0.10, (k, rho)
    assert spearmanr([FR['weights']['death'][s] for s in F.SYSTEM_PTS], [FR['weights']['severe'][s] for s in F.SYSTEM_PTS])[0] > 0.85


def test_thirteen_measured_severity_variants_none_promoted_all_below_frozen():
    v = FR['variants']
    assert len(v) == 13 and FR['n_panel_with_ihc'] == 23
    assert FR['decision']['promoted'] == []
    for k, x in v.items():
        assert x['sep'] < x['frozen_same_genes'], k                     # every variant scores below the frozen rule
        assert x['p_familywise'] > 0.05, k
    for k in ('A3_death', 'A3_severe', 'A6_death', 'A6_severe', 'A7_death', 'A7_severe', 'A5_hand'):
        assert v[k]['diff_ci'][1] < 0, k                                # clearly below: paired interval excludes 0
    for k in ('A2_death', 'A2_severe'):
        assert v[k]['diff_ci'][1] > 0, k                                # statistically tied
    for k in ('A1_death', 'A2_death', 'A3_death', 'A4_death', 'A6_death', 'A7_death'):
        assert v[k]['p_shuffled_weights'] > 0.3, k                      # measured weights look like arbitrary weights


def test_relative_exposure_predictions_p1_p2_p3():
    v, fr = RE['variants'], RE['frozen_sep']
    assert len(v) == 9 and RE['decision']['promoted'] == []
    assert abs(v['binary_sum']['sep'] - fr) < 0.05                                            # P1
    for k in ('hand_sum', 'death_sum', 'severe_sum'):                                          # P2
        assert v[k]['sep'] - v['binary_sum']['sep'] < 0.03, k
    assert max(x['sep'] for x in v.values()) < 0.944 and RE['null_family_p95'] > 0.90         # P3
    assert v['binary_sum']['sep'] == pytest.approx(0.895, abs=0.001)
    assert v['binary_sum']['diff_ci'][0] < 0 < v['binary_sum']['diff_ci'][1]                   # tied with the frozen rule
    assert v['binary_sum']['p_familywise'] > 0.05                                             # a concordance check, not a discovery


def test_worst_organ_form_is_clearly_worse_and_sum_form_is_better():
    v = RE['variants']
    assert v['death_max']['sep'] == pytest.approx(0.694, abs=0.001)
    assert v['death_max']['diff_ci'][1] < 0                                                    # my pre-declared 'most realistic' rule lost
    for w in ('binary', 'hand', 'death', 'severe'):
        assert v[w + '_sum']['sep'] - v[w + '_max']['sep'] > 0.14, w


def test_floor_t0_is_untestable_on_this_panel():
    vals = set(RE['robustness_T0'][k] for k in RE['robustness_T0'] if k.endswith('_death_max'))
    assert len(vals) == 1                                                                      # identical at T0 = 3, 10, 30


def test_system_weight_oracle_overfits():
    for form in ('sum', 'max'):
        o = RE['diagnostic_system_weight_oracle'][form]
        assert o['oracle_sep'] > 0.95 and o['nested_cv_sep'] < RE['frozen_sep'] - 0.10, form


def test_naming_test_graded_weights_help_locate_damage_but_not_significantly():
    s = RN['summary']
    assert RN['n'] == 13
    assert s['E only']['median_best_rank'] == 9.0
    for k in ('hand', 'death', 'severe'):
        assert s[k]['median_best_rank'] < 6 and s[k]['sign_p'] > 0.05, k
    assert s['binary']['median_best_rank'] >= 8.0 and s['binary']['sign_p'] == 1.0
    assert RN['per_target']['F3']['E only'] == 8.0 and RN['per_target']['F3']['death'] == 49.5
