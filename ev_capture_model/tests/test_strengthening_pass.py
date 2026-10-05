"""Tests for the September 28 strengthening pass.

These pin results -- INCLUDING the unflattering ones -- so a later change that quietly alters them
fails loudly. Where a test pins a negative or a post-hoc result, the docstring says so.
"""
import json
import os
import sys

import numpy as np
import pytest

SRC = os.path.join(os.path.dirname(__file__), '..', 'src', 'aim13_scoping')
sys.path.insert(0, SRC)

import outcome_test_robustness as OR      # noqa: E402
import candidate_hygiene as HY            # noqa: E402
import screen_rank_stability as SR        # noqa: E402
import angle5_monotone_premise as A5      # noqa: E402


def J(name):
    return json.load(open(os.path.join(SRC, name)))


# ------------------------------------------------------------- outcome-test robustness
def test_replication_groups_match_live_panel():
    (aa, ab), (ra, rb) = OR.replication_groups()
    # live panel: 9 v 9 (the results document still says 9 v 8 -- see the drift report)
    assert (len(aa), len(ab)) == (9, 9)
    assert (len(ra), len(rb)) == (11, 13)


def test_power_is_monotone_and_calibrated_at_null():
    lo = OR.power_at_auc(0.5, 9, 11, n_sim=1500)
    hi = OR.power_at_auc(0.88, 9, 11, n_sim=1500)
    assert lo < 0.10          # false-positive rate near alpha = 0.05
    assert hi > 0.90          # the control WOULD have caught an effect the size of the headline


def test_specificity_control_cannot_exclude_small_effects():
    """Pinned honest limit: at AUC 0.65 the registry control detects the effect < 40% of the time."""
    r = J('outcome_test_robustness.json')
    assert r['power']['registry_control_phase3']['power_at_auc']['0.65'] < 0.40


def test_expendability_association_is_significant_but_auc_increment_is_small():
    e = J('outcome_test_robustness.json')['expendability_direct_tests']
    assert e['fisher_exact_one_sided_p'] < 0.02
    assert e['perm_p_increment_at_least_observed'] < 0.02
    assert 0.02 < e['observed_auc_increment'] < 0.05
    assert set(e['flagged_safe_targets_that_confinement_alone_ranks_below_some_toxic_target']) == \
        {'CD22', 'CD33', 'CD38'}


def test_paired_bootstrap_ci_touches_zero():
    """Pinned: the AUC gain from the second factor is NOT established by bootstrap (CI lower = 0)."""
    p = J('outcome_test_robustness.json')['paired_two_factor_minus_one_factor']
    assert p['ci_lo'] <= 0.005
    assert p['point_difference'] == pytest.approx(0.0335, abs=2e-3)


# ------------------------------------------------------------- candidate hygiene
def _rec(cls, secretome='', rel='Approved'):
    return {'class': cls, 'secretome': secretome, 'main': 'plasma membrane', 'rel_if': rel}


def test_hygiene_flags_immunoglobulin_genes():
    assert HY.classify(_rec('immunoglobulin genes, predicted secreted proteins'))['hard_exclude']


def test_hygiene_keeps_membrane_proteins_even_if_also_labelled_intracellular():
    """The v1 rule failed exactly here (MET, MSLN, DLL3 ...). Positive membrane evidence wins."""
    c = HY.classify(_rec('predicted intracellular proteins, predicted membrane proteins'))
    assert not c['hard_exclude'] and c['tier'] == 'A'


def test_hygiene_excludes_intracellular_only_genes():
    assert HY.classify(_rec('predicted intracellular proteins'))['hard_exclude']


def test_hygiene_does_not_exclude_cd_markers():
    assert not HY.classify(_rec('cd markers, predicted intracellular proteins'))['hard_exclude']


def test_secretome_is_a_soft_flag_only():
    c = HY.classify(_rec('predicted membrane proteins', secretome='Secreted to blood'))
    assert not c['hard_exclude'] and c['secretome_annotated']


def test_hygiene_false_exclusion_rate_is_bounded():
    fx = J('candidate_hygiene.json')['false_exclusion_check']
    assert fx['false_exclusion_rate_upper_bound'] <= 0.10     # measured 6/107
    # the ones it removes must not include canonical surface antigens
    assert not {'MET', 'MSLN', 'DLL3', 'FOLH1', 'ROR1', 'PTK7', 'ERBB3'} & set(fx['fails'])


def test_positive_control_recoveries_survive_hygiene():
    assert J('candidate_hygiene.json')['positive_control']['fail_hygiene'] == {}


def test_corrected_candidate_funnel():
    h = J('candidate_hygiene.json')
    assert h['n_candidates'] == 214 and h['n_pass'] == 134
    assert h['n_pass_and_not_already_pursued'] == 128


# ------------------------------------------------------------- screen rank stability
def test_immunoglobulin_genes_leave_the_ranked_list():
    genes, _, _ = SR._load(hygiene=True)
    assert not [g for g in genes if g.startswith('IGKV')]
    all_genes, _, _ = SR._load(hygiene=False)
    assert any(g.startswith('IGKV') for g in all_genes)


def test_broad_gene_artefacts_stay_demoted_at_every_gate():
    for g, v in J('screen_rank_stability.json')['gate_sweep'].items():
        assert v['artefact_ranks']['AKAP10'] >= 200 and v['artefact_ranks']['TGFBR1'] >= 200


def test_shortlist_stability_is_limited_and_reported():
    r = J('screen_rank_stability.json')
    assert {'ERVFRD-1', 'DCT'} <= set(r['stable_and_not_already_pursued'])
    assert r['n_top25_in_>=90pct_of_draws'] < 15      # pinned: most of the top 25 is NOT robust


# ------------------------------------------------------------- Angle 2 diagnosis
def test_myeloma_clinical_units_sit_at_noise_floor():
    """Pinned negative: at 1% harm the panel's FITC/PE control channels match or beat BCMA and CD38."""
    r = J('clinical_units_noise_floor.json')
    s = r['summary']['0.01']
    assert s['best_control_ceiling'] >= r['per_marker']['CD269']['0.01']['ceiling']
    assert s['best_control_ceiling'] >= r['per_marker']['CD38']['0.01']['ceiling']
    assert all(v['n_above_floor'] <= 6 for v in r['summary'].values())


def test_dlbcl_populations_reproduce_published_auc():
    x = J('clinical_units_dlbcl.json')['crosscheck_vs_published']
    assert max(v['max_abs_diff'] for v in x.values()) < 0.005


def test_dlbcl_clean_compartment_gives_large_kill_and_contaminated_one_does_not():
    d = J('clinical_units_dlbcl.json')['per_marker']['CD19']['0.01']['per_compartment']
    assert d['T'] > 0.8 and d['O'] < 0.05


def test_decontamination_raises_other_leukocyte_kill():
    m = J('clinical_units_contamination.json')['main']['CD19']
    assert m['kill_after']['O'][0] > 10 * m['kill_before']['O'][0]
    assert m['contamination_frac_median']['O'] > 0.10


# ------------------------------------------------------------- Angle 5
def test_pav_is_monotone_and_exact_on_monotone_data():
    x = np.arange(10.0)
    ux, fit = A5.pav(x, x ** 2)
    assert np.all(np.diff(fit) >= 0) and np.allclose(fit, x ** 2)
    ux, fit = A5.pav([1, 2, 3, 4], [1, 3, 2, 4])          # violator pooled
    assert np.all(np.diff(fit) >= 0) and fit[1] == fit[2] == 2.5


def test_angle5_precommitted_wt_tests_did_not_support_premise():
    """Pinned honest result: WT->variant on-target R2 is far below the pre-stated 0.40 line."""
    t1 = J('angle5_monotone_premise.json')['T1_on_target_one_map']
    assert t1['wt->hf']['precommitted'] and t1['wt->hf']['cv_r2'] < 0.40
    assert t1['wt->lz']['precommitted'] and t1['wt->lz']['cv_r2'] < 0.40


def test_angle5_posthoc_result_is_labelled_posthoc():
    r = J('angle5_monotone_premise.json')
    assert r['T1_on_target_one_map']['hf->lz']['precommitted'] is False
    assert r['T1_on_target_one_map']['hf->lz']['cv_r2'] > 0.6
    assert abs(r['T3_on_vs_off_same_curve']['hf->lz']['median_signed_residual']) < 0.02
