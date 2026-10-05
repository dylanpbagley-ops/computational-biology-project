"""Tests for continuous expendability weights, the sensitivity analysis, and the re-ranked
discovery screen.

The most important assertion in this file is `test_continuous_is_worse_than_binary`. The
continuous scheme performs WORSE than the binary flag on the clinical prediction task, and
that negative result is pinned deliberately so a future session cannot quietly adopt the
continuous scheme as primary without noticing.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'aim13_scoping'))

import expendability_weights as W  # noqa: E402
import continuous_weight_sensitivity as S  # noqa: E402
import screen_rerank_weighted as R  # noqa: E402


# ---- the weight table itself --------------------------------------------------------
def test_weights_are_in_range():
    for ct, (w, rationale) in W.WEIGHTS.items():
        assert 0.0 <= w <= 1.0, f'{ct}: weight {w} out of range'
        assert rationale.strip(), f'{ct}: weight with no stated rationale'


def test_binary_scheme_reproduces_the_section_23_rule():
    """Under 'binary' the score must be exactly top_share + 1.0 for a haematopoietic home."""
    assert W.score(0.3, 'b-cells', 'binary') == pytest.approx(1.3)
    assert W.score(0.3, 'cardiomyocytes', 'binary') == pytest.approx(0.3)


def test_blood_lineages_are_not_treated_as_uniformly_expendable():
    """The substantive correction: T cells and microglia are not B cells."""
    assert W.weight('b-cells', 'continuous') < 0.1
    assert W.weight('t-cells', 'continuous') > 0.5
    assert W.weight('microglia', 'continuous') > 0.8
    assert W.weight('neutrophils', 'continuous') > 0.4


def test_regenerating_epithelium_is_not_treated_as_essential():
    """The other correction: gastric mucosa is not a cardiomyocyte."""
    assert W.weight('foveolar cells', 'continuous') < 0.5
    assert W.weight('melanocytes', 'continuous') < 0.3
    assert W.weight('cardiomyocytes', 'continuous') == pytest.approx(1.0)


def test_colonocytes_are_penalised_despite_regenerating():
    """CEACAM5 colitis: barrier failure during the regeneration gap is the constraint."""
    assert W.weight('colonocytes', 'continuous') >= 0.4


def test_unlisted_cell_type_defaults_to_non_expendable():
    assert W.weight('a cell type that does not exist', 'continuous') == W.DEFAULT_WEIGHT
    assert W.DEFAULT_WEIGHT >= 0.5, 'default must be conservative, not permissive'


# ---- the sensitivity analysis -------------------------------------------------------
@pytest.fixture(scope='module')
def sens():
    return S.compute_all()


def test_binary_baseline_reproduces_documented_value(sens):
    assert sens['binary']['sep'] == pytest.approx(0.9139, abs=1e-3)


def test_continuous_is_worse_than_binary(sens):
    """PINNED NEGATIVE RESULT.

    Continuous weighting degrades the clinical prediction (0.809 vs 0.914). The cause is
    identified: continuous weighting trusts the computed home cell type, and for broadly
    expressed targets that home is an artifact -- EGFR's computed home is placental tissue,
    which this scheme scores as fully expendable, so a toxic target collects a large bonus.
    The binary flag is robust to this precisely because it is coarse.

    If a future change makes continuous beat binary, that is a real improvement and this
    test SHOULD fail -- but it must fail loudly and be re-examined, not slip through.
    """
    assert sens['continuous']['sep'] < sens['binary']['sep']
    assert sens['continuous']['sep'] == pytest.approx(0.809, abs=0.01)
    assert sens['_summary']['continuous_minus_binary'] < 0


def test_every_real_scheme_remains_significant(sens):
    assert sens['_summary']['p_max'] < 0.01


def test_shuffled_null_is_reported_and_is_not_beaten_decisively(sens):
    """Honest record: the continuous scheme does not clear its own scrambled-weight null."""
    n = sens['shuffled_null']
    assert n['n_shuffles'] >= 100
    assert n['p_vs_null'] > 0.05, 'if this now clears the null, re-examine and report it'


# ---- the re-ranked screen -----------------------------------------------------------
@pytest.fixture(scope='module')
def rerank():
    return R.compute_all()


def test_broadly_expressed_artifacts_are_demoted(rerank):
    by = rerank['by_gene']
    for g in ('AKAP10', 'TGFBR1', 'SIRPA'):
        assert by[g]['rank_weighted'] > by[g]['rank_binary'], f'{g} was not demoted'
    assert by['AKAP10']['rank_weighted'] > 200


def test_tissue_restricted_candidates_are_promoted(rerank):
    by = rerank['by_gene']
    for g in ('DCT', 'ERVFRD-1', 'PSG1'):
        assert by[g]['rank_weighted'] < 10, f'{g} did not reach the top 10'
        assert by[g]['rank_weighted'] < by[g]['rank_binary']


def test_t_cell_gene_no_longer_sits_in_the_top_ten(rerank):
    """CD3D in the top 10 was the most obviously indefensible feature of the old ranking."""
    assert rerank['by_gene']['CD3D']['rank_weighted'] > 10


def test_top_candidates_are_actually_confined(rerank):
    """Nothing in the new top 10 should be a barely-expressed-anywhere gene."""
    for c in rerank['top_25_weighted'][:10]:
        assert c['confinement'] >= rerank['gate'], \
            f"{c['gene']} is in the top 10 with confinement {c['confinement']:.3f}"


def test_gate_is_documented_as_a_judgement_not_an_optimum():
    assert 'judgement' in R.__doc__.lower() or 'judgment' in R.__doc__.lower()
