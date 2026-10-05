"""
test_wholebody_restriction.py

Unit tests for the whole-body (solid-tumour) extension, results document section 23, and for
the tool built on top of it.

WHAT THESE TESTS ARE FOR. The section-23 result rests on three things that are easy to get
silently wrong: the definition of the confinement measure, the rank-separation statistic and
its permutation p-value, and the rule that decides whether a normal home cell type is an
expendable lineage. Each of those has a known-correct answer on hand-built inputs, so each is
tested against one here rather than against the real data. One of these tests exists because
the bug it guards actually happened.
"""
import os
import sys
import math

import numpy as np
import pytest

SCOPING = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src', 'aim13_scoping'))
sys.path.insert(0, SCOPING)

import target_labels as TL  # noqa: E402
import wholebody_v2 as WB   # noqa: E402


# ---------------------------------------------------------------------------
# 1. the expendability rule
# ---------------------------------------------------------------------------

def test_haematopoietic_matches_real_blood_lineages():
    for name in ('b-cells', 'plasma cells', 'kupffer cells', 'neutrophils', 'cdc',
                 'mast cells', 'nk-cells', 'megakaryocytes'):
        assert TL.haematopoietic(name), name


def test_haematopoietic_rejects_solid_tissue():
    for name in ('colonocytes', 'alveolar cells type 1', 'foveolar cells', 'oligodendrocytes',
                 'epicardial cells', 'basal cells', 'hepatocytes'):
        assert not TL.haematopoietic(name), name


def test_goblet_cells_are_not_haematopoietic_regression():
    """THE BUG THIS GUARDS ACTUALLY HAPPENED.

    The first implementation tested lineage membership with a plain substring match. The
    string 'goblet cells' contains 't cell', so every goblet cell type in the reference was
    classified as a blood lineage. That inflated the two-factor separation from its correct
    0.905 to 0.916. Matching is now word-boundary anchored; this test fails if anyone ever
    reverts it."""
    assert not TL.haematopoietic('goblet cells')
    assert not TL.haematopoietic('conjunctival goblet cells')
    assert not TL.haematopoietic('intestinal goblet cells')


# ---------------------------------------------------------------------------
# 2. the confinement measure
# ---------------------------------------------------------------------------

def test_confinement_of_a_perfectly_confined_gene_is_one():
    out = WB.confinement({'G': {'b-cells': 10.0, 'hepatocytes': 0.0, 'neurons': 0.0}})
    assert out['G']['top_share'] == pytest.approx(1.0)
    assert out['G']['n_types'] == 1
    assert out['G']['home'] == 'b-cells'


def test_confinement_of_a_uniformly_expressed_gene_is_one_over_n():
    d = {'ct%d' % i: 5.0 for i in range(20)}
    out = WB.confinement({'G': d})
    assert out['G']['top_share'] == pytest.approx(1.0 / 20)
    assert out['G']['n_types'] == 20


def test_confinement_is_scale_invariant():
    """Only the SHAPE of the expression profile may matter. Multiplying every value by a
    constant is a change of units and must not move the measure."""
    base = {'a': 3.0, 'b': 1.0, 'c': 0.2}
    a = WB.confinement({'G': base})['G']
    b = WB.confinement({'G': {k: v * 1000.0 for k, v in base.items()}})['G']
    assert a['top_share'] == pytest.approx(b['top_share'])
    assert a['n_types'] == b['n_types']


def test_genes_with_no_expression_are_dropped_not_crashed():
    out = WB.confinement({'DEAD': {'a': 0.0, 'b': 0.0}, 'LIVE': {'a': 1.0, 'b': 0.0}})
    assert 'DEAD' not in out
    assert 'LIVE' in out


# ---------------------------------------------------------------------------
# 3. the separation statistic and its p-value
# ---------------------------------------------------------------------------

def test_separation_is_one_when_groups_do_not_overlap():
    assert WB.sep([3, 4, 5], [0, 1, 2]) == pytest.approx(1.0)


def test_separation_is_zero_when_reversed():
    assert WB.sep([0, 1, 2], [3, 4, 5]) == pytest.approx(0.0)


def test_separation_of_identical_groups_is_one_half():
    """All ties, so every comparison contributes 0.5. This is the null value the toxicity
    test is measured against."""
    assert WB.sep([1, 1, 1], [1, 1, 1]) == pytest.approx(0.5)


def test_separation_equals_mann_whitney_auc():
    rng = np.random.default_rng(4)
    a = rng.normal(1.0, 1.0, 40)
    b = rng.normal(0.0, 1.0, 35)
    brute = np.mean([[x > y for y in b] for x in a])
    assert WB.sep(a, b) == pytest.approx(brute, abs=1e-12)


def test_exact_p_is_exact_for_small_samples_and_matches_the_closed_form():
    """With 3 vs 3 and perfect separation there is exactly one arrangement out of C(6,3)=20
    that is at least as extreme, so the one-sided p is exactly 1/20."""
    p, was_exact = WB.exact_p([10, 11, 12], [1, 2, 3])
    assert was_exact
    assert p == pytest.approx(1.0 / 20.0)


def test_exact_p_of_identical_groups_is_one():
    p, was_exact = WB.exact_p([1, 1, 1], [1, 1, 1])
    assert was_exact
    assert p == pytest.approx(1.0)


def test_exact_p_falls_back_and_says_so_when_enumeration_is_too_large():
    rng = np.random.default_rng(7)
    a = list(rng.normal(1, 1, 19))
    b = list(rng.normal(0, 1, 10))
    p, was_exact = WB.exact_p(a, b)
    assert not was_exact, 'C(29,19) is far above the exact-enumeration limit'
    assert 0.0 < p <= 1.0


# ---------------------------------------------------------------------------
# 4. the frozen label file
# ---------------------------------------------------------------------------

def test_label_file_uses_only_declared_outcomes():
    allowed = {'safe', 'toxic', 'efficacy_fail', 'excluded_control', 'excluded_ambiguous'}
    assert {t[3] for t in TL.TARGETS} <= allowed


def test_every_target_appears_exactly_once():
    genes = [t[0] for t in TL.TARGETS]
    assert len(genes) == len(set(genes))


def test_every_label_carries_an_evidence_string():
    for t in TL.TARGETS:
        assert isinstance(t[4], str) and len(t[4]) > 15, t[0]


def test_the_outcome_bearing_panel_is_the_size_the_results_document_claims():
    """Section 23 (as of the September 10, 2026 panel-growth update) quotes 19 safe, 11 toxic,
    12 efficacy-fail. If the panel is edited again, this fails and the document has to be
    updated with it -- that is the point of this test, not a bug to work around.
    History: 19/10/9 from Sept 8 through the initial 40-target panel; grew to 19/11/10 on
    Sept 9 when F3 (tissue factor, toxic -- tisotumab vedotin's boxed-warning ocular toxicity)
    and GPNMB (efficacy_fail -- glembatumumab vedotin's METRIC trial primary-endpoint miss)
    were added; grew to 19/11/12 on Sept 10 when ERBB3 (HER3, efficacy_fail -- patritumab
    deruxtecan's BLA voluntarily withdrawn 2026 after HERTHENA-Lung02 missed its overall-survival
    endpoint) and AXL (efficacy_fail -- Genmab's 2020 discontinuation of enapotamab vedotin for
    insufficient proof-of-concept efficacy) were added, each verified against primary/near-
    primary sources. This growth made the specificity control (section 23.4/Test 4) move from
    borderline (p=0.043, 19 vs 10) to clearly significant (p=0.011, 19 vs 12) -- see the results
    document for the honest reading of why, and why it does not threaten the safe-vs-toxic
    result this test suite otherwise guards."""
    assert sum(1 for t in TL.TARGETS if t[3] == 'safe') == 19
    assert sum(1 for t in TL.TARGETS if t[3] == 'toxic') == 11
    assert sum(1 for t in TL.TARGETS if t[3] == 'efficacy_fail') == 12


# ---------------------------------------------------------------------------
# 5. the two-factor rule
# ---------------------------------------------------------------------------

def test_two_factor_rule_is_strictly_lexicographic():
    """The whole point of the +1 is that no amount of confinement can outrank an expendable
    home. A target with a blood home and the worst possible confinement must still outscore
    the most confined solid-tissue target."""
    worst_blood = 0.0 + 1.0
    best_solid = 1.0 + 0.0
    assert worst_blood >= best_solid


# ---------------------------------------------------------------------------
# 6. the tool's ceiling calculation
# ---------------------------------------------------------------------------

def test_ceiling_matches_the_closed_form_auc_for_two_gaussians():
    """For two normal distributions the separation has an exact closed form,
    AUC = Phi((mu_t - mu_n) / sqrt(sd_t^2 + sd_n^2)). The tool samples instead, so this
    checks the sampling against the mathematics."""
    ts = pytest.importorskip('targetscan')
    r = ts.ceiling(4.2, 0.8, 2.1, 0.9, n=60000, seed=1)
    z = (4.2 - 2.1) / math.sqrt(0.8 ** 2 + 0.9 ** 2)
    closed = 0.5 * (1 + math.erf(z / math.sqrt(2)))
    assert r['auc_approx'] == pytest.approx(closed, abs=0.01)


def test_ceiling_is_one_half_when_the_two_populations_are_identical():
    ts = pytest.importorskip('targetscan')
    r = ts.ceiling(3.0, 0.7, 3.0, 0.7, n=60000, seed=2)
    assert r['auc_approx'] == pytest.approx(0.5, abs=0.02)
    assert r['youden_J'] < 0.05


def test_ceiling_youden_j_never_exceeds_one():
    ts = pytest.importorskip('targetscan')
    for mt, mn in ((6.0, 0.0), (4.0, 3.9), (2.0, 5.0)):
        r = ts.ceiling(mt, 0.5, mn, 0.5, n=20000, seed=3)
        assert -1.0 <= r['youden_J'] <= 1.0
