import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np

STAGE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STAGE / 'dependencies'))
sys.path.insert(0, str(STAGE / 'analysis'))
from load_paired_counts import (
    assert_patient_disjoint, joint_thin_counts,
    load_sample, resolve_target_features, summarize_pair, validate_csc_structure,
)
import h5py


class PairedAnalysisTests(unittest.TestCase):
    def test_joint_counts_and_directional_conditional(self):
        counts = np.array([[1,1]]*2 + [[1,0]] + [[0,1]]*3 + [[0,0]]*4, dtype=np.int32)
        got = summarize_pair(counts, 0, 1)
        self.assertEqual((got['n11_both'], got['n10_a_only'], got['n01_b_only'], got['n00_neither']), (2,1,3,4))
        self.assertAlmostEqual(got['or_coverage_fraction'], 0.6)
        self.assertAlmostEqual(got['b_given_a_undetected'], 3/7)
        self.assertAlmostEqual(got['a_given_b_undetected'], 1/5)

    def test_same_marginals_can_have_different_joint_profiles(self):
        first = np.array([[1,1]]*2 + [[1,0]] + [[0,1]]*3 + [[0,0]]*4, dtype=np.int32)
        second = np.array([[1,1]] + [[1,0]]*2 + [[0,1]]*4 + [[0,0]]*3, dtype=np.int32)
        one = summarize_pair(first, 0, 1)
        two = summarize_pair(second, 0, 1)
        self.assertEqual((one['a_detected_n'], one['b_detected_n']), (two['a_detected_n'], two['b_detected_n']))
        self.assertNotEqual(one['n11_both'], two['n11_both'])
        self.assertNotEqual(one['or_coverage_fraction'], two['or_coverage_fraction'])

    def test_zero_conditional_denominator_is_unavailable(self):
        counts = np.array([[1,1],[1,0],[1,0]], dtype=np.int64)
        got = summarize_pair(counts, 0, 1)
        self.assertIsNone(got['b_given_a_undetected'])
        self.assertAlmostEqual(got['a_given_b_undetected'], 1.0)

    def test_exact_target_symbols_resolve_under_feature_reordering(self):
        names = ['SLAMF7','GPRC5D','TNFRSF17']
        types = ['Gene Expression']*3
        self.assertEqual(resolve_target_features(names, types), [2,1,0])
        with self.assertRaises(ValueError):
            resolve_target_features(['SLAMF7','GPRC5D'], types[:2])
        with self.assertRaises(ValueError):
            resolve_target_features(['SLAMF7','GPRC5D','TNFRSF17','TNFRSF17'], types+['Gene Expression'])

    def test_csc_rejects_duplicate_row_coordinates(self):
        validate_csc_structure(np.array([0,2,3]), np.array([0,1,2]), 3, 2)
        with self.assertRaisesRegex(ValueError, 'duplicate feature-row'):
            validate_csc_structure(np.array([0,3]), np.array([1,1,2]), 3, 1)

    def test_h5_loader_uses_specimen_qualified_cells_and_feature_names(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / 'mini.h5'
            with h5py.File(path, 'w') as h5:
                m = h5.create_group('matrix')
                m.create_dataset('shape', data=np.array([3,2], dtype=np.int64))
                m.create_dataset('barcodes', data=np.array([b'cellB',b'cellA']))
                m.create_dataset('data', data=np.array([3,1,4,2], dtype=np.int32))
                m.create_dataset('indices', data=np.array([0,1,1,2], dtype=np.int32))
                m.create_dataset('indptr', data=np.array([0,2,4], dtype=np.int32))
                f = m.create_group('features')
                f.create_dataset('id', data=np.array([b'idS',b'idT',b'idG']))
                f.create_dataset('name', data=np.array([b'SLAMF7',b'TNFRSF17',b'GPRC5D']))
                f.create_dataset('feature_type', data=np.array([b'Gene Expression']*3))
            loaded = load_sample(path, 'GSM000001')
            self.assertEqual(loaded['cell_keys'], ['GSM000001|cellB','GSM000001|cellA'])
            self.assertEqual(loaded['target_counts'].tolist(), [[1,0,3],[4,2,0]])
            self.assertEqual(loaded['total_umi'].tolist(), [4,6])
            self.assertNotEqual('GSM000001|cellB', 'GSM000002|cellB')

    def test_joint_hypergeometric_preserves_joint_categories_and_seed(self):
        counts = np.array([[10,20,5],[0,0,0],[3,2,1]], dtype=np.int64)
        totals = np.array([100,10,6], dtype=np.int64)
        a = joint_thin_counts(counts, totals, 25, 9)
        b = joint_thin_counts(counts, totals, 25, 9)
        self.assertTrue(np.array_equal(a,b))
        self.assertLessEqual(int(a[0].sum()),25)
        self.assertEqual(int(a[1].sum()),0)
        self.assertEqual(a[2].tolist(), [3,2,1])
        self.assertTrue(np.all(a <= counts))

    def test_invalid_counts_and_patient_leakage_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'nonnegative integers'):
            summarize_pair(np.array([[1.5,0.0]]), 0, 1)
        assert_patient_disjoint(['P1','P2'], ['P3'])
        with self.assertRaisesRegex(ValueError, 'patient leakage'):
            assert_patient_disjoint(['P1','P2'], ['P2','P3'])


if __name__ == '__main__':
    unittest.main()



