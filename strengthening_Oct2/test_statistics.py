import unittest
import numpy as np
from diagnostics import auc, exact_p, paired_bootstrap
class StatisticsTests(unittest.TestCase):
    def test_auc_direction_and_ties(self):
        self.assertEqual(auc([2,3],[0,1]),1)
        self.assertEqual(auc([1,1],[1,1]),.5)
        self.assertEqual(auc([0],[1]),0)
    def test_exact_enumeration_including_ties(self):
        self.assertAlmostEqual(exact_p([2,3],[0,1]),1/6)
        self.assertEqual(exact_p([1,1],[1,1]),1)
        self.assertAlmostEqual(exact_p([1,2],[0,1]),2/6)
    def test_paired_not_independent_bootstrap(self):
        x=np.array([0.,.2,.6,1.]);y=np.array([0,0,1,1],bool)
        out=paired_bootstrap(x,x,y,200,7)
        self.assertEqual(out['difference'],0);self.assertEqual(out['ci95'],[0.,0.])
if __name__=='__main__':unittest.main()
