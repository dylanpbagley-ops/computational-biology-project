import unittest
import numpy as np
from frontiers import upper_at, weighted_reference, tau_from_share

class FrontierTests(unittest.TestCase):
    def test_distinct_response_classes(self):
        p=[.45,.10,.45]; q=[.10,.80,.10]
        self.assertAlmostEqual(upper_at(p,q,.2,'arbitrary'),.9)
        self.assertAlmostEqual(upper_at(p,q,.2,'monotone'),23/45)
    def test_equal_lr_auc_different_frontier(self):
        self.assertAlmostEqual(upper_at([.6,.3,.1],[1/3]*3,1/3,'arbitrary'),.6)
        self.assertAlmostEqual(upper_at([.55,.4,.05],[1/3]*3,1/3,'arbitrary'),.55)
    def test_zero_normal_mass_and_ties(self):
        self.assertAlmostEqual(upper_at([.5,.5],[0,1],0,'arbitrary'),.5)
        self.assertAlmostEqual(upper_at([.5,.5],[0,1],0,'monotone'),0)
        self.assertAlmostEqual(upper_at([.2,.8],[.2,.8],.3,'arbitrary'),.3)
    def test_mixture_and_zero_weight(self):
        p,w=weighted_reference([[1,0],[0,1]],[.5,.5],[1,0])
        np.testing.assert_allclose(p,[1,0]);self.assertEqual(w,.5)
        p,w=weighted_reference([[1,0],[0,1]],[.5,.5],[0,0])
        self.assertIsNone(p);self.assertEqual(w,0)
    def test_invalid_input(self):
        for p,q,b in [([1,1],[.5,.5],.2),([-.1,1.1],[.5,.5],.2),([.5,.5],[.5,.5],1.1),([np.nan,1],[.5,.5],.2)]:
            with self.assertRaises(ValueError):upper_at(p,q,b,'arbitrary')
        with self.assertRaises(ValueError):upper_at([1],[1],.2,'unknown')
    def test_tau_identity(self):
        x=np.array([1.,2.,7.]);s=x.max()/x.sum()
        self.assertAlmostEqual(tau_from_share(s,len(x)),sum(1-x/x.max())/(len(x)-1))
    def test_random_ordered_rules_and_arbitrary_vertex_oracle(self):
        rng=np.random.default_rng(123)
        for _ in range(30):
            p=rng.dirichlet(np.ones(5));q=rng.dirichlet(np.ones(5));f=np.sort(rng.random(5));b=float(q@f)
            self.assertLessEqual(float(p@f),upper_at(p,q,b,'monotone')+1e-12)
            self.assertLessEqual(upper_at(p,q,b,'monotone'),upper_at(p,q,b,'arbitrary')+1e-12)
            # Independent LP vertex enumeration: all coordinates binary except at most one.
            oracle=0.
            import itertools
            for fractional in range(5):
                others=[i for i in range(5) if i!=fractional]
                for bits in itertools.product([0.,1.],repeat=4):
                    g=np.zeros(5);g[others]=bits;residual=b-q@g
                    if residual>=-1e-12:
                        g[fractional]=min(1.,max(0.,residual/q[fractional]));oracle=max(oracle,float(p@g))
            self.assertAlmostEqual(upper_at(p,q,b,'arbitrary'),oracle)

if __name__=='__main__':unittest.main()
