import pathlib,sys,unittest,math
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]))
from thresholds import select_threshold, fraction_above, auc, validate_count

class ThresholdTests(unittest.TestCase):
    def test_direction_ties_and_invalid_auc(self):
        self.assertEqual(auc([1],[2]),0)
        self.assertEqual(auc([2],[1]),1)
        self.assertEqual(auc([1],[1]),.5)
        with self.assertRaises(ValueError): auc([], [1])
    def test_known_thresholds(self):
        r=select_threshold([[1,2,3]],[[0,1,2]],1/3,'equal_donor')
        self.assertEqual(r['threshold'],2)
        self.assertAlmostEqual(r['normal_mean'],1/3)
        self.assertAlmostEqual(r['tumor_mean'],2/3)
    def test_zero_budget_and_none(self):
        self.assertEqual(select_threshold([[1]],[[1]],0,'donor_envelope')['tumor_worst'],0)
        self.assertEqual(select_threshold([[1]],[[0]],0,'donor_envelope')['tumor_worst'],1)
        self.assertEqual(fraction_above([1,2],math.inf),0)
    def test_donor_weights_and_envelope(self):
        # One large low-expression donor must not swamp a high-expression donor.
        normal=[[0]*99,[2]]; tumor=[[1,3]]
        a=select_threshold(tumor,normal,.1,'pooled')
        b=select_threshold(tumor,normal,.1,'donor_envelope')
        self.assertEqual(a['tumor_mean'],1)
        self.assertEqual(a['normal_worst'],1)
        self.assertEqual(b['tumor_worst'],.5)
        self.assertEqual(b['normal_worst'],0)
    def test_exhaustive_small_envelope(self):
        import itertools
        for a,b,c,d in itertools.product(range(3),repeat=4):
            tumor=[[a,b]]; normal=[[c],[d]]
            r=select_threshold(tumor,normal,.5,'donor_envelope')
            feasible=[fraction_above(tumor[0],t) for t in [0,1,2,3] if max(fraction_above(n,t) for n in normal)<=.5]
            self.assertEqual(r['tumor_worst'],max(feasible))
    def test_invalid_measurements(self):
        for x in [-1,float('nan')]:
            with self.assertRaises(ValueError): select_threshold([[x]],[[0]],.1,'pooled')
        with self.assertRaises(ValueError):select_threshold([[1]],[],.1,'pooled')
        with self.assertRaises(ValueError):select_threshold([[1]],[[0]],1.1,'pooled')
    def test_count_denominators(self):
        with self.assertRaises(ValueError):validate_count(-1,4,'patients')
        with self.assertRaises(ValueError):validate_count(5,4,'patients')
        with self.assertRaises(ValueError):validate_count('NR',4,'patients')
        self.assertEqual(validate_count(5,4,'events'),5)

if __name__=='__main__': unittest.main()
