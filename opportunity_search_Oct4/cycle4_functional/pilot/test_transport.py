import unittest
import numpy as np
from transport_helpers import normalized_log_frequency,summarize_envelope,guide_filter_reasons,exact_source_rows,count_only_positions

class TransportHelpers(unittest.TestCase):
    def test_identical_frequencies_zero(self):
        counts=np.array([[100,100],[200,200],[300,300]],dtype=np.int64)
        logf=normalized_log_frequency(counts)
        np.testing.assert_allclose(logf[:,0]-logf[:,1],0,atol=1e-12)
    def test_depletion_negative_and_enrichment_positive(self):
        before=np.array([200,200],dtype=np.int64)
        depleted=np.array([50,350],dtype=np.int64)
        enriched=np.array([350,50],dtype=np.int64)
        d=normalized_log_frequency(np.column_stack([before,depleted]))[:,1]-normalized_log_frequency(np.column_stack([before,before]))[:,1]
        e=normalized_log_frequency(np.column_stack([before,enriched]))[:,1]-normalized_log_frequency(np.column_stack([before,before]))[:,1]
        self.assertLess(float(d[0]),0)
        self.assertGreater(float(e[0]),0)
    def test_column_scaling_invariant_up_to_pseudocount(self):
        base=np.array([[10000,20000],[30000,40000],[50000,60000]],dtype=np.int64)
        scaled=base.copy(); scaled[:,1]*=2
        a=normalized_log_frequency(base)[:,1]
        b=normalized_log_frequency(scaled)[:,1]
        np.testing.assert_allclose(a,b,atol=1e-5)
    def test_nine_value_envelope(self):
        med,lo,hi,cls=summarize_envelope(np.array([-3,-2,-1,-4,-2,-1,-3,-2,-1],float))
        self.assertEqual((med,lo,hi,cls),(-2.0,-4.0,-1.0,'sign_stable_sensitizing'))
        self.assertEqual(summarize_envelope(np.array([-2,-1,0,-1,-2,-1,-3,-1,-2],float))[3],'mixed')
    def test_low_missing_duplicate_keys_never_pass(self):
        self.assertTrue(guide_filter_reasons('g1','A','A',19,[25,25,25]))
        self.assertIn('duplicate_guide_id',guide_filter_reasons('g1','A','A',30,[25,25,25],duplicate_id=True))
        self.assertIn('missing_T0_count',guide_filter_reasons('g1','A','A',None,[25,25,25],missing_t0=True))
        self.assertIn('missing_primary_control_count',guide_filter_reasons('g1','A','A',30,[25,25],missing_control=True))
    def test_annotation_prefixed_header_uses_count_only_indices(self):
        full=['Guide.ID','Chr','Start','End','Strand','Gene.symbol','Guide.Seq','MUS_15','MUS_25']
        names,index=count_only_positions(full)
        self.assertEqual(names,['MUS_15','MUS_25'])
        self.assertEqual(index['MUS_15'],0)
        self.assertEqual(index['MUS_25'],1)
        counts=np.array([[11,23]],dtype=np.int64)
        self.assertEqual(int(counts[0,index['MUS_15']]),11)
        self.assertEqual(int(counts[0,index['MUS_25']]),23)
    def test_count_symbol_remains_authoritative_over_date_like_library_label(self):
        reasons=guide_filter_reasons('g1','March1','01-Mar',30,[25,25,25])
        self.assertEqual(reasons,[])
    def test_exact_source_keeps_renca_ha_separate(self):
        rows=[{'source':'Renca cell line + haemagglutinin antigen','treatment':'Control'},
              {'source':'Renca cell line','treatment':'IFN-gamma'}]
        picked=exact_source_rows(rows,'Renca cell line + haemagglutinin antigen')
        self.assertEqual(len(picked),1)
        self.assertEqual(picked[0]['treatment'],'Control')

if __name__=='__main__': unittest.main()
