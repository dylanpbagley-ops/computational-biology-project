import csv,json,pathlib,sys,unittest
OUT=pathlib.Path(__file__).resolve().parents[2]

def read(p):
    with p.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

class EvidenceTests(unittest.TestCase):
    def test_review_initial_judgments_blank(self):
        for filename in ['operational_index.tsv','batch_studies.tsv']:
            rows=read(OUT/'evidence/reviewer_packet'/filename)
            self.assertTrue(rows)
            for r in rows:
                self.assertFalse(any(x in r for x in ['score','historical_outcome','assistant_judgment']))
                self.assertTrue(all(v=='' for k,v in r.items() if k.startswith('initial_') or k in ['reviewer','review_date']))
    def test_missing_and_zero_are_distinct(self):
        rows=read(OUT/'evidence/observations.tsv')
        missing=next(r for r in rows if r['observation_id']=='OBS-CD19-DEPLETION')
        zero=next(r for r in rows if r['observation_id']=='OBS-CA9-NO-LIVER')
        self.assertEqual(missing['event_count'],'NR')
        self.assertEqual(missing['observation_state'],'indeterminate')
        self.assertEqual(zero['event_count'],'0')
        self.assertEqual(zero['denominator'],'4')
        self.assertEqual(zero['observation_state'],'explicit_nonobservation')
        self.assertFalse((OUT/'evidence/amended_analysis_labels.tsv').exists())
    def test_all_operational_records_preserved_in_screening(self):
        rows=read(OUT/'evidence/screening.tsv')
        self.assertEqual(len(rows),44);self.assertEqual(len({r['gene'] for r in rows}),44)
        fap=next(r for r in rows if r['gene']=='FAP')
        self.assertEqual(fap['evidence_state'],'ineligible')
    def test_trial_reports_not_independent_units(self):
        rs=read(OUT/'evidence/overlaps.tsv')
        r=next(r for r in rs if r['record_a']=='E003' and r['record_b']=='E004')
        self.assertIn('first three',r['relationship'])
        self.assertEqual(r['overlap_group'],'CAIX12')
    def test_every_held_out_sample_excluded(self):
        rows=read(OUT/'analysis/results/b_held_out.tsv')
        self.assertEqual(len(rows),3*3*8*9)
        for r in rows:
            t=set(r['train_tumor_ids'].split(','));n=set(r['train_normal_ids'].split(','))
            self.assertEqual(len(t),7);self.assertEqual(len(n),8)
            self.assertNotIn(r['held_tumor'],t)
            self.assertNotIn(r['held_normal'],n)
            self.assertFalse(t&n)
    def test_empirical_envelope_constraints(self):
        rows=read(OUT/'analysis/results/b_summary.tsv')
        for r in rows:
            if r['policy']=='donor_envelope':self.assertLessEqual(float(r['normal_worst']),float(r['budget'])+1e-12)
        # Counterexample to calling training-envelope protection universal.
        held=read(OUT/'analysis/results/b_held_out_summary.tsv')
        self.assertTrue(any(float(r['normal_worst'])>.1 for r in held if r['policy']=='donor_envelope'))
    def test_protected_inputs_unchanged(self):
        import hashlib
        hashes=json.loads((OUT/'input_hashes.json').read_text())
        for path,expected in hashes.items():
            self.assertEqual(hashlib.sha256((OUT.parent/path).read_bytes()).hexdigest(),expected,path)

if __name__=='__main__':unittest.main()
