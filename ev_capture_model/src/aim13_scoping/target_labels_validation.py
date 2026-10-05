"""FROZEN LABEL TABLE for the VALIDATION set: targets found after the main panel, labelled BEFORE any expression score was computed.

Same rules as target_labels.py:  safe = approval or documented clinical tolerability, no dose-limiting on-target normal-tissue toxicity;
toxic = documented severe / dose-limiting / fatal toxicity attributed (by the authors) to the target on a normal tissue.
Found by two independent literature searches (Sept 29, 2026) that were not shown any expression data; every source below was opened
by the searcher (PubMed abstract or full text, or the FDA label). Confidence and hedging are recorded honestly, because almost
every clinical on-target attribution is hedged by the authors (payload, lymphodepletion, CRS, or off-target cross-reactivity).

tier  'primary'   = pre-specified for the validation test (confidence medium or higher, not flagged ambiguous)
      'hedged'    = attribution or severity is disputed; reported separately and descriptively, never pooled into the primary test
      'excluded'  = ambiguous under the frozen rules (same gene has an approval AND a toxicity record), non-target, or out of scope
Columns: gene, common, agent, modality, outcome, tier, confidence, evidence, source
"""
VALIDATION = [
 ('MAGEA4','MAGE-A4','afamitresgene autoleucel (Tecelra)','TCR-T','safe','primary','high','FDA accelerated approval 2024; adverse events are lymphodepletion cytopenias and CRS, none attributed to normal-tissue on-target effects; no treatment-related deaths in SPEARHEAD-1','PMID 38554725; DailyMed ab24631f-3364-46e1-8074-7244863bcbab'),
 ('SSTR2','somatostatin receptor 2','177Lu-dotatate (Lutathera)','radioligand','safe','primary','medium-high','FDA 2018/EMA 2017; NETTER-1 no evidence of renal toxicity, grade 3-4 neutropenia 1 percent; kidney dose is tubular reabsorption, not receptor binding','PMID 28076709, 34793718; DailyMed 72d1a024-00b7-418a-b36e-b2cb48f2ab55'),
 ('CCR4','CCR4','mogamulizumab (Poteligeo)','depleting mAb','safe','primary','medium','FDA 2018; Treg/T-cell depletion, rash 25-35 percent manageable, no boxed warning','PMID 30100375; DailyMed e53960ab-42a1-40d1-9c7d-eb013fe7f18f'),
 ('CTAG1B','NY-ESO-1','NY-ESO-1 TCR-T incl. letetresgene autoleucel','TCR-T','safe','primary','medium','Well tolerated across myeloma, synovial sarcoma and MRCLS series; adverse events are CRS and cytopenias; not approved (breakthrough designation)','PMID 26193344, 29891538, 39836945'),
 ('IL2RA','CD25','denileukin diftitox-cxdl (Lymphir)','immunotoxin','safe','primary','medium','FDA 2024 (CTCL); Treg/activated-T depletion; hallmark toxicity is capillary leak (boxed), a toxin-class effect not attributed to the target','PMID 39700456; DailyMed ac3613d7-4304-40e7-a3f5-adb8a8adca5d'),
 ('GUCY2C','guanylyl cyclase C','IM96 CAR-T; PF-07062119 (GUCY2C x CD3)','CAR-T; TCE','toxic','primary','medium-high','CAR-T: 55 percent grade 3 diarrhoea, authors call it a notable on-target off-tumour toxicity; TCE: dose-limiting colitis/diarrhoea, plus CRS','PMID 42241671 (full text); PMID 42696650 (abstract)'),
 ('PSCA','prostate stem cell antigen','PSCA CAR-T (Dorff 2024)','CAR-T','toxic','primary','medium','Two dose-limiting grade 3 cystitis events needing steroids, authors: likely on-target off-tumour with a contribution from cyclophosphamide','PMID 38867077 (full text)'),
 ('TFRC','transferrin receptor (CD71)','CX-2029','ADC (MMAE)','toxic','hedged','low-medium','Dose-limiting anaemia/neutropenia; authors: etiology likely multifactorial; MMAE overlap','PMID 34083236'),
 ('SLC44A4','SLC44A4','ASG-5ME','ADC (MMAE)','toxic','hedged','low-medium','Two grade 5 events, narrow therapeutic index; authors: both on-target and off-target effects','PMID 30725389'),
 ('MLANA','MART-1','anti-MART-1 TCR-T (Johnson 2009)','TCR-T','toxic','hedged','high attribution, low severity','Melanocyte destruction in skin, eye, ear (uveitis, hearing loss) sometimes needing steroids; no grade/DLT data seen','PMID 19451549 (abstract)'),
 ('MAGEA12','MAGE-A12','anti-MAGE-A3 TCR (Morgan 2013)','TCR-T','toxic','hedged','low-medium','Two deaths with brain necrosis; authors suggest unrecognised MAGE-A12 in brain; off-target cannot be excluded; label gene disputed','PMID 23377668 (full text)'),
 ('EPHA2','EphA2','MEDI-547','ADC (auristatin)','toxic','hedged','low','Bleeding/coagulation events, trial stopped after 6 patients; attribution unresolved','PMID 22370972'),
 ('PMEL','gp100','tebentafusp (approved); gp100 TCR-T (Johnson 2009)','ImmTAC; TCR-T','excluded','excluded','ambiguous','Same gene has an approval AND an on-target melanocyte toxicity record, ambiguous by the frozen rule (as for IL3RA)','PMID 34551229, 19451549'),
 ('SLC6A2','norepinephrine transporter','iobenguane I-131 (Azedra)','radioligand','excluded','excluded','low-medium','Transporter substrate not an antigen; product withdrawn 2023; myelosuppression is a radionuclide effect','PMID 30291194'),
 ('CD276','B7-H3','ifinatamab deruxtecan','ADC','excluded','excluded','low','Regulatory status unresolved; early-phase data only','PMID 41926962'),
 ('PRAME','PRAME','IMA203','TCR-T','excluded','excluded','low','Phase 1 only','PMID 40205198'),
]
PRIMARY_SAFE = [t[0] for t in VALIDATION if t[5] == 'primary' and t[4] == 'safe']
PRIMARY_TOXIC = [t[0] for t in VALIDATION if t[5] == 'primary' and t[4] == 'toxic']
HEDGED_TOXIC = [t[0] for t in VALIDATION if t[5] == 'hedged' and t[4] == 'toxic']
