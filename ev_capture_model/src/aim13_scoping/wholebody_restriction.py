"""
FEASIBILITY TEST for extending the project beyond blood cancers.

The outcome test (section 22) showed that the measure which actually predicts whether a
drug succeeded is RESTRICTION to the cancer's own normal lineage, not the tumor-vs-normal
ceiling. That measure only needs to know how confined a protein is across normal cell types
-- which means it can be computed for ANY cancer, not just ones whose normal counterpart
sits in blood, as long as a whole-body normal reference exists.

Reference: Human Protein Atlas single-cell consensus, 154 human cell types spanning every
major organ (cardiomyocytes, alveolar cells, hepatocytes, pancreatic cells, neurons, ...),
~20,000 genes, 3,087,080 rows.

Measure (aggregate per cell type, so this is a restriction score, not a per-cell AUC):
   top_share  = expression in the single highest cell type / total across all 154
   n_expr     = number of cell types above 10% of that maximum ("how many places is it on?")
   home       = the highest-expressing cell type
A target is 'lineage-confined' when top_share is high and n_expr is small.

BLOOD POSITIVE CONTROLS are included so the measure can be checked against the result it
must reproduce: the approved B-cell targets should come out confined to B cells.
"""
import os, io, json, zipfile, collections, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
D=os.path.abspath(os.path.join(HERE,'..','..','..','ev_capture_data','normal_tissue_atlas'))
Z=os.path.join(D,'rna_single_cell_type.tsv.zip')

GENES={
 # --- blood positive controls: the measure must reproduce the known answer ---
 'CD19':'B approved','MS4A1':'B approved (CD20)','CD22':'B approved','CD79B':'B approved',
 'TNFRSF13C':'B, BAFF-R, trials','CD72':'B, trials','FCER2':'B FAILED (CD23/lumiliximab)',
 'CD40':'B FAILED (dacetuzumab)','CD80':'B FAILED (galiximab)','HLA-DRA':'B FAILED (apolizumab)',
 'PTPRC':'pan-leukocyte, not targetable','CD44':'broad, never worked',
 'IL3RA':'AML FAILED (talacotuzumab)','CD33':'AML approved (gemtuzumab)','CD47':'AML FAILED (magrolimab)',
 # --- myeloma ---
 'TNFRSF17':'myeloma approved (BCMA)','SDC1':'myeloma FAILED (CD138/indatuximab)',
 'SLAMF7':'myeloma approved (elotuzumab)','CD38':'myeloma approved (daratumumab)',
 # --- SOLID TUMOUR targets, the actual test ---
 'ERBB2':'HER2 - CAR-T caused a FATAL on-target event (Morgan 2010)',
 'EGFR':'EGFR - CAR-T skin/lung toxicity',
 'CEACAM5':'CEA - CAR-T severe colitis',
 'CA9':'CAIX - CAR-T liver toxicity',
 'MSLN':'mesothelin - CAR-T tolerable; a PANCREATIC target',
 'CLDN18':'claudin-18.2 - zolbetuximab APPROVED 2024; CAR-T in gastric/PANCREATIC',
 'MUC1':'MUC1 - PANCREATIC target, no approval',
 'GPC3':'glypican-3 - liver cancer, tolerable',
 'FOLH1':'PSMA - approved radioligand',
 'TACSTD2':'TROP2 - sacituzumab APPROVED',
 'NECTIN4':'Nectin-4 - enfortumab APPROVED',
 'CD276':'B7-H3 - trials',
 'IL13RA2':'IL13Ra2 - glioma trials',
 'ROR1':'ROR1 - trials',
 'DLL3':'DLL3 - tarlatamab APPROVED (small cell lung)',
}

z=zipfile.ZipFile(Z); name=z.namelist()[0]
vals=collections.defaultdict(dict)
with z.open(name) as f:
    fh=io.TextIOWrapper(f,'utf-8'); hdr=fh.readline().rstrip('\n').split('\t')
    gi,ci,vi=hdr.index('Gene name'),hdr.index('Cell type'),hdr.index('nCPM')
    for line in fh:
        p=line.rstrip('\n').split('\t')
        if p[gi] in GENES:
            try: vals[p[gi]][p[ci]]=float(p[vi])
            except ValueError: pass
print('cell types in reference: %d ; genes found: %d of %d\n' % (
    len(next(iter(vals.values()))), len(vals), len(GENES)))

rows=[]
for g,note in GENES.items():
    d=vals.get(g)
    if not d: print('%-10s NOT FOUND' % g); continue
    v=np.array(list(d.values())); ks=list(d.keys())
    tot=v.sum()
    if tot<=0: print('%-10s no expression' % g); continue
    top=int(np.argmax(v))
    share=float(v[top]/tot); nexpr=int((v>0.10*v[top]).sum())
    rows.append((g,share,nexpr,ks[top],note))
rows.sort(key=lambda r:-r[1])
print('%-10s %10s %7s  %-32s %s' % ('GENE','top_share','n_types','highest cell type','clinical note'))
for g,s,n,h,note in rows:
    print('%-10s %10.3f %7d  %-32s %s' % (g,s,n,h[:32],note))
json.dump({g:{'top_share':s,'n_types':n,'home':h,'note':note} for g,s,n,h,note in rows},
          open(os.path.join(HERE,'wholebody_restriction.json'),'w'), indent=1)
