"""
Load the healthy-donor reference: GSE164378 (Hao et al., Cell 2021).
161,764 PBMCs from 8 HEALTHY donors, 228 surface proteins per cell, with the authors'
own expert cell-type annotations (so the normal compartments are theirs, not mine).

This is the first genuinely healthy, multi-lineage human reference in the project, and it
covers the targets of BOTH diseases studied here (CD19/CD20/CD22/CD79b for lymphoma;
CD269/BCMA, CD38, CD319/SLAMF7, CD138, CD307/FcRH5 for myeloma), plus four isotype
controls and eleven antigens measured with two independent antibody clones.
"""
import os, sys, gzip, numpy as np, pandas as pd
HERE=os.path.dirname(os.path.abspath(__file__))
D=os.path.abspath(os.path.join(HERE,'..','..','..','ev_capture_data','healthy_pbmc_citeseq'))
NAMES = """CD39,Rat-IgG1-1,CD107a,CD62P,TCR-2,CD30,CD31,CD34,CD35,CD36,CD223,TIGIT,TCR-V-9,CD226,CD178,CD319,CD171,Siglec-8,CD340,Rat-IgG2b,VEGFR-3,CD29,CD62E,CD4-2,CD4-1,CD22,CD3-1,CD20,CD27,CD45RB,CD25,CD24,CD146,Galectin-9,CD142,CD141,CD294,Rat-IgG1-2,CD45RA,CX3CR1,CD56-2,CD56-1,CD45RO,CD303,GP130,CD253,CD357,CD11b-1,CD354,CD11b-2,CLEC12A,CD38-2,CD38-1,Folate,Rag-IgG2c,CD209,CD152,CD154,CD155,Cadherin,CD201,CD204,CD205,CD206,CD207,CD1d,CD284,CD1c,Podoplanin,CD1a,CD366,IgD,IgM,CD66a/c/e,CD49d,LOX-1,TIM-4,CD98,CD370,CD49a,CD44-2,C5L2,CD44-1,CD158e1,CD124,CD127,CD126,CD279,CD278,CD123,CD122,CD96,CD274,CD95,CD271,CD270,CD90,CD272,CD16,CD14,CD15,CD13,CD267,CD26-2,CD200,CD26-1,CD18,CD19,CD194,TCR-1,TCR-V-7.2,CD70,CD71,CD72,CD73,TCR-V-2,CD177,CD301,CD140a,CD140b,CD305,CD304,CD2,CD309,CD85g,CD110,CD8,CD9,HLA-DR,CD137,CD134,CD135,CD61,CD192,CD268,CD269,CD81,CD80,CD83,CD193,TSLPR,CD86,CCR10,Notch-1,Notch-2,CD337,CD79b,CD275-2,CD275-1,CD79a,CD49b,CD64,CD63,CD45-2,CD45-1,CD133-2,CD133-1,CD69,CD68,CD314,CD11a/CD18,CD186,CD185,CD184,CD103,CD102,CD106,CD105,CD66b,CD252,CD109,TCR-V-24-J-18,Integrin-7,CD158b,CD158f,CD8a,CD203c,CD52,CD195,CD196,CD57,CD54,CD55,CD99,CD59,CD199,CD93,CD244,CD158,CD235ab,CD273,CD243,CD325,CD324,CD307e,CD172a,CD307d,CD42b,CD115,CD117,XCR1,CD112,MERTK,B7-H4,CD21,CD307c/FcRL3,CLEC2,CD48,CD47,CD46,CD41,CD40,CD43,CD338,CD235a,CD335,CD3-2,CD119,CD169,CD28,CD161,CD163,CD138-1,CD164,CD138-2,CD144,CD202b,CD11c""".split(',')

def main():
    mtx=os.path.join(D,'GSM5008738_ADT_3P-matrix.mtx.gz')
    with gzip.open(mtx,'rt') as fh:
        for line in fh:
            if line.startswith('%'): continue
            nrow,ncol,nnz=(int(x) for x in line.split()); break
    print('matrix: %d proteins x %d cells, %d nonzeros' % (nrow,ncol,nnz))
    assert nrow==len(NAMES), (nrow,len(NAMES))
    df=pd.read_csv(mtx, sep=' ', skiprows=3, header=None, names=['r','c','v'],
                   dtype={'r':np.int32,'c':np.int32,'v':np.float32}, engine='c')
    print('parsed %d entries' % len(df))
    M=np.zeros((ncol,nrow), dtype=np.float32)
    M[df['c'].values-1, df['r'].values-1] = df['v'].values
    del df
    meta=pd.read_csv(os.path.join(D,'GSE164378_sc.meta.data_3P.csv.gz'))
    meta.rename(columns={meta.columns[0]:'cell'}, inplace=True)
    print('meta: %d rows; donors %s' % (len(meta), sorted(meta['donor'].unique())))
    print('celltype.l1 counts:'); print(meta['celltype.l1'].value_counts().to_string())
    assert len(meta)==ncol, (len(meta), ncol)
    np.savez_compressed(os.path.join(D,'healthy_adt.npz'), adt=M, names=np.array(NAMES),
                        celltype_l1=meta['celltype.l1'].astype(str).values,
                        celltype_l2=meta['celltype.l2'].astype(str).values,
                        donor=meta['donor'].astype(str).values)
    print('saved healthy_adt.npz')

if __name__=='__main__': main()
