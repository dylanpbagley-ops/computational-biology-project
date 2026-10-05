import os, sys
sys.path.insert(0, os.path.expanduser('~/mnt/outputs/ev_capture_model/src'))
import aim12_single_cell_cytof_reordering as a12
D=os.path.expanduser('~/mnt/outputs/ev_capture_data')
def pan(p):
    raw=open(p,'rb').read(400000); h=a12.parse_fcs3_header(raw); tk,st=a12.parse_fcs3_text_segment(raw,h)
    return tk, {i:st[i] for i in sorted(st)}
refs = {
 'DLB_CytB2': f'{D}/raw_cytof/CytB2_DLB_P5904.fcs',
 'DLB_CytB4': f'{D}/raw_cytof/CytB4_DLB_O5296.fcs',
 'DLB_CytB6': f'{D}/raw_cytof/CytB6_DLB_R5642.fcs',
 'RLN_CytB2': f'{D}/raw_cytof_rln/CytB2_RLN_O8067.fcs',
 'RLN_CytB4': f'{D}/raw_cytof_rln/CytB4_RLN_K5851.fcs',
 'RLN_CytB6': f'{D}/raw_cytof_rln/CytB6_RLN_G6768.fcs',
}
pans={}
for k,p in refs.items():
    tk,st=pan(p); pans[k]=st
    print('%-11s PAR=%s  n_stains=%d' % (k, tk['PAR'], len(st)))
a=set(pans['DLB_CytB2'].values()); 
for k in pans:
    s=set(pans[k].values())
    print('\n%-11s missing vs DLB_CytB2: %s' % (k, sorted(a-s)))
    print('%-11s extra   vs DLB_CytB2: %s' % ('', sorted(s-a)))
# check index alignment for the antigens we use
ANT=['CD19','CD20','CD22','CD79B','CD38','CD3','CD45','193Ir-DNA','195Pt-Viability','IgKap','IgLam','HLA-DR','CD44','CD27','CD48','CD72','CD40','CD21','CXCR5']
print('\nindex alignment (find_parameter_index_by_stain_substring):')
print('%-16s' % 'antigen', ''.join('%-12s'%k for k in refs))
for m in ANT:
    row=[]
    for k in refs:
        try: row.append(str(a12.find_parameter_index_by_stain_substring(pans[k], m)))
        except Exception: row.append('ABSENT')
    mark='' if len(set(row))==1 else '   <== MISALIGNED'
    print('%-16s' % m, ''.join('%-12s'%x for x in row), mark)
