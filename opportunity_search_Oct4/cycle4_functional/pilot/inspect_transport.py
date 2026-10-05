from __future__ import annotations
import csv, gzip, hashlib, json, math, re, sys, collections
from pathlib import Path
import numpy as np
import openpyxl
from transport_helpers import normalized_log_frequency,summarize_envelope,guide_filter_reasons

STAGE=Path(__file__).resolve().parents[1]
PILOT=Path(__file__).resolve().parent
ROOT=STAGE.parents[1]
DATA=STAGE/'data'

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):
            h.update(chunk)
    return h.hexdigest()
def fail(message):
    raise SystemExit('STOP: '+message)
def read_json(path):
    def reject(v): raise ValueError('non-finite JSON '+v)
    return json.loads(path.read_text(encoding='utf-8'),parse_constant=reject)
def write_tsv(path,header,rows):
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n',quoting=csv.QUOTE_MINIMAL)
        w.writerow(header); w.writerows(rows)
def text_value(v):
    if v is None: return ''
    if isinstance(v,float):
        if not math.isfinite(v): return ''
        return format(v,'.17g')
    return str(v)
def finite_num(v):
    try:
        x=float(v)
        return x if math.isfinite(x) else None
    except (TypeError,ValueError):
        return None

def verify_protocol_and_inputs():
    protocol=read_json(PILOT/'protocol.json')
    check=dict(protocol); digest=check.pop('protocol_sha256'); check.pop('protocol_hash_definition',None)
    canonical=json.dumps(check,sort_keys=True,separators=(',',':'),ensure_ascii=False,allow_nan=False).encode('utf-8')
    if hashlib.sha256(canonical).hexdigest()!=digest: fail('protocol SHA mismatch')
    pmeta=read_json(PILOT/'preservation_before.json')
    if pmeta.get('protocol',{}).get('sha256')!=digest: fail('protocol hash differs from pilot preservation record')
    if pmeta.get('reviewed_plan',{}).get('sha256')!=protocol['reviewed_plan']['sha256']: fail('reviewed plan hash differs from pilot preservation record')
    for item in protocol['inputs']:
        path=STAGE/item['path']
        if not path.is_file() or path.stat().st_size!=item['bytes'] or sha(path)!=item['sha256']:
            fail('protected data input changed: '+item['path'])
    for item in protocol['frozen_prior_art']:
        path=STAGE/item['path']
        if not path.is_file() or path.stat().st_size!=item['bytes'] or sha(path)!=item['sha256']:
            fail('frozen prior-art input changed: '+item['path'])
    plan=ROOT/protocol['reviewed_plan']['path']
    if not plan.is_file() or sha(plan)!=protocol['reviewed_plan']['sha256']:
        fail('reviewed plan changed')
    return protocol,digest

def read_metadata(protocol,count_names):
    path=DATA/'sample_mapping.tsv'
    with path.open('r',encoding='utf-8',newline='') as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    count_col_to_index={name:i for i,name in enumerate(count_names)}
    context_rows=[]
    for spec in protocol['sample_mapping']['source_contexts']:
        exact=[r for r in rows if r['source']==spec['source']]
        if not exact: fail('source context absent: '+spec['source'])
        t0=[r for r in exact if r['treatment']=='none' and r['time_point']==spec['t0_timepoint']]
        if len(t0)!=1: fail('expected one exact model T0 for '+spec['model'])
        def group(treatment,tp):
            selected=[r for r in exact if r['treatment']==treatment and r['time_point']==tp]
            selected=sorted(selected,key=lambda r:r['replicate_label'])
            if [r['replicate_label'] for r in selected]!=['A','B','C']:
                fail(f'{spec["model"]} {treatment} {tp}: expected unique A/B/C cultures; got {[r["replicate_label"] for r in selected]}')
            idx=[]
            for r in selected:
                c=r['count_column']
                if c not in count_col_to_index: fail('sample map column absent from count header: '+c)
                idx.append(count_col_to_index[c])
            return idx
        one={'spec':spec,'t0_index':count_col_to_index[t0[0]['count_column']], 't0_column':t0[0]['count_column'],'early_treated':group(spec['treated_treatment'],spec['early']['treated_timepoint']),'early_control':group(spec['control_treatment'].replace(' (WT-CD8 per source paper)',''),spec['early']['control_timepoint']),'late_treated':group(spec['treated_treatment'],spec['late']['treated_timepoint']),'late_control':group(spec['control_treatment'].replace(' (WT-CD8 per source paper)',''),spec['late']['control_timepoint'])}
        # Replicate culture labels and source context are recorded before using the count values.
        context_rows.append(one)
    return rows,context_rows,count_col_to_index

def main():
    protocol,protocol_sha=verify_protocol_and_inputs()
    # Open only the count header and exact metadata mapping first; protocol already frozen.
    count_path=DATA/'GSE149933_mTKO_readCounts_4May20.txt.gz'
    with gzip.open(count_path,'rt',encoding='utf-8',newline='') as f:
        reader=csv.reader(f,delimiter='\t'); header=next(reader)
        required=['Guide.ID','Gene.symbol']
        if any(k not in header for k in required): fail('required named count field missing')
        guide_i=header.index('Guide.ID'); gene_i=header.index('Gene.symbol')
        count_cols=[i for i,name in enumerate(header) if name.startswith('MUS_')]
        count_names=[header[i] for i in count_cols]
        if len(count_cols)!=115 or len(set(count_names))!=115: fail('count matrix must have 115 unique MUS columns')
        metadata_rows,contexts,count_col_index=read_metadata(protocol,count_names)
        expected_guides=94528
        counts=np.empty((expected_guides,len(count_cols)),dtype=np.int64)
        guide_ids=[]; gene_symbols=[]; total_rows=0
        for row in reader:
            if len(row)!=len(header): fail(f'malformed count row {total_rows+2}: width={len(row)}')
            if total_rows>=expected_guides: fail('unexpected extra guide row')
            guide_ids.append(row[guide_i]); gene_symbols.append(row[gene_i])
            try: counts[total_rows,:]=[int(row[i]) for i in count_cols]
            except ValueError: fail(f'noninteger read count at source row {total_rows+2}')
            total_rows+=1
        if total_rows!=expected_guides: fail(f'expected {expected_guides} guide rows, observed {total_rows}')
    if (counts<0).any(): fail('negative raw read count')
    count_totals=counts.sum(axis=0,dtype=np.int64)
    logfreq=normalized_log_frequency(counts)
    guide_id_counts=collections.Counter(guide_ids)
    symbol_to_indices=collections.defaultdict(list)
    for i,g in enumerate(gene_symbols):
        if g: symbol_to_indices[g].append(i)

    lib_path=DATA/'GSE149933_mTKO_library.txt.gz'
    with gzip.open(lib_path,'rt',encoding='utf-8',newline='') as f:
        reader=csv.DictReader(f,delimiter='\t')
        lib_header=reader.fieldnames or []
        if lib_header!=['SEQUENCE','GUIDE_ID','GENE','EXON','SeqScore']: fail('library annotation header changed')
        library_rows=list(reader)
    if len(library_rows)!=expected_guides: fail('guide annotation row count changed')
    libid_counts=collections.Counter(r['GUIDE_ID'] for r in library_rows)
    lib_by_id={}
    for r in library_rows:
        if libid_counts[r['GUIDE_ID']]==1: lib_by_id[r['GUIDE_ID']]=r['GENE']

    workbook=openpyxl.load_workbook(DATA/'Ramos_SupplementaryData1.xlsx',read_only=True,data_only=True)
    ws=workbook.active
    xrows=ws.iter_rows(values_only=True)
    xheader=list(next(xrows))
    columns=protocol['ramos_columns']['primary_in_vitro']
    needed=['Gene']+columns
    if any(c not in xheader for c in needed): fail('exact Ramos workbook header changed')
    pos={name:xheader.index(name) for name in needed}
    panel=[]
    for i,row in enumerate(xrows,2):
        if len(row)!=len(xheader): fail('Ramos workbook row width changed')
        panel.append({'panel_row':i,'gene':row[pos['Gene']],'raw_1to2_lfc':row[pos[columns[0]]],'raw_1to2_fdr':row[pos[columns[1]]],'raw_1to10_lfc':row[pos[columns[2]]],'raw_1to10_fdr':row[pos[columns[3]]]})
    workbook.close()
    if len(panel)!=1021: fail(f'expected 1021 selected validation rows, saw {len(panel)}')
    panel_counts=collections.Counter(p['gene'] for p in panel if p['gene'] not in (None,''))
    panel_symbols=sorted(panel_counts)
    panel_genes=set(panel_symbols)

    # Validate the exact six contexts and source-defined early/latest pairs before calculations.
    for one in contexts:
        spec=one['spec']
        for key in ['early_treated','early_control','late_treated','late_control']:
            if len(one[key])!=3: fail(f'{spec["model"]}: missing matched culture replicate group {key}')
    actual_renca=[r for r in metadata_rows if r['source']=='Renca cell line' and r['treatment'] in ('IFN-gamma','TNF-alpha')]
    # This exact source is deliberately never joined into the Renca-HA CTL context.

    guide_eligibility=[]
    effects=[]
    results_by_gene=collections.defaultdict(dict)
    eligibility_totals=collections.Counter()
    guide_i_by_symbol={g:symbol_to_indices.get(g,[]) for g in panel_genes}
    # Every exact panel gene is evaluated in each model; missing genes retain explicit status rows.
    for one in contexts:
        spec=one['spec']; model=spec['model']
        t0_idx=one['t0_index']; late_ctrl_idx=one['late_control']
        retained_by_gene={}; filter_by_guide={}
        for gene in panel_symbols:
            keep=[]; indices=guide_i_by_symbol[gene]
            for idx in indices:
                gid=guide_ids[idx]
                library_gene=lib_by_id.get(gid)
                reasons=guide_filter_reasons(gid,gene,library_gene,int(counts[idx,t0_idx]),[int(counts[idx,j]) for j in late_ctrl_idx],duplicate_id=(guide_id_counts[gid]!=1),missing_t0=False,missing_control=False)
                if not reasons: keep.append(idx)
                filter_by_guide[idx]=(reasons,library_gene)
                control_counts=[int(counts[idx,j]) for j in late_ctrl_idx]
                guide_eligibility.append([gene,model,spec['source'],gid,gene_symbols[idx],library_gene or '',int(counts[idx,t0_idx]),*control_counts,';'.join(reasons) if reasons else 'eligible',not bool(reasons)])
                eligibility_totals.update(reasons if reasons else ['eligible'])
            retained_by_gene[gene]=keep
        for p in panel:
            gene=p['gene']
            if gene in (None,''): continue
            duplicate_panel=(panel_counts[gene]!=1)
            indices=guide_i_by_symbol.get(gene,[])
            keep=retained_by_gene.get(gene,[])
            if duplicate_panel:
                gene_status='nonunique_panel_symbol'
            elif not indices:
                gene_status='no_exact_count_gene_symbol_match'
            elif len(keep)<3:
                gene_status='fewer_than_three_retained_unique_guides'
            else:
                gene_status='eligible'
            for time_key,role in [('early','early'),('late','primary_latest')]:
                t_idx=one[time_key+'_treated']; c_idx=one[time_key+'_control']
                pair_values=[]; med=lo=hi=None; cls='not_estimable'
                if gene_status=='eligible':
                    for ti in t_idx:
                        for ci in c_idx:
                            d=logfreq[keep,ti]-logfreq[keep,ci]
                            pair_values.append(float(np.median(d)))
                    med,lo,hi,cls=summarize_envelope(pair_values)
                baseline=None
                if gene_status=='eligible':
                    baseline_differences=logfreq[np.ix_(keep,late_ctrl_idx)]-logfreq[np.ix_(keep,[t0_idx])]
                    baseline=float(np.median(baseline_differences))
                model_result={'panel_row':p['panel_row'],'model':model,'source':spec['source'],'gene_status':gene_status,'n_guides_in_count_matrix':len(indices),'n_retained_unique_guide_ids':len(keep),'median':med,'min':lo,'max':hi,'class':cls,'pair_values':pair_values,'baseline_control_vs_T0_median':baseline,'treated_timepoint':spec[time_key]['treated_timepoint'],'control_timepoint':spec[time_key]['control_timepoint']}
                results_by_gene[gene][(model,role)]=model_result
                effects.append([p['panel_row'],gene,model,spec['source'],role,spec[time_key]['treated_timepoint'],spec[time_key]['control_timepoint'],gene_status,len(indices),len(keep),len(set(guide_ids[i] for i in keep)),*(pair_values if pair_values else ['']*9),med if med is not None else '',lo if lo is not None else '',hi if hi is not None else '',cls,baseline if baseline is not None else '','single T0; baseline diagnostic only; not subtracted from paired CTL effect'])

    effect_header=['panel_row','gene','model','source_context','timepoint_role','treated_timepoint','control_timepoint','gene_status','n_count_rows_for_symbol','n_retained_unique_guides','n_unique_guide_ids','pair_1_treatedA_controlA','pair_2_treatedA_controlB','pair_3_treatedA_controlC','pair_4_treatedB_controlA','pair_5_treatedB_controlB','pair_6_treatedB_controlC','pair_7_treatedC_controlA','pair_8_treatedC_controlB','pair_9_treatedC_controlC','median_of_9','minimum_of_9','maximum_of_9','direction_class','primary_control_vs_T0_baseline_median','baseline_note']
    write_tsv(PILOT/'gene_effects.tsv',effect_header,effects)
    guide_header=['gene','model','source_context','guide_id','count_matrix_gene_symbol','library_gene_symbol','T0_raw_count','primary_control_A_raw_count','primary_control_B_raw_count','primary_control_C_raw_count','exclusion_reasons','eligible']
    write_tsv(PILOT/'guide_eligibility.tsv',guide_header,guide_eligibility)

    flag_rows=[]; primary_hits=[]; symmetric_hits=[]; anchor_results=[]
    anchor_set=set(protocol['anchors'])
    for p in panel:
        gene=p['gene']; gene_ok=gene not in (None,'') and panel_counts[gene]==1
        late=[]; early=[]
        for one in contexts:
            model=one['spec']['model']
            late.append(results_by_gene.get(gene,{}).get((model,'primary_latest')))
            early.append(results_by_gene.get(gene,{}).get((model,'early')))
        late_sens=[r['model'] for r in late if r and r['gene_status']=='eligible' and r['class']=='sign_stable_sensitizing']
        late_res=[r['model'] for r in late if r and r['gene_status']=='eligible' and r['class']=='sign_stable_resistant']
        late_mixed=[r['model'] for r in late if r and r['gene_status']=='eligible' and r['class']=='mixed']
        early_sens=[r['model'] for r in early if r and r['gene_status']=='eligible' and r['class']=='sign_stable_sensitizing']
        consistent=[]; flips=[]
        for l,e in zip(late,early):
            if not l or not e or l['gene_status']!='eligible' or e['gene_status']!='eligible': continue
            stable={'sign_stable_sensitizing','sign_stable_resistant'}
            if l['class'] in stable and e['class']==l['class']: consistent.append(l['model'])
            if (l['class'],e['class']) in [('sign_stable_sensitizing','sign_stable_resistant'),('sign_stable_resistant','sign_stable_sensitizing')]: flips.append(l['model'])
        lfc2=finite_num(p['raw_1to2_lfc']); fdr2=finite_num(p['raw_1to2_fdr']); lfc10=finite_num(p['raw_1to10_lfc']); fdr10=finite_num(p['raw_1to10_fdr'])
        res2=(lfc2 is not None and fdr2 is not None and lfc2>0 and fdr2<=0.05)
        res10=(lfc10 is not None and fdr10 is not None and lfc10>0 and fdr10<=0.05)
        sens2=(lfc2 is not None and fdr2 is not None and lfc2<0 and fdr2<=0.05)
        sens10=(lfc10 is not None and fdr10 is not None and lfc10<0 and fdr10<=0.05)
        primary=bool(gene_ok and len(late_sens)>=3 and res2 and res10)
        symmetric=bool(gene_ok and len(late_res)>=3 and sens2 and sens10)
        if primary: primary_hits.append(gene)
        if symmetric: symmetric_hits.append(gene)
        anchor=(gene in anchor_set)
        if anchor:
            anchor_results.append({'gene':gene,'late_sensitizing_models':late_sens,'late_resistant_models':late_res,'late_mixed_models':late_mixed,'early_late_same_stable_models':consistent,'early_late_stable_sign_flip_models':flips,'ramos_1to2_lfc':lfc2,'ramos_1to2_fdr':fdr2,'ramos_1to10_lfc':lfc10,'ramos_1to10_fdr':fdr10,'ramos_supported_resistance_both':bool(res2 and res10),'ramos_supported_sensitization_both':bool(sens2 and sens10)})
        notes=[]
        if not gene_ok: notes.append('nonunique_or_missing_panel_symbol')
        if not symbol_to_indices.get(gene): notes.append('no_exact_count_gene_symbol_match')
        if len(late_sens)<3: notes.append('fewer_than_three_late_sign_stable_sensitizing_models')
        if not (res2 and res10): notes.append('Ramos_resistance_rule_not_met_in_both_in_vitro_conditions')
        flag_rows.append([p['panel_row'],text_value(gene),gene_ok,len(symbol_to_indices.get(gene,[])),len(late_sens),';'.join(late_sens),len(late_res),';'.join(late_res),len(late_mixed),';'.join(late_mixed),len(consistent),';'.join(consistent),';'.join(flips),text_value(p['raw_1to2_lfc']),text_value(p['raw_1to2_fdr']),text_value(p['raw_1to10_lfc']),text_value(p['raw_1to10_fdr']),res2,res10,sens2,sens10,primary,symmetric,anchor,';'.join(notes)])
    flag_header=['panel_row','gene','unique_panel_symbol','count_matrix_rows_for_exact_symbol','late_sign_stable_sensitizing_model_count','late_sign_stable_sensitizing_models','late_sign_stable_resistant_model_count','late_sign_stable_resistant_models','late_mixed_model_count','late_mixed_models','early_late_same_stable_sign_model_count','early_late_same_stable_sign_models','early_late_stable_sign_flip_models','Ramos_1to2_LFC_raw','Ramos_1to2_FDR_raw','Ramos_1to10_LFC_raw','Ramos_1to10_FDR_raw','Ramos_supported_resistance_1to2','Ramos_supported_resistance_1to10','Ramos_supported_sensitization_1to2','Ramos_supported_sensitization_1to10','primary_transport_flag','symmetric_descriptive_flag','predefined_anchor','notes']
    write_tsv(PILOT/'transport_flags.tsv',flag_header,flag_rows)

    # Verify expected source/model and control mapping; report raw-data quality and exact exclusions.
    missing_symbol=sum(1 for g in gene_symbols if not g)
    duplicate_guide_rows=sum(n for n in guide_id_counts.values() if n>1)
    library_mismatch_rows=sum(1 for gid,g in zip(guide_ids,gene_symbols) if gid in lib_by_id and g and lib_by_id[gid]!=g)
    panel_match_count=sum(1 for g,n in panel_counts.items() if n==1 and symbol_to_indices.get(g))
    eligible_gene_model_pairs=sum(1 for x in results_by_gene.values() for r in x.values() if r['gene_status']=='eligible')
    early_unstable=[]
    for gene,mp in results_by_gene.items():
        for one in contexts:
            m=one['spec']['model']; e=mp.get((m,'early')); l=mp.get((m,'primary_latest'))
            if e and l and e['gene_status']=='eligible' and l['gene_status']=='eligible' and (e['class'],l['class']) in [('sign_stable_sensitizing','sign_stable_resistant'),('sign_stable_resistant','sign_stable_sensitizing')]:
                early_unstable.append({'gene':gene,'model':m,'early':e['class'],'latest':l['class']})
    summary={
      'schema_version':1,'pilot_start_utc':protocol['actual_pilot_start_utc'],'computed_at_utc':__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(timespec='seconds').replace('+00:00','Z'),
      'protocol_sha256':protocol_sha,'reviewed_plan_sha256':protocol['reviewed_plan']['sha256'],
      'scientific_reproduction_checks':{'protocol_hash_verified':True,'reviewed_plan_hash_verified':True,'protected_raw_input_hashes_verified':True,'protected_prior_art_hashes_verified':True,'raw_input_hashes':{x['path']:x['sha256'] for x in protocol['inputs']},'prior_art_hashes':{x['path']:x['sha256'] for x in protocol['frozen_prior_art']},'count_matrix_sha256':sha(count_path),'library_sha256':sha(DATA/'GSE149933_mTKO_library.txt.gz'),'ramos_workbook_sha256':sha(DATA/'Ramos_SupplementaryData1.xlsx')},
      'input_counts':{'count_matrix_guide_rows':int(total_rows),'count_matrix_sample_columns':len(count_cols),'library_annotation_rows':len(library_rows),'sample_mapping_rows':len(metadata_rows),'mapped_count_columns_used':len(set(sum([one['early_treated']+one['early_control']+one['late_treated']+one['late_control']+[one['t0_index']] for one in contexts],[]))),'ramos_selected_panel_rows':len(panel),'unique_nonempty_panel_symbols':len(panel_symbols),'unique_panel_symbols_exactly_matched_to_count_matrix':int(panel_match_count)},
      'mapping':{'contexts':[{**one['spec'],'t0_column':one['t0_column'],'early_treated_columns':[count_names[i] for i in one['early_treated']],'early_control_columns':[count_names[i] for i in one['early_control']],'late_treated_columns':[count_names[i] for i in one['late_treated']],'late_control_columns':[count_names[i] for i in one['late_control']]} for one in contexts],'antigen_free_renca_cytokine_sample_rows_excluded':len(actual_renca)},
      'quality':{'missing_count_matrix_gene_symbols':int(missing_symbol),'duplicate_guide_id_rows':int(duplicate_guide_rows),'library_gene_symbol_mismatch_rows':int(library_mismatch_rows),'guide_filter_reason_counts':dict(eligibility_totals),'eligible_gene_model_timepoint_rows':int(eligible_gene_model_pairs),'early_latest_stable_direction_flips':early_unstable},
      'primary_test':{'primary_flag_count':len(primary_hits),'primary_flag_genes':primary_hits,'symmetric_descriptive_flag_count':len(symmetric_hits),'symmetric_descriptive_flag_genes':symmetric_hits,'rule':'>=3/6 latest CTL contexts sign-stable sensitizing AND positive Ramos LFC/FDR<=0.05 in both in-vitro conditions'},
      'predefined_anchor_results':anchor_results,
      'interpretive_limits':['culture-pair envelopes are observed descriptive minima/maxima, not confidence intervals','six cell lines and one CAR-T leukemia model are the cross-study units; guides/genes/cultures are not independent tumors','source-defined control harvest days sometimes differ from CTL days (for example EMT6 T15/T19 and Renca-HA T22/T24), so growth/timing confounding remains','Renca-HA CTL rows never pool with antigen-free Renca cytokine rows','selection frequency is not direct percent killing and normalization differs from published drugZ/MAGeCK','selected panel, single B-ALL model, one study contrast and mouse-only data preclude clinical inference','cross-study differences cannot identify whether modality, tumor background, antigen, timing, or another study factor caused a pattern'],
      'decision':'Pilot only; parent must determine whether any exact pattern is distinct and biologically substantial. No automatic project selection.'
    }
    # JSON must be finite; values were normalized to Python builtins above.
    (PILOT/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({'protocol_sha256':protocol_sha,'primary_flags':len(primary_hits),'symmetric_flags':len(symmetric_hits),'panel_rows':len(panel),'exact_symbol_matches':panel_match_count,'eligible_gene_model_timepoint_rows':eligible_gene_model_pairs,'missing_gene_symbols':missing_symbol,'duplicate_guide_rows':duplicate_guide_rows,'library_mismatch_rows':library_mismatch_rows,'antigen_free_renca_excluded_rows':len(actual_renca)},ensure_ascii=True))

if __name__=='__main__': main()
