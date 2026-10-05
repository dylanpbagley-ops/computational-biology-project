"""
Cross-check: does every headline number quoted in the results document actually appear in the
result JSON the analysis wrote?

This exists because the project has already been bitten once by prose drifting away from code
(Aim 12's counting rule, recorded in §20). Rather than trusting that the write-up was updated,
this script re-reads the JSON and re-asserts every number sections 23-25 quote. It fails loudly
on any mismatch.
"""
import json, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DOCS = os.path.abspath(os.path.join(HERE, '..', '..', 'docs'))
DOC = os.path.join(DOCS, 'Aim13_Results_Two_Comparator_Ceiling.md')

W = json.load(open(os.path.join(HERE, 'wholebody_v2.json')))
P = json.load(open(os.path.join(HERE, 'pancreatic_thread.json')))
N = json.load(open(os.path.join(HERE, 'two_antigen_null_control.json')))
text = open(DOC, encoding='utf-8').read()

fails, checks = [], 0


def check(label, value, fmt='%.3f', in_doc=True):
    """Assert the value, formatted as the document formats it, appears in the document."""
    global checks
    checks += 1
    s = fmt % value
    ok = (s in text) if in_doc else True
    if not ok:
        fails.append('%-52s expected to find "%s" in the results document' % (label, s))
    print('  %-52s %-12s %s' % (label, s, 'ok' if ok else 'MISSING FROM DOC'))


print('=' * 78)
print('Section 23 — whole-body restriction')
print('=' * 78)
t = W['tests']
check('one-factor separation (top_share)', t['A_one_factor_top_share'])
check('two-factor separation', t['A_two_factor'])
check('two-factor p', t['A_two_factor_p'], '%.5f')
check('leave-one-out minimum separation', t['loo_sep_min'])
check('leave-one-out worst p', t['loo_p_max'], '%.5f')
check('solid-tumour-only separation', W['solid_only']['sep'])
check('solid-tumour-only exact p', W['solid_only']['p'], '%.4f')
check('germ-removed two-factor', W['germ_removed']['two_factor'])
check('reference B two-factor separation', t['B_two_factor'])
check('Spearman between references', t['spearman_A_vs_B'])

lab = W['labels']
n_safe = sum(1 for g in lab if lab[g] == 'safe')
n_tox = sum(1 for g in lab if lab[g] == 'toxic')
n_eff = sum(1 for g in lab if lab[g] == 'efficacy_fail')
print('\n  panel composition: %d safe, %d toxic, %d efficacy-fail' % (n_safe, n_tox, n_eff))
for n, want in ((n_safe, 19), (n_tox, 11), (n_eff, 12)):  # updated 2026-09-10: ERBB3, AXL added
    checks += 1
    if n != want:
        fails.append('panel composition changed: got %d, document says %d' % (n, want))

# section 23.4 -- the specificity control (safe vs efficacy-fail). Not expected to be near
# chance any more (see the results document's honest discussion of why); the checker's job is
# only to catch it drifting from whatever the document currently states, not to enforce a
# particular value.
cA = W['reference_A_cell_type']
safe_ts = [cA[g]['top_share'] for g in cA if lab[g] == 'safe']
eff_ts = [cA[g]['top_share'] for g in cA if lab[g] == 'efficacy_fail']
def _sep(a, b):
    gt = sum(1 for x in a for y in b if x > y)
    eq = sum(1 for x in a for y in b if x == y)
    return (gt + 0.5 * eq) / (len(a) * len(b))
sep_ef = _sep(safe_ts, eff_ts)
check('section 23.4 specificity: safe vs efficacy-fail separation', sep_ef)

# the theorem-adjacent invariant: the two-factor rule must never rank a solid-tissue target
# above a blood-lineage one, whatever its confinement.
cA = W['reference_A_cell_type']; haem = W['haem_home_A']
sc = {g: cA[g]['top_share'] + (1.0 if haem[g] else 0.0) for g in cA}
blood_min = min(sc[g] for g in sc if haem[g])
solid_max = max(sc[g] for g in sc if not haem[g])
checks += 1
print('\n  lexicographic invariant: worst blood %.3f >= best solid %.3f' % (blood_min, solid_max))
if blood_min < solid_max:
    fails.append('two-factor rule is not lexicographic: %.3f < %.3f' % (blood_min, solid_max))

print('\n' + '=' * 78)
print('Section 24 — pancreatic thread')
print('=' * 78)
check('CLDN18.2 calibration Spearman r', P['cldn182']['spearman_r'])
check('CLDN18.2 RMSE', P['cldn182']['rmse'])
check('CPTAC log-rank p', P['cptac']['logrank_p'], '%.3f')
checks += 1
if P['cptac']['n'] != 135:
    fails.append('CPTAC n changed from 135 to %d' % P['cptac']['n'])
rank = P['ranking']
checks += 1
if rank[0][1] != 'claudin-18.2':
    fails.append('CLDN18.2 no longer ranks first in the PDAC panel')
msln_rank = [i for i, r in enumerate(rank, 1) if r[0] == 'MSLN'][0]
print('  mesothelin PDAC rank: %d of %d  (document says 5th — the corrected claim)'
      % (msln_rank, len(rank)))
checks += 1
if msln_rank != 5:
    fails.append('mesothelin PDAC rank changed from 5 to %d; §24.3 must be updated' % msln_rank)

print('\n' + '=' * 78)
print('Section 25 — shuffled-label null control')
print('=' * 78)
for k, lbl in (('single', 'best single-antigen J'),
               ('stat_gain', 'statistical gain'),
               ('eng_gain', 'engineering gain')):
    r = float(np.median(N['real'][k])); n_ = float(np.median(N['null'][k]))
    print('  %-30s real %+.4f   shuffled %+.4f' % (lbl, r, n_))
checks += 1
if float(np.median(N['null']['single'])) > 0.05:
    fails.append('the null control is broken: shuffled single-antigen J is not near zero')
checks += 1
if float(np.median(N['real']['eng_gain'])) <= 10 * float(np.median(N['null']['eng_gain']) + 1e-9):
    fails.append('engineering gain is no longer clearly above the shuffled floor')

print('\n' + '=' * 78)
print('September 2026 expansion angles — new headline numbers')
print('=' * 78)

_here = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    with open(os.path.join(_here, name)) as fh:
        return json.load(fh)


# --- Angle 3: bootstrap confidence intervals ---------------------------------------
try:
    CI = _load('auc_confidence_intervals.json')
    for _k, _pt in (('A_one_factor_top_share', 0.8804), ('A_two_factor', 0.9139),
                    ('B_two_factor', 0.7512), ('solid_only_one_factor', 0.8125)):
        checks += 1
        _v = CI[_k]
        print('  CI %-26s %.3f  [%.3f, %.3f]' % (_k, _v['point'], _v['ci_lo'], _v['ci_hi']))
        if abs(_v['point'] - _pt) > 2e-3:
            fails.append('CI point estimate for %s drifted: %.4f vs %.4f'
                         % (_k, _v['point'], _pt))
        if not (_v['ci_lo'] <= _v['point'] <= _v['ci_hi']):
            fails.append('CI for %s does not bracket its point estimate' % _k)
    checks += 1
    if CI['B_two_factor']['ci_lo'] >= 0.60:
        fails.append('reference B lower bound no longer reaches toward chance; §23.7 '
                     'must stop being described as weakly constrained')
except FileNotFoundError:
    fails.append('auc_confidence_intervals.json missing — run auc_confidence_intervals.py')

# --- Angle 3b: label sensitivity ----------------------------------------------------
try:
    LS = _load('label_sensitivity.json')
    _s = LS['_summary']
    checks += 1
    print('  label sensitivity range %.3f-%.3f, worst p %.5f'
          % (_s['sep_min'], _s['sep_max'], _s['p_max']))
    if abs(_s['sep_min'] - 0.835) > 0.005 or abs(_s['sep_max'] - 0.926) > 0.005:
        fails.append('label-sensitivity range moved from 0.835-0.926')
    checks += 1
    if _s['p_max'] >= 0.01:
        fails.append('a label relabeling now breaks significance')
except FileNotFoundError:
    fails.append('label_sensitivity.json missing — run label_sensitivity.py')

# --- Angle 1: continuous weights (a PINNED NEGATIVE result) -------------------------
try:
    CW = _load('continuous_weight_sensitivity.json')
    checks += 1
    print('  binary %.3f vs continuous %.3f (continuous is expected to be WORSE)'
          % (CW['binary']['sep'], CW['continuous']['sep']))
    if CW['continuous']['sep'] >= CW['binary']['sep']:
        fails.append('continuous weighting now BEATS binary — this reverses a documented '
                     'negative result and must be re-examined, not silently adopted')
    checks += 1
    if abs(CW['binary']['sep'] - 0.9139) > 2e-3:
        fails.append('binary baseline drifted from the section 23 value')
except FileNotFoundError:
    fails.append('continuous_weight_sensitivity.json missing')

# --- Angle 1b: re-ranked screen + widened exclusion set ------------------------------
try:
    RR = _load('screen_rerank_weighted.json')
    checks += 1
    print('  screen re-rank: %d candidates, %d flagged by widened exclusion, %d new'
          % (RR['n_candidates'], RR['n_flagged_by_widened_exclusion'], RR['n_genuinely_new']))
    for _g, _lim in (('DCT', 10), ('ERVFRD-1', 10), ('PSG1', 10)):
        if RR['by_gene'][_g]['rank_weighted'] > _lim:
            fails.append('%s fell out of the re-ranked top %d' % (_g, _lim))
    checks += 1
    if RR['by_gene']['CD3D']['rank_weighted'] <= 10:
        fails.append('CD3D (a T-cell gene) is back in the re-ranked top 10')
    checks += 1
    if RR['n_flagged_by_widened_exclusion'] < 6:
        fails.append('widened exclusion set no longer catches the 6 known-target '
                     'candidates; novelty count would be overstated again')
except FileNotFoundError:
    fails.append('screen_rerank_weighted.json missing')

# --- September 28 strengthening pass ------------------------------------------------
print('\n' + '=' * 78)
print('September 28 strengthening pass')
print('=' * 78)
try:
    OR_ = _load('outcome_test_robustness.json')
    # section 27 quotes these; the document was written on an older panel and is stale
    check('section 27 replication: confinement separation (live panel)',
          OR_['autopsy_replication']['point'])
    check('section 27.3 registry control separation (live registry)',
          OR_['registry_control_phase3']['point'])
    checks += 1
    _e = OR_['expendability_direct_tests']
    print('  expendability Fisher p %.4f, perm p %.4f, paired AUC gain %+.3f [%.3f, %.3f]'
          % (_e['fisher_exact_one_sided_p'], _e['perm_p_increment_at_least_observed'],
             OR_['paired_two_factor_minus_one_factor']['point_difference'],
             OR_['paired_two_factor_minus_one_factor']['ci_lo'],
             OR_['paired_two_factor_minus_one_factor']['ci_hi']))
    if _e['fisher_exact_one_sided_p'] >= 0.05:
        fails.append('haematopoietic-home association with safety is no longer significant')
    checks += 1
    if OR_['paired_two_factor_minus_one_factor']['ci_lo'] > 0.005:
        fails.append('paired AUC gain from the expendability term is now clearly positive; '
                     'revisit the wording that it is small and bootstrap-fragile')
except FileNotFoundError:
    fails.append('outcome_test_robustness.json missing')

try:
    CV_ = _load('continuous_weights_v2.json')
    for _k in ('E1', 'E2', 'E3', 'E4'):
        check('continuous v2 scheme %s separation' % _k, CV_['schemes'][_k]['sep'])
    check('continuous v2 null p for E3 (shuffled weights)', CV_['schemes']['E3']['p_vs_shuffled_weights'], '%.2f')
    checks += 1
    if not CV_['decision']['binary_stays_primary']:
        fails.append('a continuous scheme now meets the adoption rule; revisit the binary-primary wording')
except FileNotFoundError:
    fails.append('continuous_weights_v2.json missing')

# --- third continuous attempt: measured cell-turnover weights (a PINNED NEGATIVE result) ---
try:
    ET_ = _load('external_turnover_weights_test.json')
    for _k in ('X1', 'X2', 'X3', 'X4', 'X5', 'X6'):
        check('external turnover scheme %s separation' % _k, ET_['schemes'][_k]['sep'])
    check('external turnover mapped cell types', ET_['n_types_mapped'], '%d')
    if not ET_['decision']['binary_stays_primary']:
        fails.append('a turnover-weight scheme now meets the adoption rule; revisit the binary-primary wording')
except FileNotFoundError:
    fails.append('external_turnover_weights_test.json missing')

# --- round 4: blood-flag cleanup and the search for a more effective rule (PINNED NEGATIVE/NEUTRAL results) ---
try:
    FCL_ = _load('flag_cleanup.json'); RIB_ = _load('rule_improvement.json'); RIC_ = _load('rule_improvement_c.json')
    check('flag cleanup M1 separation', FCL_['schemes']['M1']['sep'])
    check('flag cleanup M2 separation', FCL_['schemes']['M2']['sep'])
    for _k in ('C1', 'C2', 'C4', 'C5', 'C6', 'C7', 'C8'):
        check('rule search 4B %s alone' % _k, RIB_['schemes'][_k]['sep'])
    for _k in ('V1', 'V2', 'V3', 'V4', 'D1', 'D2', 'D3'):
        check('rule search 4C %s alone' % _k, RIC_['schemes'][_k]['sep'])
    if RIB_['decision']['promoted'] or RIC_['decision']['promoted']:
        fails.append('a search scheme was promoted; revisit the section 37 wording')
except FileNotFoundError:
    fails.append('round 4 result files missing')

# --- round 5: validation set, flip sensitivity, independent naming test, threshold-form rules ---
try:
    VT_ = _load('validation_test.json'); VA_ = _load('validation_addendum.json'); IT_ = _load('independent_tests.json'); TH_ = _load('threshold_harm.json')
    check('validation primary separation', VT_['primary']['sep_rule'])
    check('validation primary exact p', VT_['primary']['p_exact_permutation'])
    check('validation with absent gene as safest', VA_['sep_with_absent_gene_as_safest'])
    check('validation hedged toxic vs safe', VT_['hedged']['sep_hedged_toxic_vs_primary_safe'])
    check('validation pooled (not independent)', VT_['pooled_not_independent']['sep'])
    for _k in ('F1 NECTIN4->toxic', 'F2 F3,FOLR1->safe', 'F3 F1+F2', 'F4 TNFRSF17->toxic'):
        check('label flip %s' % _k, VT_['flip_sensitivity'][_k]['sep'])
    for _m in ('N2', 'N3', 'N4'):
        check('naming sign-test p %s' % _m, IT_['naming']['summary'][_m]['sign_test_p'])
    check('vital-share rule R1', IT_['rules']['R1']['sep']); check('vital-share rule R2', IT_['rules']['R2']['sep'])
    for _k in ('TH1(1)', 'TH2(1)', 'TH1(10)', 'TH2(10)', 'TH1(50)', 'TH2(50)'):
        check('threshold harm %s' % _k, TH_['schemes'][_k]['sep'])
    if TH_['decision']['promoted'] or IT_['decision']['R1']['adopt'] or IT_['decision']['R2']['adopt']:
        fails.append('a round 5 scheme was promoted/adopted; revisit section 38')
except FileNotFoundError:
    fails.append('round 5 result files missing')

# --- round 6 (section 39): weights diagnostic, measured FAERS severity weights, relative-exposure rules, naming test ---
try:
    WD_ = _load('weights_diagnostic.json'); FR_ = _load('faers_rules.json'); RX_ = _load('relative_exposure.json'); RN_ = _load('relative_naming.json')
    check('blood flag alone', WD_['blood_flag_only_sep']); check('frozen rule within non-blood targets', WD_['frozen_within_nonblood_sep'])
    for _k in ('F', 'G'):
        check('family weights %s in-sample best' % _k, WD_[_k]['oracle_sep']); check('family weights %s nested CV' % _k, WD_[_k]['nested_cv_sep'])
        check('family weights %s shuffled-label mean' % _k, WD_[_k]['null_oracle_mean']); check('family weights %s shuffled-label p95' % _k, WD_[_k]['null_oracle_p95'])
    for _k, _x in FR_['variants'].items():
        check('FAERS %s separation' % _k, _x['sep']); check('FAERS %s paired diff (abs)' % _k, abs(_x['diff_vs_frozen']))
        check('FAERS %s shuffled-weight p' % _k, _x['p_shuffled_weights'])
    for _k, _x in RX_['variants'].items():
        check('rel-exposure %s separation' % _k, _x['sep']); check('rel-exposure %s paired diff (abs)' % _k, abs(_x['diff_vs_frozen']))
        check('rel-exposure %s shuffled-weight p' % _k, _x['p_shuffled_weights']); check('rel-exposure %s family-wise p' % _k, _x['p_familywise'])
        check('rel-exposure %s non-blood separation' % _k, _x['nonblood_sep'])
    check('rel-exposure family-wise null p95', RX_['null_family_p95'])
    _o = RX_['diagnostic_system_weight_oracle']['sum']
    check('system-weight oracle (sum) in-sample', _o['oracle_sep']); check('system-weight oracle (sum) nested CV', _o['nested_cv_sep'])
    check('system-weight oracle (sum) shuffled-label p95', RX_['diagnostic_system_weight_oracle']['sum']['null_oracle_p95'])
    for _m in ('hand', 'death', 'severe'):
        check('naming sign-test p %s' % _m, RN_['summary'][_m]['sign_p'])
        check('naming median best rank %s' % _m, RN_['summary'][_m]['median_best_rank'], fmt='%.1f')
    check('naming F3 rank with death weights', RN_['per_target']['F3']['death'], fmt='%.1f')
    if 'now number 74' not in text or '52 alternatives had been tried at this point' not in text:
        fails.append('the count of alternatives tried (74) must be stated in sections 38.5 and 39')
    if FR_['decision']['promoted'] or RX_['decision']['promoted']:
        fails.append('a round 6 rule was promoted; revisit section 39')
except FileNotFoundError:
    fails.append('round 6 result files missing')

try:
    HY_ = _load('candidate_hygiene.json')
    checks += 1
    print('  candidate funnel: %d -> %d pass hygiene -> %d also not already pursued; '
          'false-exclusion upper bound %.0f%%' % (HY_['n_candidates'], HY_['n_pass'],
          HY_['n_pass_and_not_already_pursued'],
          100 * HY_['false_exclusion_check']['false_exclusion_rate_upper_bound']))
    if HY_['false_exclusion_check']['false_exclusion_rate_upper_bound'] > 0.10:
        fails.append('hygiene filter now falsely excludes >10% of known surface targets')
    checks += 1
    if HY_['positive_control']['fail_hygiene']:
        fails.append('a positive-control recovery fails hygiene: %s' % HY_['positive_control']['fail_hygiene'])
except FileNotFoundError:
    fails.append('candidate_hygiene.json missing')

try:
    NF = _load('clinical_units_noise_floor.json')
    checks += 1
    _s = NF['summary']['0.01']
    print('  myeloma clinical units @1%% harm: best control %.3f vs BCMA %.3f, %d/%d markers above floor'
          % (_s['best_control_ceiling'], NF['per_marker']['CD269']['0.01']['ceiling'],
             _s['n_above_floor'], _s['n_markers']))
    if _s['n_above_floor'] > 6:
        fails.append('myeloma clinical-units metric now clears the noise floor for many markers; '
                     'the Angle 2 diagnosis (floor-limited) must be revisited')
    DL = _load('clinical_units_dlbcl.json')
    checks += 1
    _x = max(v['max_abs_diff'] for v in DL['crosscheck_vs_published'].values())
    print('  DLBCL populations vs published AUC: max abs diff %.4f' % _x)
    if _x > 0.005:
        fails.append('DLBCL clinical-units populations no longer reproduce the published AUCs')
except FileNotFoundError:
    fails.append('clinical-units noise/DLBCL json missing')

try:
    A5_ = _load('angle5_monotone_premise.json')
    checks += 1
    _r = A5_['T1_on_target_one_map']
    print('  Angle 5 T1 R2: wt->hf %.3f (precommitted), wt->lz %.3f (precommitted), hf->lz %.3f (POST-HOC)'
          % (_r['wt->hf']['cv_r2'], _r['wt->lz']['cv_r2'], _r['hf->lz']['cv_r2']))
    if _r['wt->hf']['cv_r2'] >= 0.40:
        fails.append('Angle 5 precommitted WT test now clears the 0.40 line; reread the precommitment')
except FileNotFoundError:
    fails.append('angle5_monotone_premise.json missing')


try:
    LFR_ = _load('label_fragility.json'); BBL_ = _load('breadth_bound_lemma.json')
    for _k in range(1, 7):
        check('fragility k=%d lowest separation' % _k, LFR_['F1'][str(_k)]['sep_min'])
    for _k in (4, 5, 6):
        check('fragility k=%d share still p<0.05 (%%)' % _k, 100 * LFR_['F1'][str(_k)]['frac_p_lt_0.05'], fmt='%.1f')
    check('fragility k=3 worst p', LFR_['F1']['3']['worst_p']); check('fragility k=4 worst p', LFR_['F1']['4']['worst_p'])
    _s = LFR_['F2_summary']
    check('weak-label all-eight separation', _s['all_eight']['sep']); check('weak-label all-eight p', _s['all_eight']['p'])
    check('weak-label subsets min separation', _s['sep_min']); check('weak-label subsets max separation', _s['sep_max'])
    check('weak-label subsets share p<0.05 (%)', 100 * _s['frac_p_lt_0.05'], fmt='%.1f')
    check('weak-label 4 swaps max p', _s['by_n_flipped_max_p']['4'])
    for _r in ('0.03', '0.1', '0.3'):
        check('Spearman(top share, exposed types) r=%s (abs)' % _r, abs(BBL_['spearman_top_share_vs_E'][_r]), fmt='%.2f')
    check('median bound safe (r=0.1)', BBL_['safe_vs_toxic_median_bound_r0.1']['safe'], fmt='%.0f')
    check('median bound toxic (r=0.1)', BBL_['safe_vs_toxic_median_bound_r0.1']['toxic'], fmt='%.0f')
    checks += 1
    if LFR_['fragility_index'] != 4 or LFR_['k4_sets_losing_significance']['n'] != 349:
        fails.append('fragility index or count of significance-losing 4-swap sets changed; revise section 40.1')
    if BBL_['violations_panel'] or BBL_['violations_random'] or not BBL_['tightness_all_within_1']:
        fails.append('breadth lemma violated or no longer tight; revise section 40.2')
    if 'now number 74' not in text:
        fails.append('section 40 must not change the count of alternatives (74)')
except FileNotFoundError:
    fails.append('label_fragility.json or breadth_bound_lemma.json missing')

print('\n' + '=' * 78)
if fails:
    print('%d CHECKS RUN, %d FAILED' % (checks, len(fails)))
    for f in fails:
        print('  FAIL: %s' % f)
    sys.exit(1)
print('%d checks run, all passed. Document and code agree.' % checks)
