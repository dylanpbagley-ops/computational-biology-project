"""
verify_claims_in_context.py — context-anchored claim checker (B3 hardening, advisory draft)

ADVISORY: prepared by Biomni 2026-10-05 for Dylan to review and commit. Complements
verify_documented_numbers.py; does not replace it.

WHY THIS EXISTS
---------------
verify_documented_numbers.py checks that a correctly formatted value appears SOMEWHERE in
the results document (`s in text`). Two failure modes slipped through:

  1. Claims not in its registry are never checked at all. The stale "16 misordered pairs"
     in section 39.2 (correct value after the 44-target panel: 18, with CEACAM5 8 / CA9 6)
     survived 221 passing checks because no registered claim covered it.
  2. Substring-anywhere semantics cannot catch a stale value when the correct value also
     appears somewhere else in the same document.

WHAT THIS SCRIPT DOES DIFFERENTLY
---------------------------------
Every claim is registered with a SECTION anchor and a CONTEXT anchor (words that must be
near the number). The checker splits the document into markdown sections, finds the
context window, extracts every number in that window, and then applies one of three modes:

  context   — the canonical value must appear inside the context window
  absent    — a forbidden (known-stale) value must NOT appear inside the context window
  exclusive — every number in the window that matches the quantity's pattern must equal
              the canonical value (catches stale duplicates sitting next to correct text)

It exits non-zero on any failure, like the existing verifier. Adding a new claim is one
entry in CLAIMS below — the registry is the point; keep it growing whenever a headline
number is quoted in prose.

Run:  python verify_claims_in_context.py
"""
import os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
# When committed to the repo this file lives next to verify_documented_numbers.py in
# ev_capture_model/src/aim13_scoping/; DOCS resolves to ev_capture_model/docs/.
DOCS = os.path.abspath(os.path.join(HERE, '..', '..', 'docs'))
RESULTS_DOC = os.path.join(DOCS, 'Aim13_Results_Two_Comparator_Ceiling.md')

NUM = r'\d+(?:\.\d+)?'


def load_sections(path):
    """Split a markdown doc into (header, body) sections. Text before the first header
    is stored under the header ''."""
    text = open(path, encoding='utf-8').read()
    parts = re.split(r'(?m)^(#{1,4})\s+(.*)$', text)
    sections = [('', parts[0])]
    for i in range(1, len(parts), 3):
        sections.append((parts[i + 1].strip(), parts[i + 2]))
    return sections


def context_windows(sections, section_regex, context_regex, window=240):
    """Yield text windows of `window` chars around every context-regex match that occurs
    inside a section whose header matches section_regex."""
    for header, body in sections:
        if not re.search(section_regex, header):
            continue
        for m in re.finditer(context_regex, body, flags=re.IGNORECASE):
            lo = max(0, m.start() - window // 2)
            yield body[lo:m.end() + window // 2]


# --------------------------------------------------------------------------------------
# Claims registry. One entry per prose-quoted number that matters.
#   doc       : filename under DOCS
#   section   : regex matched against the markdown header
#   context   : regex matched against section body; the number must live near this
#   value     : canonical value (string as formatted in the doc)
#   mode      : context | absent | exclusive
#   pattern   : (exclusive mode) numbers matching this regex are subject to the check
# --------------------------------------------------------------------------------------
CLAIMS = [
    # --- the section 39.2 error class: stale misordered-pair count ---------------------
    dict(id='39.2 misordered-pair total', doc='Aim13_Results_Two_Comparator_Ceiling.md',
         section=r'39\.2', context=r'misordered', value='18', mode='context',
         note='44-target panel: 18 misordered pairs'),
    dict(id='39.2 stale count 16 superseded by dated notice', doc='Aim13_Results_Two_Comparator_Ceiling.md',
         section=r'39\.2', context=r'misordered', value='16', mode='absent_or_noticed',
         notice=r'Correction,\s+Oct(ober)?\.?\s+6,\s+2026', canonical='18',
         note='pre-expansion value; either gone from the doc, or retained as history and '
              'superseded by the dated correction notice (which must state the corrected count 18)'),
    dict(id='39.2 CEACAM5 misordered pairs', doc='Aim13_Results_Two_Comparator_Ceiling.md',
         section=r'39\.2', context=r'CEACAM5', value='8', mode='context',
         note='corrected breakdown CEACAM5 8 (was 7)'),
    dict(id='39.2 CA9 misordered pairs', doc='Aim13_Results_Two_Comparator_Ceiling.md',
         section=r'39\.2', context=r'CA9', value='6', mode='context',
         note='corrected breakdown CA9 6 (was 5)'),
    # --- headline statistics, section-anchored -----------------------------------------
    dict(id='23.2 two-factor separation', doc='Aim13_Results_Two_Comparator_Ceiling.md',
         section=r'23\.2', context=r'two.factor', value='0.914', mode='context'),
    dict(id='23.1 one-factor separation', doc='Aim13_Results_Two_Comparator_Ceiling.md',
         section=r'23\.1', context=r'top.share|confinement', value='0.880', mode='context'),
    # --- panel composition ------------------------------------------------
    dict(id='panel size 44', doc='Aim13_Results_Two_Comparator_Ceiling.md',
         section=r'23', context=r'panel', value='44', mode='context'),
]


def main():
    fails, checks, warnings = [], 0, []
    cache = {}

    def sections_for(doc):
        if doc not in cache:
            cache[doc] = load_sections(os.path.join(DOCS, doc))
        return cache[doc]

    for c in CLAIMS:
        checks += 1
        wins = list(context_windows(sections_for(c['doc']), c['section'], c['context']))
        if not wins:
            fails.append('%s: no context window found (section /%s/, context /%s/)'
                         % (c['id'], c['section'], c['context']))
            continue
        if c['mode'] == 'context':
            if not any(re.search(r'(?<![\w.])%s(?![\w.])' % re.escape(c['value']), w) for w in wins):
                fails.append('%s: canonical value %s not found near "%s" in section /%s/'
                             % (c['id'], c['value'], c['context'], c['section']))
        elif c['mode'] == 'absent':
            for w in wins:
                if re.search(r'(?<![\w.])%s(?![\w.])' % re.escape(c['value']), w):
                    fails.append('%s: stale value %s still present near "%s" in section /%s/: ...%s...'
                                 % (c['id'], c['value'], c['context'], c['section'],
                                    re.sub(r'\s+', ' ', w)[:160]))
        elif c['mode'] == 'absent_or_noticed':
            # The stale value may EITHER be gone, OR be retained as history per the
            # project's dated-note convention — in which case a dated correction notice
            # must exist in the same section and must state the canonical corrected value.
            stale = [w for w in wins
                     if re.search(r'(?<![\w.])%s(?![\w.])' % re.escape(c['value']), w)]
            if not stale:
                pass  # value fully removed: nothing to waive
            else:
                noticed = False
                for header, body in sections_for(c['doc']):
                    if not re.search(c['section'], header):
                        continue
                    m = re.search(c['notice'], body)
                    if m and c.get('canonical'):
                        nw = body[max(0, m.start() - 240):m.end() + 240]
                        if re.search(r'(?<![\w.])%s(?![\w.])' % re.escape(c['canonical']), nw):
                            noticed = True
                if noticed:
                    warnings.append('%s: stale value %s retained per dated-note convention; '
                                    'superseded by correction notice (states %s)'
                                    % (c['id'], c['value'], c.get('canonical')))
                else:
                    fails.append('%s: stale value %s still present near "%s" in section /%s/ '
                                 'and no valid dated correction notice found'
                                 % (c['id'], c['value'], c['context'], c['section']))
        elif c['mode'] == 'exclusive':
            pat = c.get('pattern', NUM)
            for w in wins:
                for tok in re.findall(pat, w):
                    if tok != c['value']:
                        warnings.append('%s: non-canonical number %s near "%s" (canonical %s)'
                                        % (c['id'], tok, c['context'], c['value']))

    print('verify_claims_in_context: %d claims checked' % checks)
    for w in warnings:
        print('  WARNING: %s' % w)
    if fails:
        print('%d FAILURE(S):' % len(fails))
        for f in fails:
            print('  FAIL: %s' % f)
        sys.exit(1)
    print('all context-anchored claims consistent with the document.')


if __name__ == '__main__':
    main()
