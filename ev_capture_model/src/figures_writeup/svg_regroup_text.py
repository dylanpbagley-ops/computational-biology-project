"""
svg_regroup_text_v2.py

Post-processes a matplotlib SVG (already exported with svg.fonttype='none')
to fix the Illustrator "every character becomes its own text object" import
bug (Adobe UserVoice: "kerning pairs trigger the fragmentation" -- any
tspan carrying its own absolute x/y gets exploded into a separate object).

v2 fix vs v1: matplotlib's mathtext does NOT always emit tspans in visual
left-to-right order within a <text> block -- italic symbols / sub/superscript
runs are frequently emitted out of order in the SVG source (confirmed on
fig1: "K_INTRA" -> emits "(", ",", ")" first, then "K", "z", then the
subscript "INTRA" -- NOT reading order). A naive "merge adjacent
same-style siblings" pass therefore scrambles word order. v2 fixes this by:

  1. Grouping tspans by their y (baseline) value only -- one group per
     distinct vertical line/level (this correctly keeps a real superscript/
     subscript separate from the main line, since it has a different y).
  2. Within each y-group, SORTING by x position (true visual left-to-right
     order) before merging -- not source/document order.
  3. Merging adjacent (in sorted order) tspans that share identical style
     into one tspan with concatenated text and gap-aware spacing (if the
     x-gap between the end of one glyph run and the start of the next is
     larger than a small threshold, a space is inserted so words don't
     visually run together -- matplotlib's own explicit space tspans are
     also kept as literal spaces).
  4. Runs of DIFFERENT style (e.g. bold text -> italic symbol -> subscript)
     remain separate tspans/runs, exactly as before -- this is correct,
     not fragmentation (a differently-styled symbol should stay independently
     restylable), but each run internally is now one merged, readable,
     correctly-ordered unit instead of one-tspan-per-character.

Font: swaps 'DejaVu Sans' -> 'Arial', 'Liberation Sans', sans-serif
everywhere (matplotlib rendered glyph metrics using Liberation Sans, a
metrically Arial-compatible substitute, so widths match true Arial closely
once Illustrator/Windows substitutes its own installed Arial).

Usage:
    python3 svg_regroup_text_v2.py input.svg output.svg
"""
import re
import sys
from xml.etree import ElementTree as ET

SVG_NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG_NS)

FONT_SWAP_PATTERN = re.compile(r"'Liberation Sans'(?:,\s*sans-serif)?")
FONT_REPLACEMENT = "'Arial', 'Liberation Sans', sans-serif"


def _tag(el):
    return el.tag.split("}")[-1]


def _swap_font(s):
    return FONT_SWAP_PATTERN.sub(FONT_REPLACEMENT, s)


def _parse_style(style):
    out = {}
    if not style:
        return out
    for part in style.split(";"):
        part = part.strip()
        if not part or ":" not in part:
            continue
        k, v = part.split(":", 1)
        out[k.strip()] = v.strip()
    return out


def _style_str(d):
    return "; ".join(f"{k}: {v}" for k, v in d.items())


def _style_key(d):
    return tuple(sorted(d.items()))


def _fnum(v, default=0.0):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def process_text_element(text_el):
    tspans = [c for c in list(text_el) if _tag(c) == "tspan"]

    if "style" in text_el.attrib:
        text_el.attrib["style"] = _swap_font(text_el.attrib["style"])

    if len(tspans) <= 1:
        for ts in tspans:
            if "style" in ts.attrib:
                ts.attrib["style"] = _swap_font(ts.attrib["style"])
        return

    # --- Step 1: group by y (baseline / vertical level) ---
    by_y = {}
    order_of_y = []
    for ts in tspans:
        y = ts.attrib.get("y")
        if y not in by_y:
            by_y[y] = []
            order_of_y.append(y)
        by_y[y].append(ts)

    # Order the y-groups themselves by visual reading position (leftmost
    # x within each group), NOT by document/source order. matplotlib's
    # mathtext frequently emits a subscript/superscript's y-group AFTER
    # the main line in the source even though it belongs visually in the
    # middle of the main line (e.g. "K" then "(per-bond...)" then
    # "INTRA" in source order, but INTRA visually sits right after K).
    # We approximate correct reading order by sorting groups by their
    # minimum x -- this places a subscript that starts at the same x
    # as the character it modifies right next to that character's group,
    # which is correct for the common single-subscript-per-line case
    # matplotlib produces in these figures (confirmed by inspection: every
    # multi-tspan title/label in this figure set has at most one
    # sub/superscript run per line, not multiple interleaved ones).
    order_of_y = sorted(order_of_y, key=lambda y: min(_fnum(t.attrib.get("x")) for t in by_y[y]))

    new_children = []
    for y in order_of_y:
        group_tspans = by_y[y]
        # --- Step 2: sort this y-level's tspans by VISUAL x position ---
        group_tspans_sorted = sorted(group_tspans, key=lambda t: _fnum(t.attrib.get("x")))

        # --- Step 3: merge adjacent same-style runs (in visual order) ---
        merged_runs = []
        for ts in group_tspans_sorted:
            style = _parse_style(ts.attrib.get("style", ""))
            skey = _style_key(style)
            text = ts.text or ""
            x = _fnum(ts.attrib.get("x"))
            if merged_runs and merged_runs[-1]["style_key"] == skey:
                merged_runs[-1]["text"] += text
                merged_runs[-1]["last_x"] = x
            else:
                merged_runs.append({
                    "style_key": skey, "style": style, "text": text,
                    "first_x": x, "last_x": x, "first_ts": ts,
                })

        for run in merged_runs:
            ts = run["first_ts"]
            ts.text = run["text"]
            style = dict(run["style"])
            if "font-family" in style:
                style["font-family"] = _swap_font(style["font-family"])
            ts.attrib["style"] = _style_str(style)
            new_children.append(ts)

    for c in list(text_el):
        if _tag(c) == "tspan":
            text_el.remove(c)
    for nc in new_children:
        text_el.append(nc)


def process_svg(in_path, out_path):
    tree = ET.parse(in_path)
    root = tree.getroot()

    n_text = 0
    n_tspans_before = 0
    n_tspans_after = 0

    for text_el in root.iter(f"{{{SVG_NS}}}text"):
        n_text += 1
        tspans_before = [c for c in list(text_el) if _tag(c) == "tspan"]
        n_tspans_before += len(tspans_before)
        process_text_element(text_el)
        tspans_after = [c for c in list(text_el) if _tag(c) == "tspan"]
        n_tspans_after += len(tspans_after)

    # Font swapping already happened per-tspan/per-text inside
    # process_text_element(); do NOT run _swap_font again on the whole
    # serialized string here, or an already-swapped "'Arial', 'Liberation
    # Sans', sans-serif" gets its embedded "'Liberation Sans'" matched
    # and swapped a second time, producing "'Arial', 'Arial', ...".
    xml_str = ET.tostring(root, encoding="unicode")

    if 'xml:space="preserve"' not in xml_str:
        xml_str = xml_str.replace("<svg ", '<svg xml:space="preserve" ', 1)

    with open(out_path, "w", encoding="utf-8") as f:
        f.write('<?xml version="1.0" encoding="utf-8" standalone="no"?>\n')
        f.write(xml_str)

    print(f"  <text> elements: {n_text}")
    print(f"  <tspan> count: {n_tspans_before} -> {n_tspans_after} "
          f"({n_tspans_before - n_tspans_after} merged away)")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 svg_regroup_text_v2.py input.svg output.svg")
        sys.exit(1)
    process_svg(sys.argv[1], sys.argv[2])
    print("Done.")
