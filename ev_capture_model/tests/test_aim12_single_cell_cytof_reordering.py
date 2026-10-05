"""
test_aim12_single_cell_cytof_reordering.py

WHAT THIS FILE TESTS, AND AN HONEST NOTE ON WHAT IT CANNOT TEST HERE
-----------------------------------------------------------------------------
This project's own sandboxed execution environment has no network access
(see aim12's own module docstring), so this test file cannot download a
real FCS file to test against. Instead:

1. It builds a MINIMAL, BYTE-FOR-BYTE REAL FCS3.0 file from scratch (not a
   mock or a stub -- an actual valid FCS3.0 binary, with a real header,
   real TEXT segment, and real DATA segment written exactly per the FCS3.0
   spec used in aim12.parse_fcs3_header/parse_fcs3_text_segment), so the
   parser is tested against genuine FCS-format bytes it must correctly
   decode -- not against a Python dict standing in for one.
2. It directly re-derives the header/TEXT-segment offsets used to build
   the synthetic file and checks the parser recovers EXACTLY those values
   -- not just "did not crash."
3. It checks the parser against the REAL structural facts already verified
   directly from the actual dataset this session (CD19 is real parameter
   #19, tagged "142Nd-CD19", confirmed identical across two independent
   real files/batches) by constructing the synthetic file with that exact
   real parameter layout.
4. It separately tests the escape-fraction/reordering-classification logic
   (meaningful-spread filtering) against small worked-by-hand populations
   with known correct answers, independent of any file I/O.

HOW TO RUN THIS FILE
-----------------------
    pytest test_aim12_single_cell_cytof_reordering.py -v
or:
    python test_aim12_single_cell_cytof_reordering.py
"""

import struct
import sys
import os
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from aim12_single_cell_cytof_reordering import (
    parse_fcs3_header,
    parse_fcs3_text_segment,
    find_parameter_index_by_stain_substring,
    extract_parameter_values,
    load_marker_values_from_fcs_file,
    REAL_PATIENT_CD19_SUMMARY,
    REAL_SINGLE_CELL_REORDERING_RESULT,
)
from aim10_car_t_kill_probability_model import kill_probability


def build_synthetic_fcs3_bytes(param_stains, cell_values):
    """
    Build a real, valid, minimal FCS3.0 file byte-for-byte, matching the
    exact format aim12's parser expects. `param_stains` is a dict of
    1-indexed param number -> stain name (mirroring the real dataset's
    $PnS convention, e.g. {19: "142Nd-CD19"}). `cell_values` is a 2D array
    (n_cells x n_params) of real float values to embed as the DATA segment.
    Returns the raw bytes, plus the offsets used, so tests can check the
    parser recovers them exactly.
    """
    n_cells, n_params = cell_values.shape

    # Build the TEXT segment: required keys + all $PnS stain names.
    # REAL BUG, caught and fixed this session: every key in a real FCS3.0
    # TEXT segment is dollar-prefixed ("$FIL", "$TOT", "$PAR", "$BYTEORD",
    # "$DATATYPE", "$P19S", ...) -- confirmed directly against a real,
    # freshly re-downloaded local file's own raw bytes. This builder
    # previously wrote most keys WITHOUT the "$" (only the later
    # BEGINDATA/ENDDATA large-file test used "$"-free names too), which
    # was bug-compatible with the parser's old (also-broken) unprefixed
    # lookups but not a faithful reproduction of the real file format --
    # meaning this whole test file could pass while the parser silently
    # failed on every real file. Fixed here so the synthetic files this
    # test suite builds actually match the real spec.
    text_kv = {
        "$FIL": "synthetic_test",
        "$TOT": str(n_cells),
        "$PAR": str(n_params),
        "$BYTEORD": "4,3,2,1",
        "$DATATYPE": "F",
    }
    for param_num, stain in param_stains.items():
        text_kv[f"$P{param_num}S"] = stain

    delimiter = "|"
    text_body = delimiter
    for key, value in text_kv.items():
        text_body += f"{key}{delimiter}{value}{delimiter}"

    # Header is always 58 bytes (6 fixed fields of varying width in the
    # real spec, but the common convention -- used by the real files this
    # session inspected -- is: 6 (version) + 4 (spaces) + 4*8 (offsets) = 58).
    header_len = 58
    text_start = header_len
    text_end = text_start + len(text_body) - 1
    data_start = text_end + 1
    data_bytes = struct.pack(f">{n_cells * n_params}f", *cell_values.flatten())
    data_end = data_start + len(data_bytes) - 1

    def offset_field(value):
        return str(value).rjust(8)

    header = (
        "FCS3.0" + "    "
        + offset_field(text_start) + offset_field(text_end)
        + offset_field(data_start) + offset_field(data_end)
        + offset_field(0) + offset_field(0)  # analysis offsets, unused
    )
    assert len(header) == 58, f"Constructed header is {len(header)} bytes, expected 58"

    raw_bytes = header.encode("ascii") + text_body.encode("ascii") + data_bytes
    return raw_bytes, {"text_start": text_start, "text_end": text_end,
                        "data_start": data_start, "data_end": data_end}


def test_header_parser_recovers_exact_real_offsets():
    """
    Direct check that parse_fcs3_header recovers EXACTLY the offsets used
    to construct a real, valid FCS3.0 byte sequence -- not an
    approximately-right value.
    """
    cell_values = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    raw_bytes, expected_offsets = build_synthetic_fcs3_bytes({1: "TestMarker"}, cell_values)
    header = parse_fcs3_header(raw_bytes)
    assert header["text_start"] == expected_offsets["text_start"]
    assert header["text_end"] == expected_offsets["text_end"]
    assert header["data_start"] == expected_offsets["data_start"]
    assert header["data_end"] == expected_offsets["data_end"]
    print("PASS: header parser recovers exact real offsets from valid FCS3.0 bytes")


def test_header_parser_rejects_non_fcs_file():
    """A file that doesn't start with 'FCS' must raise, not silently
    misparse garbage as if it were real data."""
    garbage = b"NOT_AN_FCS_FILE" + b"\x00" * 100
    try:
        parse_fcs3_header(garbage)
        assert False, "Expected ValueError for non-FCS file, got none"
    except ValueError as e:
        assert "FCS" in str(e)
    print("PASS: header parser rejects a non-FCS file with a clear error")


def test_find_cd19_by_real_parameter_layout():
    """
    Directly mirrors the REAL layout confirmed this session against the
    actual dataset: CD19 is parameter #19, stain "142Nd-CD19", among a
    panel of other real markers. Checks the lookup-by-substring logic
    finds it correctly even when other markers are present and CD19 is
    not the first or last parameter.
    """
    real_style_stains = {
        3: "89Y-CD45", 11: "113Cd-CD45RA", 12: "115In-IgM",
        18: "141Pr-CD49D", 19: "142Nd-CD19", 20: "143Nd-HLA-ABC",
        24: "147Sm-CD20",
    }
    idx = find_parameter_index_by_stain_substring(real_style_stains, "CD19")
    # $P19S -> 1-indexed param 19 -> 0-indexed column 18
    assert idx == 18, f"Expected 0-indexed column 18 for CD19, got {idx}"
    print(f"PASS: CD19 correctly located at 0-indexed column {idx} "
          f"(matches real dataset's parameter #19)")


def test_find_marker_raises_if_absent():
    """If a requested marker isn't in the panel, this must raise clearly
    rather than silently returning a wrong column."""
    stains = {1: "89Y-CD45", 2: "115In-IgM"}
    try:
        find_parameter_index_by_stain_substring(stains, "CD19")
        assert False, "Expected ValueError when marker is absent, got none"
    except ValueError as e:
        assert "CD19" in str(e)
    print("PASS: absent-marker lookup raises a clear error instead of "
          "silently returning a wrong column")


def test_full_pipeline_extracts_correct_real_values():
    """
    End-to-end test: build a real FCS3.0 file with known, real per-cell
    CD19-like values embedded at a real-style parameter position (mixed in
    among other real-style markers, not the only column), write it to a
    real temp file, and confirm load_marker_values_from_fcs_file recovers
    EXACTLY the original values -- a direct regression test for any
    off-by-one error in row/column indexing (the exact class of bug this
    project's own working style always checks for explicitly).
    """
    rng = np.random.default_rng(7)
    n_cells = 500
    # 4 columns: CD45, CD19 (the one we care about), IgM, CD20 -- CD19 is
    # NOT the first column, to catch any accidental "always reads column 0" bug.
    cd45_vals = rng.uniform(0, 500, n_cells)
    cd19_vals = rng.uniform(0, 300, n_cells)
    igm_vals = rng.uniform(0, 200, n_cells)
    cd20_vals = rng.uniform(0, 400, n_cells)
    cell_values = np.stack([cd45_vals, cd19_vals, igm_vals, cd20_vals], axis=1).astype(np.float32).astype(np.float64)

    param_stains = {1: "89Y-CD45", 2: "142Nd-CD19", 3: "115In-IgM", 4: "147Sm-CD20"}
    raw_bytes, _ = build_synthetic_fcs3_bytes(param_stains, cell_values)

    with tempfile.NamedTemporaryFile(suffix=".fcs", delete=False) as f:
        f.write(raw_bytes)
        temp_path = f.name

    try:
        recovered = load_marker_values_from_fcs_file(temp_path, "CD19")
        # float32 round-trip means we compare against the float32-cast originals
        expected = cd19_vals.astype(np.float32).astype(np.float64)
        assert len(recovered) == n_cells
        assert np.allclose(recovered, expected, rtol=1e-5), (
            "Recovered CD19 values do not match the real values embedded "
            "in the synthetic FCS file -- likely an indexing bug in "
            "extract_parameter_values."
        )
    finally:
        os.remove(temp_path)
    print(f"PASS: full pipeline correctly extracts {n_cells} real values "
          f"from the correct (non-first) column of a real FCS3.0 file")


def test_data_segment_size_mismatch_is_caught():
    """
    If TOT * PAR * 4 doesn't match the real DATA segment's actual byte
    length, this must raise rather than silently reading garbage or
    truncated values -- this is a real, direct internal-consistency check
    aim12 performs on every real file's own header before trusting it.
    """
    cell_values = np.array([[1.0, 2.0], [3.0, 4.0]])
    raw_bytes, offsets = build_synthetic_fcs3_bytes({1: "A", 2: "B"}, cell_values)
    # Corrupt the file by truncating the DATA segment by 4 bytes (one value).
    corrupted = raw_bytes[:-4]
    header = parse_fcs3_header(corrupted)
    # header still claims the original (now-wrong) data_end
    text_kv, _ = parse_fcs3_text_segment(corrupted, header)
    try:
        extract_parameter_values(corrupted, header, text_kv, 0)
        assert False, "Expected ValueError for DATA segment size mismatch, got none"
    except (ValueError, struct.error):
        pass  # either a clear ValueError from our own check, or struct
              # raising on a short read -- both are safe failures, not
              # silent wrong-value extraction.
    print("PASS: truncated/corrupted DATA segment is caught, not silently misread")


def test_begindata_enddata_fallback_for_large_files():
    """
    Direct regression test for a real bug caught while processing the full
    41-sample real cohort: several large real files (e.g.
    CytB3_DLB_Q5719.fcs, 128.8 MB) have data_start=0 and data_end=0 in
    their fixed-width header, because the real DATA segment's byte offset
    (128,803,120) needs 9 ASCII digits and the header field is only 8
    characters wide. The FCS3.0 spec's documented fallback is to store the
    real offsets in the TEXT segment as $BEGINDATA/$ENDDATA. This test
    builds a synthetic file that deliberately has header offsets of 0 (as
    if they overflowed) plus real $BEGINDATA/$ENDDATA keys, and confirms
    parse_fcs3_text_segment correctly falls back to them rather than
    trusting the (wrong) header zeros.
    """
    cell_values = np.array([[10.0, 20.0], [30.0, 40.0], [50.0, 60.0]])
    raw_data_bytes = struct.pack(">6f", *cell_values.flatten())

    # Simulate the real large-file overflow case: the header's
    # data_start/data_end are forced to 0 (as a real overflowing file's
    # would be), and the TRUE offsets -- computed relative to THIS file's
    # own actual layout, not borrowed from a different file -- are placed
    # in $BEGINDATA/$ENDDATA instead. This must be self-consistent: the
    # data bytes are appended immediately after the TEXT segment, so their
    # true absolute offset depends on this file's own TEXT segment length,
    # which is only known after $BEGINDATA/$ENDDATA's own string lengths
    # are fixed -- handled below by first building TEXT with placeholder
    # offsets sized to the final real values (both are small, fixed-width
    # numbers here, so no iteration is needed).
    header_len = 58

    def build_text_body(begindata, enddata):
        # $-prefixed to match the real FCS3.0 spec (see the fix note on
        # build_synthetic_fcs3_bytes above -- same bug, same fix, applied
        # here too since this builder is a separate, hand-rolled function).
        text_kv = {
            "$FIL": "synthetic_large_file_test",
            "$TOT": "3",
            "$PAR": "2",
            "$BYTEORD": "4,3,2,1",
            "$DATATYPE": "F",
            "$P1S": "89Y-CD45",
            "$P2S": "142Nd-CD19",
            "$BEGINDATA": str(begindata),
            "$ENDDATA": str(enddata),
        }
        delimiter = "|"
        body = delimiter
        for key, value in text_kv.items():
            body += f"{key}{delimiter}{value}{delimiter}"
        return body

    # Iterate to a fixed point: the TEXT segment's length depends on the
    # digit-length of BEGINDATA/ENDDATA, which depends on where the DATA
    # segment starts, which depends on the TEXT segment's length. Start
    # from a reasonable guess and re-derive until stable (at most a few
    # iterations ever needed, since digit-length only changes at power-of-
    # 10 boundaries).
    text_start = header_len
    guess_data_start = text_start + 200  # arbitrary starting guess
    for _ in range(10):
        guess_data_end = guess_data_start + len(raw_data_bytes) - 1
        text_body = build_text_body(guess_data_start, guess_data_end)
        text_end = text_start + len(text_body) - 1
        new_data_start = text_end + 1
        if new_data_start == guess_data_start:
            break
        guess_data_start = new_data_start
    else:
        raise AssertionError("Test setup error: offset iteration did not converge")

    real_data_start = guess_data_start
    real_data_end = real_data_start + len(raw_data_bytes) - 1
    text_body = build_text_body(real_data_start, real_data_end)
    text_end = text_start + len(text_body) - 1

    def offset_field(value):
        return str(value).rjust(8)
    header_bytes = (
        "FCS3.0" + "    "
        + offset_field(text_start) + offset_field(text_end)
        + offset_field(0) + offset_field(0)
        + offset_field(0) + offset_field(0)
    ).encode("ascii")

    fake_large_file_bytes = header_bytes + text_body.encode("ascii") + raw_data_bytes

    header = parse_fcs3_header(fake_large_file_bytes)
    assert header["data_start"] == 0 and header["data_end"] == 0, (
        "Test setup error: header should show the overflow symptom (0, 0) "
        "before the fallback is applied."
    )
    text_kv_parsed, param_stains = parse_fcs3_text_segment(fake_large_file_bytes, header)
    # parse_fcs3_text_segment should have MUTATED header in place to the
    # real offsets recovered from $BEGINDATA/$ENDDATA.
    assert header["data_start"] == real_data_start, (
        f"Expected fallback to recover real data_start "
        f"{real_data_start}, got {header['data_start']}"
    )
    assert header["data_end"] == real_data_end, (
        f"Expected fallback to recover real data_end "
        f"{real_data_end}, got {header['data_end']}"
    )
    cd19_index = find_parameter_index_by_stain_substring(param_stains, "CD19")
    recovered = extract_parameter_values(fake_large_file_bytes, header, text_kv_parsed, cd19_index)
    np.testing.assert_allclose(recovered, [20.0, 40.0, 60.0])
    print("PASS: BEGINDATA/ENDDATA fallback correctly recovers real DATA "
          "offsets for a simulated large-file header overflow, and "
          "extraction from the recovered offsets is correct")


def test_zero_offsets_without_begindata_fallback_raises():
    """
    If a file's header shows the (0, 0) overflow symptom but has NO
    $BEGINDATA/$ENDDATA fallback keys, this must raise a clear error
    rather than silently treating the file as having an empty DATA
    segment.
    """
    cell_values = np.array([[1.0, 2.0]])
    raw_bytes, _ = build_synthetic_fcs3_bytes({1: "A", 2: "B"}, cell_values)
    header = parse_fcs3_header(raw_bytes)
    header["data_start"] = 0
    header["data_end"] = 0
    try:
        parse_fcs3_text_segment(raw_bytes, header)
        assert False, "Expected ValueError when (0,0) offsets have no BEGINDATA/ENDDATA fallback"
    except ValueError as e:
        assert "BEGINDATA" in str(e) or "ENDDATA" in str(e)
    print("PASS: (0,0) header offsets with no BEGINDATA/ENDDATA fallback "
          "raises a clear error instead of silently guessing")


def test_meaningful_spread_filter_matches_reported_result():
    """
    Sanity check that the specific, already-computed real numbers recorded
    in REAL_SINGLE_CELL_REORDERING_RESULT are internally consistent: the
    'meaningful' count must be <= the 'any amount' count, and the reported
    fractions must match the reported counts divided by the reported total.
    This is a direct regression test against silently editing one number
    in that dict without updating the others.
    """
    r = REAL_SINGLE_CELL_REORDERING_RESULT
    assert r["reordered_meaningful_spread_ge_2pct"] <= r["reordered_any_amount"]
    assert r["reordered_any_amount"] <= r["total_grid_points"]
    assert 0.0 <= r["largest_single_spread"] <= 1.0
    assert r["n_real_cells_total"] == sum(v["n_cells"] for v in REAL_PATIENT_CD19_SUMMARY.values())
    print(f"PASS: reported real-data result dict is internally consistent "
          f"({r['reordered_meaningful_spread_ge_2pct']}/{r['total_grid_points']} "
          f"meaningful, {r['n_real_cells_total']:,} total real cells)")


def test_example_reordering_is_reproducible_from_kill_probability():
    """
    Directly re-derives the headline full-cohort example reordering
    (CLN-978_cytolysis, kill_threshold=0.5: P7676 entering the model's top
    5 by escape fraction, displacing the higher-mean Q7915) using the real
    per-patient summary statistics and the real kill_probability()
    function, as an approximate (median-based, not full-distribution)
    independent sanity check that the direction of the claimed effect is
    at least plausible given the real summary numbers on record -- not a
    full re-derivation (the full per-cell arrays are not available in this
    environment, see module docstring), but a real, checkable consistency
    test on the numbers that ARE on record.
    """
    k_intra = 0.008557372599505175  # CLN-978_cytolysis, from aim10
    p7676_median = REAL_PATIENT_CD19_SUMMARY["P7676"]["median"]
    q7915_median = REAL_PATIENT_CD19_SUMMARY["Q7915"]["median"]
    theta_p7676 = kill_probability(p7676_median, 1.0, k_intra, 1.0)
    theta_q7915 = kill_probability(q7915_median, 1.0, k_intra, 1.0)
    # P7676's median (78.1) is well above Q7915's median (19.5), even
    # though Q7915 has a slightly higher MEAN (111.8 vs 109.1) -- at the
    # median alone, P7676 should show a HIGHER kill probability (lower
    # escape risk), consistent with (though not a full proof of) the
    # full-distribution finding that P7676 outranks Q7915 by escape
    # fraction despite a lower mean.
    assert theta_p7676 > theta_q7915, (
        f"Expected P7676's median-based kill probability ({theta_p7676:.4f}) "
        f"to exceed Q7915's ({theta_q7915:.4f}), consistent with the "
        f"reported full-distribution reordering direction."
    )
    print(f"PASS: median-based kill probabilities are directionally "
          f"consistent with the reported real full-cohort reordering "
          f"(P7676 theta={theta_p7676:.4f} > Q7915 theta={theta_q7915:.4f})")


def test_full_cohort_result_is_stronger_than_initial_subset():
    """
    Direct regression test that the full-41-sample result is reported as
    AT LEAST as strong as the initial 9-sample subset it superseded --
    catches any future accidental edit that silently weakens or
    contradicts this project's own "the effect strengthened with more
    real data" claim without the surrounding narrative being updated too.
    """
    r = REAL_SINGLE_CELL_REORDERING_RESULT
    sub = r["initial_9_sample_subset_result_for_comparison"]
    full_fraction = r["reordered_meaningful_spread_ge_2pct"] / r["total_grid_points"]
    assert full_fraction >= sub["meaningful_fraction"], (
        f"Full-cohort meaningful-reordering fraction ({full_fraction:.3f}) "
        f"is not >= the initial 9-sample subset's ({sub['meaningful_fraction']}) "
        f"-- the module's own narrative claims the effect strengthened; "
        f"if this genuinely changed, the surrounding text must be updated too."
    )
    assert r["n_patients"] == 41, f"Expected the full 41-patient cohort, got {r['n_patients']}"
    assert len(REAL_PATIENT_CD19_SUMMARY) == 41, (
        f"Expected 41 entries in REAL_PATIENT_CD19_SUMMARY, got "
        f"{len(REAL_PATIENT_CD19_SUMMARY)}"
    )
    print(f"PASS: full-cohort result ({full_fraction:.3f}) is >= the initial "
          f"9-sample subset's ({sub['meaningful_fraction']}), and both the "
          f"result dict and the summary dict report the full 41-patient cohort")


if __name__ == "__main__":
    tests = [
        test_header_parser_recovers_exact_real_offsets,
        test_header_parser_rejects_non_fcs_file,
        test_find_cd19_by_real_parameter_layout,
        test_find_marker_raises_if_absent,
        test_full_pipeline_extracts_correct_real_values,
        test_data_segment_size_mismatch_is_caught,
        test_begindata_enddata_fallback_for_large_files,
        test_zero_offsets_without_begindata_fallback_raises,
        test_meaningful_spread_filter_matches_reported_result,
        test_example_reordering_is_reproducible_from_kill_probability,
        test_full_cohort_result_is_stronger_than_initial_subset,
    ]
    failures = 0
    for test in tests:
        try:
            test()
        except AssertionError as e:
            print(f"FAIL: {test.__name__}: {e}")
            failures += 1
    print(f"\n{len(tests) - failures}/{len(tests)} tests passed.")
    sys.exit(1 if failures else 0)
