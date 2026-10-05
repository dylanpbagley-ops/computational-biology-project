"""
export_svgs_fixed.py

Corrected re-exporter for all write-up figures. The original export_svgs.py
used a hardcoded, incomplete list of "entry point" function names
('build', 'plot_reordering_example', 'plot_skew_scatter', 'plot_comparison')
and silently SKIPPED any figure module whose real entry point wasn't one of
those four -- which is exactly why some figures (fig7 in particular, plus
several others) never actually got regenerated with svg.fonttype='none' and
were still carrying matplotlib's default path-outlined text.

This version calls each figure's ACTUAL save function with its real
signature (confirmed by reading every fig*.py script directly), for both:
  - the .svg version (live, editable text; standardized font)
  - the matching .png version (same code path, run twice -- every one of
    these scripts is deterministic, confirmed by diffing repeated runs, so
    there is no drift risk between the two output files)

Font standardization: DejaVu Sans (matplotlib's own default, and the only
family every original figure was already actually using -- confirmed by
grepping the existing SVGs) is set explicitly via rcParams so it's a
deliberate, documented choice rather than an implicit default that could
silently change on a different machine/matplotlib version.

Run from ev_capture_model/src/figures_writeup/:
    python3 export_svgs_fixed.py
"""
import os
import sys
import importlib
import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"       # keep <text>, never convert to <path>
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["DejaVu Sans"]

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
FIGURES_DIR = os.path.abspath(os.path.join(HERE, "..", "..", "..", "Poster, Paper, and Write Up", "Figures"))


def _run_pair(stem, runner):
    """Call `runner(path)` once for the .svg output and once for the
    identical .png output, then report the resulting <text> vs <path>
    element counts in the SVG as a live sanity check."""
    svg_path = os.path.join(FIGURES_DIR, stem + ".svg")
    png_path = os.path.join(FIGURES_DIR, stem + ".png")
    runner(svg_path)
    runner(png_path)
    svg_content = open(svg_path, encoding="utf-8").read()
    n_text = svg_content.count("<text")
    n_path = svg_content.count("<path")
    print("  %-48s <text>: %-4d <path>: %-5d" % (stem, n_text, n_path))


def main():
    # --- Simple one-arg build/plot(save_path) figures ---
    simple_positional = {
        "fig0_graphical_abstract": "build",
        "fig1_superselectivity_heatmap_poster": "build",
        "fig5_aim12_cytof_reordering_example_poster": "build",
        "fig8_evidence_strength_comparison_poster": "build",
    }
    simple_kwarg = {
        "fig2_ca199_egfr_ranking_flip": "plot_flip_comparison",
        "fig5_aim12_cytof_reordering_example": "plot_reordering_example",
        "fig6_aim12_cytof_distribution_skew": "plot_skew_scatter",
        "fig7_conceptual_schematic": "plot_schematic",
    }

    for modname, fn_name in simple_positional.items():
        mod = importlib.import_module(modname)
        fn = getattr(mod, fn_name)
        _run_pair(modname, lambda path, fn=fn: fn(path))

    for modname, fn_name in simple_kwarg.items():
        mod = importlib.import_module(modname)
        fn = getattr(mod, fn_name)
        _run_pair(modname, lambda path, fn=fn: fn(save_path=path))

    # --- fig1: needs build_alpha_grid() first, then plot_heatmap(...) ---
    def run_fig1(path):
        mod = importlib.import_module("fig1_superselectivity_heatmap")
        k_intra_values, z_values, grid, excluded = mod.build_alpha_grid()
        mod.plot_heatmap(k_intra_values, z_values, grid, excluded, save_path=path)
    _run_pair("fig1_superselectivity_heatmap", run_fig1)

    # --- fig3: needs compute_marker_level_yields() first ---
    def run_fig3(path):
        mod = importlib.import_module("fig3_stemcell_marker_crossover")
        k_intra_values, marker_yields = mod.compute_marker_level_yields()
        mod.plot_crossover(k_intra_values, marker_yields, save_path=path)
    _run_pair("fig3_stemcell_marker_crossover", run_fig3)

    # --- fig4: needs build_status_grid() first ---
    def run_fig4(path):
        mod = importlib.import_module("fig4_aim11_lbcl_reordering")
        fit_names, thresholds, status, result = mod.build_status_grid()
        mod.plot_status_grid(fit_names, thresholds, status, result, save_path=path)
    _run_pair("fig4_aim11_lbcl_reordering", run_fig4)

    # --- fig8 (paper): needs build_comparison() first ---
    def run_fig8(path):
        mod = importlib.import_module("fig8_evidence_strength_comparison")
        bars = mod.build_comparison()
        mod.plot_comparison(bars, save_path=path)
    _run_pair("fig8_evidence_strength_comparison", run_fig8)

    # --- fig9, fig10: build(out_path=None) ---
    def run_fig9(path):
        mod = importlib.import_module("fig9_wholebody_restriction")
        mod.build(out_path=path)
    _run_pair("fig9_wholebody_restriction", run_fig9)

    def run_fig10(path):
        mod = importlib.import_module("fig10_pancreatic_thread")
        mod.build(out_path=path)
    _run_pair("fig10_pancreatic_thread", run_fig10)


if __name__ == "__main__":
    print("Regenerating all figures with svg.fonttype='none' and standardized DejaVu Sans...")
    main()
    print("Done.")
