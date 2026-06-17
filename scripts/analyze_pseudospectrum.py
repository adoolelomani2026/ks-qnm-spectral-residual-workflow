#!/usr/bin/env python3
"""Generate scalar KS QNM pseudospectrum diagnostics from the validated catalogue."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from qnm.pseudospectrum import (
    DEFAULT_QUANTILES,
    DEFAULT_THRESHOLDS,
    compute_pseudospectrum_analysis,
    compute_resolution_check,
    plot_pseudospectrum_contours,
    plot_pseudospectrum_sensitivity,
    plot_resolution_check,
    write_grid_csv,
    write_report,
    write_summary_csv,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--catalogue",
        type=Path,
        default=ROOT_DIR / "outputs" / "results" / "qnm_catalogue.csv",
        help="validated QNM catalogue CSV",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=ROOT_DIR / "outputs" / "results",
        help="directory for generated pseudospectrum tables and report",
    )
    parser.add_argument(
        "--figures-dir",
        type=Path,
        default=ROOT_DIR / "outputs" / "figures",
        help="directory for generated pseudospectrum figures",
    )
    parser.add_argument("--spectral-n", type=int, default=64, help="Chebyshev size for main contours")
    parser.add_argument("--grid-size", type=int, default=81, help="main contour grid size per axis")
    parser.add_argument("--half-width", type=float, default=0.025, help="half-width of local complex-frequency window")
    parser.add_argument(
        "--resolution-grid-size",
        type=int,
        default=41,
        help="grid size per axis for finite-N resolution checks",
    )
    args = parser.parse_args()

    args.results_dir.mkdir(parents=True, exist_ok=True)
    args.figures_dir.mkdir(parents=True, exist_ok=True)

    grids, summaries = compute_pseudospectrum_analysis(
        catalogue_path=args.catalogue,
        spectral_n=args.spectral_n,
        grid_size=args.grid_size,
        half_width=args.half_width,
    )
    resolution_summaries = compute_resolution_check(
        catalogue_path=args.catalogue,
        grid_size=args.resolution_grid_size,
        half_width=args.half_width,
    )

    grid_csv = args.results_dir / "scalar_l2_pseudospectrum_grid.csv"
    summary_csv = args.results_dir / "scalar_l2_pseudospectrum_summary.csv"
    resolution_csv = args.results_dir / "scalar_l2_pseudospectrum_resolution_check.csv"
    report = args.results_dir / "scalar_pseudospectrum_report.md"
    contour_figure = args.figures_dir / "scalar_l2_pseudospectrum_contours.png"
    sensitivity_figure = args.figures_dir / "scalar_l2_pseudospectrum_sensitivity.png"
    resolution_figure = args.figures_dir / "scalar_l2_pseudospectrum_resolution_check.png"

    write_grid_csv(grid_csv, grids)
    write_summary_csv(summary_csv, summaries, thresholds=DEFAULT_THRESHOLDS, quantiles=DEFAULT_QUANTILES)
    write_summary_csv(
        resolution_csv,
        resolution_summaries,
        thresholds=DEFAULT_THRESHOLDS,
        quantiles=DEFAULT_QUANTILES,
    )
    plot_pseudospectrum_contours(contour_figure, grids)
    plot_pseudospectrum_sensitivity(sensitivity_figure, summaries, thresholds=DEFAULT_THRESHOLDS)
    plot_resolution_check(resolution_figure, resolution_summaries)
    write_report(report, summaries, resolution_summaries, thresholds=DEFAULT_THRESHOLDS)

    endpoint = sorted(summaries, key=lambda row: row.a)[-1]
    schwarzschild = sorted(summaries, key=lambda row: row.a)[0]
    threshold = DEFAULT_THRESHOLDS[0]
    area_factor = endpoint.threshold_areas[threshold] / schwarzschild.threshold_areas[threshold]
    q10_gain = (-endpoint.quantiles[0.10]) - (-schwarzschild.quantiles[0.10])

    print(f"Read catalogue: {args.catalogue}")
    print(f"Wrote grid CSV: {grid_csv}")
    print(f"Wrote summary CSV: {summary_csv}")
    print(f"Wrote resolution CSV: {resolution_csv}")
    print(f"Wrote report: {report}")
    print(f"Wrote contour figure: {contour_figure}")
    print(f"Wrote sensitivity figure: {sensitivity_figure}")
    print(f"Wrote resolution figure: {resolution_figure}")
    print(f"Q10 susceptibility gain a/M=0 -> 1: {q10_gain:.3f}")
    print(f"Area factor for log10(eta)<={threshold:g}: {area_factor:.2f}")


if __name__ == "__main__":
    main()
