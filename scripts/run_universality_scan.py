#!/usr/bin/env python3
"""Compare KS and Hayward scalar pseudospectral sensitivity trends."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from qnm.universality import (
    compute_barrier_metrics,
    compute_correlations,
    compute_mode_pair_diagnostics,
    plot_barrier_correlation,
    plot_model_heatmap,
    plot_softening_vs_sensitivity,
    run_universality_scan,
    summarize_model_branches,
    write_assessment,
    write_barrier_csv,
    write_dict_csv,
    write_scan_csv,
    write_verdict_csv,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=ROOT_DIR / "outputs" / "results")
    parser.add_argument("--figures-dir", type=Path, default=ROOT_DIR / "outputs" / "figures")
    parser.add_argument("--main-grid-size", type=int, default=25)
    parser.add_argument("--resolution-grid-size", type=int, default=19)
    args = parser.parse_args()

    args.results_dir.mkdir(parents=True, exist_ok=True)
    args.figures_dir.mkdir(parents=True, exist_ok=True)

    print("Running KS/Hayward universality scan...")
    scan_rows = run_universality_scan(
        main_grid_size=args.main_grid_size,
        resolution_grid_size=args.resolution_grid_size,
    )
    verdicts = summarize_model_branches(scan_rows)
    barriers = compute_barrier_metrics()
    pairs = compute_mode_pair_diagnostics(scan_rows)
    correlations = compute_correlations(verdicts, barriers)

    scan_csv = args.results_dir / "universality_model_scan.csv"
    verdict_csv = args.results_dir / "universality_branch_verdicts.csv"
    barrier_csv = args.results_dir / "universality_barrier_metrics.csv"
    pair_csv = args.results_dir / "universality_mode_pair_diagnostics.csv"
    corr_csv = args.results_dir / "universality_correlations.csv"
    report = args.results_dir / "universality_assessment.md"

    write_scan_csv(scan_csv, scan_rows)
    write_verdict_csv(verdict_csv, verdicts)
    write_barrier_csv(barrier_csv, barriers)
    write_dict_csv(pair_csv, pairs)
    write_dict_csv(corr_csv, correlations)
    write_assessment(report, verdicts, correlations, pairs)

    heatmap = args.figures_dir / "universality_endpoint_gain_heatmap.png"
    scatter = args.figures_dir / "universality_softening_vs_sensitivity.png"
    barrier = args.figures_dir / "universality_barrier_correlation.png"
    plot_model_heatmap(heatmap, verdicts)
    plot_softening_vs_sensitivity(scatter, verdicts)
    plot_barrier_correlation(barrier, verdicts, barriers)

    support = [verdict for verdict in verdicts if verdict.universality_support]
    print(f"Wrote scan CSV: {scan_csv}")
    print(f"Wrote verdict CSV: {verdict_csv}")
    print(f"Wrote report: {report}")
    print(f"Wrote heatmap: {heatmap}")
    print(f"Universality-supporting branches: {len(support)} / {len(verdicts)}")
    print("Verdict: Result is PRD/CQG quality but not PRL quality.")


if __name__ == "__main__":
    main()
