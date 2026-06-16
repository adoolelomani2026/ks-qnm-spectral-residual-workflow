#!/usr/bin/env python3
"""Analyze Schwarzschild-relative physics trends in the generated QNM catalogue."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from qnm.analysis import read_catalogue_csv, write_physics_analysis


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--catalogue",
        type=Path,
        default=ROOT_DIR / "outputs" / "results" / "qnm_catalogue.csv",
        help="input catalogue CSV",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=ROOT_DIR / "outputs" / "results",
        help="directory for generated analysis tables and report",
    )
    parser.add_argument(
        "--figures-dir",
        type=Path,
        default=ROOT_DIR / "outputs" / "figures",
        help="directory for generated analysis figures",
    )
    args = parser.parse_args()

    rows = read_catalogue_csv(args.catalogue)
    outputs = write_physics_analysis(args.results_dir, args.figures_dir, rows)

    print(f"Read catalogue: {args.catalogue}")
    print("Wrote catalogue physics analysis:")
    for label, path in outputs.items():
        print(f"  {label}: {path}")


if __name__ == "__main__":
    main()
