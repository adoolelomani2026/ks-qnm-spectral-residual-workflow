#!/usr/bin/env python3
"""Generate continued-fraction truncation checks for CQG revision 1."""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from qnm.leaver import (  # noqa: E402
    DEFAULT_ROOT_MAX_FUNCTION_EVALUATIONS,
    DEFAULT_ROOT_RESIDUAL_TOLERANCE,
    DEFAULT_ROOT_TOLERANCE,
    DEFAULT_TAYLOR_ORDER,
    solve_leaver_mode,
    spectral_mode_targets,
)


DEFAULT_DEPTHS = (60, 90, 120, 180, 240, 320)


@dataclass(frozen=True)
class ConvergenceRow:
    case: str
    a: float
    mode: str
    overtone: int
    spectral_n: int
    taylor_order: int
    cf_depth: int
    omega: complex
    continued_fraction_abs: float
    final_omega: complex
    difference_from_final: float
    relative_difference_from_final: float


def compute_rows(
    depths: tuple[int, ...] = DEFAULT_DEPTHS,
    taylor_order: int = DEFAULT_TAYLOR_ORDER,
    spectral_n: int = 32,
) -> list[ConvergenceRow]:
    if tuple(sorted(set(depths))) != depths:
        raise ValueError("Continued-fraction depths must be unique and increasing.")

    deformation_grid = [0.0, 0.2, 0.5, 1.0]
    targets = spectral_mode_targets(deformation_grid, spectral_n)
    cases = [
        ("schwarzschild_scalar_l2_n0", 0.0, "fundamental", 0),
        ("ks_a1_scalar_l2_n0", 1.0, "fundamental", 0),
        ("ks_a1_scalar_l2_n1", 1.0, "first_overtone", 1),
    ]

    rows: list[ConvergenceRow] = []
    for case, a_value, mode, overtone in cases:
        initial_guess = targets[a_value][mode]
        solutions: list[tuple[int, complex, float]] = []
        for depth in depths:
            omega, residual_abs = solve_leaver_mode(
                a_value,
                initial_guess,
                depth=depth,
                order=taylor_order,
            )
            solutions.append((depth, omega, residual_abs))

        final_omega = solutions[-1][1]
        for depth, omega, residual_abs in solutions:
            difference = abs(omega - final_omega)
            rows.append(
                ConvergenceRow(
                    case=case,
                    a=a_value,
                    mode=mode,
                    overtone=overtone,
                    spectral_n=spectral_n,
                    taylor_order=taylor_order,
                    cf_depth=depth,
                    omega=omega,
                    continued_fraction_abs=residual_abs,
                    final_omega=final_omega,
                    difference_from_final=difference,
                    relative_difference_from_final=difference / abs(final_omega),
                )
            )
    return rows


def write_csv(path: Path, rows: list[ConvergenceRow]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "case",
                "a_over_M",
                "mode",
                "overtone",
                "spectral_initial_guess_N",
                "taylor_order",
                "continued_fraction_depth",
                "omega_real",
                "omega_imag",
                "continued_fraction_abs",
                "final_depth_omega_real",
                "final_depth_omega_imag",
                "absolute_difference_from_final_depth",
                "relative_difference_from_final_depth",
                "root_method",
                "root_tolerance",
                "root_max_function_evaluations",
                "arithmetic",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.case,
                    row.a,
                    row.mode,
                    row.overtone,
                    row.spectral_n,
                    row.taylor_order,
                    row.cf_depth,
                    row.omega.real,
                    row.omega.imag,
                    row.continued_fraction_abs,
                    row.final_omega.real,
                    row.final_omega.imag,
                    row.difference_from_final,
                    row.relative_difference_from_final,
                    "scipy.optimize.root(method=hybr)",
                    DEFAULT_ROOT_TOLERANCE,
                    DEFAULT_ROOT_MAX_FUNCTION_EVALUATIONS,
                    "IEEE-754 double-precision complex arithmetic",
                ]
            )


def write_report(path: Path, rows: list[ConvergenceRow]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Continued-Fraction Truncation Check",
        "",
        f"- Taylor order: `{rows[0].taylor_order}`.",
        f"- Root solver: SciPy MINPACK hybrid method with tolerance `{DEFAULT_ROOT_TOLERANCE:.1e}`.",
        f"- Maximum residual-function evaluations: `{DEFAULT_ROOT_MAX_FUNCTION_EVALUATIONS}`.",
        f"- Accepted root residual threshold: `{DEFAULT_ROOT_RESIDUAL_TOLERANCE:.1e}`.",
        "- Arithmetic: IEEE-754 double precision.",
        "- Initial guesses: branch-tracked `N=32` Chebyshev eigenvalues.",
        "",
        "| case | CF depth | Re(M omega) | Im(M omega) | difference from depth 320 | |CF| |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            f"| {row.case} | {row.cf_depth} | {row.omega.real:.12f} | "
            f"{row.omega.imag:.12f} | {row.difference_from_final:.3e} | "
            f"{row.continued_fraction_abs:.3e} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT_DIR / "outputs" / "results" / "leaver_convergence.csv",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=ROOT_DIR / "outputs" / "results" / "leaver_convergence_report.md",
    )
    args = parser.parse_args()

    rows = compute_rows()
    write_csv(args.output, rows)
    write_report(args.report, rows)
    print(f"Wrote {args.output}")
    print(f"Wrote {args.report}")


if __name__ == "__main__":
    main()
