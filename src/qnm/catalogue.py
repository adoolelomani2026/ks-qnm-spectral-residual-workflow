#!/usr/bin/env python3
"""Build scalar and axial-gravitational KS QNM catalogues.

The catalogue uses the Chebyshev spectral solver as the primary algorithm and
keeps Leaver-style continued-fraction cross-validation enabled for every row.
It does not alter the spectral residual framework.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from .common import (
    A_VALUES,
    CATALOGUE_ELL_VALUES,
    CATALOGUE_OVERTONES,
    CATALOGUE_SPECTRAL_N,
    PERTURBATION_TYPES,
    SCHWARZSCHILD_REFERENCES,
    select_physical_mode,
)
from .leaver import (
    DEFAULT_CF_DEPTH,
    DEFAULT_TAYLOR_ORDER,
    solve_leaver_mode,
)
from .spectral import build_spectral_problem, generalized_eigenvalues


MODE_LABELS = {
    0: "fundamental",
    1: "first_overtone",
    2: "second_overtone",
}

LITERATURE_SOURCE = (
    "Berti-Cardoso-Starinets review/ringdown tables; "
    "rounded scalar references also agree with recent Schwarzschild tables"
)
CATALOGUE_CROSS_VALIDATION_THRESHOLD = 1.0e-4
GRAVITATIONAL_LITERATURE_TOLERANCE = 5.0e-5
SCALAR_LITERATURE_TOLERANCE = 5.0e-5


@dataclass
class CatalogueRow:
    perturbation_type: str
    ell: int
    overtone: int
    mode: str
    a: float
    spectral_n: int
    omega_spectral: complex
    omega_leaver: complex
    spectral_leaver_relative_difference: float
    continued_fraction_abs: float
    schwarzschild_literature: complex | None
    literature_relative_error: float | None
    literature_tolerance: float | None
    spectral_leaver_validation_passed: bool
    literature_validation_passed: bool | None


def reference_targets(perturbation_type: str, ell: int) -> list[complex]:
    targets: list[complex] = []
    for overtone in CATALOGUE_OVERTONES:
        key = (perturbation_type, ell, overtone)
        if key not in SCHWARZSCHILD_REFERENCES:
            raise KeyError(f"Missing Schwarzschild reference target for {key}")
        targets.append(SCHWARZSCHILD_REFERENCES[key])
    return targets


def select_modes(values: np.ndarray, targets: list[complex]) -> list[complex]:
    selected: list[complex] = []
    for target in targets:
        selected.append(select_physical_mode(values, target, exclude=selected))
    return selected


def run_catalogue(
    a_values: list[float] | None = None,
    perturbation_types: list[str] | None = None,
    ell_values: list[int] | None = None,
    spectral_n: int = CATALOGUE_SPECTRAL_N,
    leaver_depth: int = DEFAULT_CF_DEPTH,
    taylor_order: int = DEFAULT_TAYLOR_ORDER,
    validation_threshold: float = CATALOGUE_CROSS_VALIDATION_THRESHOLD,
    leaver_residual_tolerance: float = 1.0e-4,
) -> list[CatalogueRow]:
    a_values = A_VALUES if a_values is None else a_values
    perturbation_types = PERTURBATION_TYPES if perturbation_types is None else perturbation_types
    ell_values = CATALOGUE_ELL_VALUES if ell_values is None else ell_values

    rows: list[CatalogueRow] = []
    for perturbation_type in perturbation_types:
        literature_tolerance = (
            GRAVITATIONAL_LITERATURE_TOLERANCE
            if perturbation_type == "gravitational"
            else SCALAR_LITERATURE_TOLERANCE
        )
        for ell in ell_values:
            targets = reference_targets(perturbation_type, ell)
            for a in a_values:
                problem = build_spectral_problem(
                    a,
                    spectral_n,
                    ell=ell,
                    perturbation_type=perturbation_type,
                )
                values = generalized_eigenvalues(problem)
                spectral_modes = select_modes(values, targets)
                next_targets: list[complex] = []

                for overtone, omega_spectral in zip(CATALOGUE_OVERTONES, spectral_modes):
                    omega_leaver, cf_abs = solve_leaver_mode(
                        a,
                        omega_spectral,
                        depth=leaver_depth,
                        order=taylor_order,
                        ell=ell,
                        perturbation_type=perturbation_type,
                        residual_tolerance=leaver_residual_tolerance,
                    )
                    relative_difference = float(abs(omega_leaver - omega_spectral) / abs(omega_spectral))
                    literature = None
                    literature_error = None
                    literature_passed = None
                    if a == 0.0:
                        literature = SCHWARZSCHILD_REFERENCES[(perturbation_type, ell, overtone)]
                        literature_error = float(abs(omega_leaver - literature) / abs(literature))
                        literature_passed = literature_error < literature_tolerance

                    rows.append(
                        CatalogueRow(
                            perturbation_type=perturbation_type,
                            ell=ell,
                            overtone=overtone,
                            mode=MODE_LABELS[overtone],
                            a=a,
                            spectral_n=spectral_n,
                            omega_spectral=omega_spectral,
                            omega_leaver=omega_leaver,
                            spectral_leaver_relative_difference=relative_difference,
                            continued_fraction_abs=cf_abs,
                            schwarzschild_literature=literature,
                            literature_relative_error=literature_error,
                            literature_tolerance=literature_tolerance if a == 0.0 else None,
                            spectral_leaver_validation_passed=relative_difference < validation_threshold,
                            literature_validation_passed=literature_passed,
                        )
                    )
                    next_targets.append(omega_leaver)
                targets = next_targets
    return rows


def write_catalogue(output: Path, rows: list[CatalogueRow]) -> None:
    output.parent.mkdir(exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "perturbation_type",
                "ell",
                "overtone",
                "mode",
                "a_over_M",
                "spectral_N",
                "spectral_real",
                "spectral_imag",
                "leaver_real",
                "leaver_imag",
                "spectral_leaver_relative_difference",
                "continued_fraction_abs",
                "schwarzschild_literature_real",
                "schwarzschild_literature_imag",
                "literature_relative_error",
                "literature_tolerance",
                "spectral_leaver_validation_passed",
                "literature_validation_passed",
                "literature_source",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.perturbation_type,
                    row.ell,
                    row.overtone,
                    row.mode,
                    row.a,
                    row.spectral_n,
                    row.omega_spectral.real,
                    row.omega_spectral.imag,
                    row.omega_leaver.real,
                    row.omega_leaver.imag,
                    row.spectral_leaver_relative_difference,
                    row.continued_fraction_abs,
                    "" if row.schwarzschild_literature is None else row.schwarzschild_literature.real,
                    "" if row.schwarzschild_literature is None else row.schwarzschild_literature.imag,
                    "" if row.literature_relative_error is None else row.literature_relative_error,
                    "" if row.literature_tolerance is None else row.literature_tolerance,
                    row.spectral_leaver_validation_passed,
                    "" if row.literature_validation_passed is None else row.literature_validation_passed,
                    LITERATURE_SOURCE if row.schwarzschild_literature is not None else "",
                ]
            )


def write_catalogue_report(output: Path, rows: list[CatalogueRow]) -> None:
    output.parent.mkdir(exist_ok=True)
    worst_cross = max(rows, key=lambda row: row.spectral_leaver_relative_difference)
    literature_rows = [row for row in rows if row.literature_relative_error is not None]
    worst_literature = max(literature_rows, key=lambda row: row.literature_relative_error or 0.0)
    spectral_n = rows[0].spectral_n if rows else CATALOGUE_SPECTRAL_N
    lines = [
        "# QNM Catalogue Report",
        "",
        "This is the Leaver-validated catalogue. It extends the scalar Chebyshev spectral",
        "workflow to the axial gravitational Regge-Wheeler-type sector while preserving",
        "the scalar sector.",
        "",
        "## Scope",
        "",
        "- Perturbation types: scalar, gravitational.",
        "- Multipoles: ell = 2, 3, 4.",
        "- Modes: fundamental, first overtone, second overtone.",
        "- Deformations: a/M = " + ", ".join(f"{a:g}" for a in A_VALUES) + ".",
        "- Spectral comparison size: N = " + str(spectral_n) + ".",
        "",
        "For a/M > 0 the gravitational potential is the KS-lapse-deformed axial",
        "Regge-Wheeler model, `V=f_a(r)[ell(ell+1)/r^2 - 6M/r^3]`.",
        "",
        "## Validation",
        "",
        f"- Worst spectral/Leaver relative difference: `{worst_cross.spectral_leaver_relative_difference:.3e}` "
        f"({worst_cross.perturbation_type}, ell={worst_cross.ell}, n={worst_cross.overtone}, a/M={worst_cross.a:g}).",
        f"- Worst Schwarzschild literature relative error: `{(worst_literature.literature_relative_error or 0.0):.3e}` "
        f"({worst_literature.perturbation_type}, ell={worst_literature.ell}, n={worst_literature.overtone}).",
        "",
        "The automated catalogue validation fails if any spectral/Leaver relative",
        f"difference exceeds `{CATALOGUE_CROSS_VALIDATION_THRESHOLD:.1e}`. Literature checks use",
        "rounded table tolerances because several source tables report six significant figures.",
        "The Leaver solver is independent of Chebyshev collocation, matrix-pencil data,",
        "and residual minimization, but it intentionally shares the same perturbation",
        "equation, compact coordinate, endpoint factorization, and potential model.",
        "The continued-fraction residual is reported row-by-row; high-deformation second",
        "overtones can have larger CF residuals because the finite Frobenius recurrence",
        "reduction is less well conditioned there.",
    ]
    output.write_text("\n".join(lines) + "\n")


def assert_catalogue_validation(rows: list[CatalogueRow]) -> None:
    cross_failed = [row for row in rows if not row.spectral_leaver_validation_passed]
    literature_failed = [row for row in rows if row.literature_validation_passed is False]
    if cross_failed or literature_failed:
        messages = []
        if cross_failed:
            messages.append(
                "spectral/Leaver: "
                + ", ".join(
                    f"{row.perturbation_type} ell={row.ell} n={row.overtone} a={row.a:g} "
                    f"diff={row.spectral_leaver_relative_difference:.3e}"
                    for row in cross_failed
                )
            )
        if literature_failed:
            messages.append(
                "literature: "
                + ", ".join(
                    f"{row.perturbation_type} ell={row.ell} n={row.overtone} "
                    f"err={row.literature_relative_error:.3e}"
                    for row in literature_failed
                    if row.literature_relative_error is not None
                )
            )
        raise AssertionError("Catalogue validation failed: " + "; ".join(messages))


def plot_mode_trajectories(rows: list[CatalogueRow], figures_dir: Path) -> list[Path]:
    figures_dir.mkdir(exist_ok=True)
    outputs: list[Path] = []
    for perturbation_type in PERTURBATION_TYPES:
        for ell in CATALOGUE_ELL_VALUES:
            subset = [
                row
                for row in rows
                if row.perturbation_type == perturbation_type and row.ell == ell
            ]
            if not subset:
                continue
            fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3))
            for overtone in CATALOGUE_OVERTONES:
                mode_rows = sorted(
                    [row for row in subset if row.overtone == overtone],
                    key=lambda row: row.a,
                )
                a_values = np.array([row.a for row in mode_rows])
                omegas = np.array([row.omega_leaver for row in mode_rows])
                label = MODE_LABELS[overtone]
                axes[0].plot(a_values, omegas.real, marker="o", label=label)
                axes[1].plot(a_values, -omegas.imag, marker="o", label=label)
            axes[0].set_xlabel(r"$a/M$")
            axes[0].set_ylabel(r"$\mathrm{Re}(M\omega)$")
            axes[1].set_xlabel(r"$a/M$")
            axes[1].set_ylabel(r"$-\mathrm{Im}(M\omega)$")
            for axis in axes:
                axis.grid(alpha=0.25)
                axis.legend(frameon=False, fontsize=8)
            fig.tight_layout()
            output = figures_dir / f"mode_trajectories_{perturbation_type}_l{ell}.png"
            fig.savefig(output, dpi=180)
            plt.close(fig)
            outputs.append(output)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("outputs/results/qnm_catalogue.csv"))
    parser.add_argument("--figures-dir", type=Path, default=Path("outputs/figures"))
    parser.add_argument("--report", type=Path, default=Path("outputs/results/qnm_catalogue_report.md"))
    parser.add_argument("--spectral-n", type=int, default=CATALOGUE_SPECTRAL_N)
    args = parser.parse_args()

    rows = run_catalogue(spectral_n=args.spectral_n)
    write_catalogue(args.output, rows)
    write_catalogue_report(args.report, rows)
    outputs = plot_mode_trajectories(rows, args.figures_dir)
    assert_catalogue_validation(rows)

    worst = max(rows, key=lambda row: row.spectral_leaver_relative_difference)
    print(f"Wrote catalogue: {args.output}")
    print(f"Wrote report: {args.report}")
    print("Wrote trajectory figures:")
    for output in outputs:
        print(f"  {output}")
    print(
        "Worst spectral/Leaver relative difference: "
        f"{worst.spectral_leaver_relative_difference:.3e} "
        f"({worst.perturbation_type}, ell={worst.ell}, n={worst.overtone}, a/M={worst.a:g})"
    )


if __name__ == "__main__":
    main()
