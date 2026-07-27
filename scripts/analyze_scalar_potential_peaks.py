#!/usr/bin/env python3
"""Analyze scalar-potential peaks for the CQG multipole-sensitivity revision."""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar


ROOT_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from qnm.common import A_VALUES, f_ks, horizon_radius, scalar_potential  # noqa: E402


@dataclass(frozen=True)
class PeakRow:
    ell: int
    a: float
    horizon_radius: float
    peak_radius: float
    peak_height: float
    tortoise_curvature: float
    fractional_peak_radius_shift: float
    fractional_peak_height_shift: float
    fractional_abs_curvature_shift: float
    sqrt_peak_height_shift: float


def _potential(r: float, ell: int, a: float) -> float:
    return float(scalar_potential(np.array([r]), ell, a)[0])


def _second_radial_derivative(r: float, ell: int, a: float) -> float:
    """Five-point radial second derivative at the optimized peak."""

    step = max(1.0e-5, 2.0e-4 * r)
    values = [_potential(r + offset * step, ell, a) for offset in (-2, -1, 0, 1, 2)]
    return (-values[4] + 16.0 * values[3] - 30.0 * values[2] + 16.0 * values[1] - values[0]) / (
        12.0 * step * step
    )


def compute_peak_rows(
    a_values: list[float] | None = None,
    ell_values: tuple[int, ...] = (2, 3, 4),
) -> list[PeakRow]:
    a_values = A_VALUES if a_values is None else a_values
    raw: dict[tuple[int, float], tuple[float, float, float, float]] = {}
    for ell in ell_values:
        for a_value in a_values:
            horizon = horizon_radius(a_value)
            result = minimize_scalar(
                lambda radius: -_potential(radius, ell, a_value),
                bounds=(horizon + 1.0e-6, 20.0),
                method="bounded",
                options={"xatol": 1.0e-13},
            )
            peak_radius = float(result.x)
            peak_height = _potential(peak_radius, ell, a_value)
            radial_curvature = _second_radial_derivative(peak_radius, ell, a_value)
            tortoise_curvature = float(f_ks(peak_radius, a_value) ** 2 * radial_curvature)
            raw[(ell, a_value)] = (horizon, peak_radius, peak_height, tortoise_curvature)

    rows: list[PeakRow] = []
    for ell in ell_values:
        _, reference_radius, reference_height, reference_curvature = raw[(ell, 0.0)]
        for a_value in a_values:
            horizon, peak_radius, peak_height, tortoise_curvature = raw[(ell, a_value)]
            rows.append(
                PeakRow(
                    ell=ell,
                    a=a_value,
                    horizon_radius=horizon,
                    peak_radius=peak_radius,
                    peak_height=peak_height,
                    tortoise_curvature=tortoise_curvature,
                    fractional_peak_radius_shift=(peak_radius - reference_radius) / reference_radius,
                    fractional_peak_height_shift=(peak_height - reference_height) / reference_height,
                    fractional_abs_curvature_shift=(
                        abs(tortoise_curvature) - abs(reference_curvature)
                    )
                    / abs(reference_curvature),
                    sqrt_peak_height_shift=(peak_height / reference_height) ** 0.5 - 1.0,
                )
            )
    return rows


def write_csv(path: Path, rows: list[PeakRow]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "ell",
                "a_over_M",
                "horizon_radius_over_M",
                "peak_radius_over_M",
                "M2_peak_height",
                "M4_tortoise_curvature",
                "fractional_peak_radius_shift",
                "fractional_peak_height_shift",
                "fractional_abs_tortoise_curvature_shift",
                "fractional_sqrt_peak_height_shift",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.ell,
                    row.a,
                    row.horizon_radius,
                    row.peak_radius,
                    row.peak_height,
                    row.tortoise_curvature,
                    row.fractional_peak_radius_shift,
                    row.fractional_peak_height_shift,
                    row.fractional_abs_curvature_shift,
                    row.sqrt_peak_height_shift,
                ]
            )


def plot_rows(path: Path, rows: list[PeakRow]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.7))
    for ell in (2, 3, 4):
        subset = sorted((row for row in rows if row.ell == ell), key=lambda row: row.a)
        a_values = [row.a for row in subset]
        axes[0].plot(a_values, [100.0 * row.fractional_peak_radius_shift for row in subset], "o-", label=fr"$\ell={ell}$")
        axes[1].plot(a_values, [100.0 * row.fractional_peak_height_shift for row in subset], "o-")
        axes[2].plot(a_values, [100.0 * row.fractional_abs_curvature_shift for row in subset], "o-")
    axes[0].set_ylabel(r"$\Delta r_{\rm peak}/r_{\rm peak,0}$ (\%)")
    axes[1].set_ylabel(r"$\Delta V_{\rm peak}/V_{\rm peak,0}$ (\%)")
    axes[2].set_ylabel(r"$\Delta |V''_{*,\rm peak}|/|V''_{*,\rm peak,0}|$ (\%)")
    for axis in axes:
        axis.set_xlabel(r"$a/M$")
        axis.axhline(0.0, color="black", linewidth=0.8, alpha=0.35)
        axis.grid(alpha=0.25)
    axes[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def write_report(path: Path, rows: list[PeakRow]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    endpoint_rows = [row for row in rows if row.a == 1.0]
    lines = [
        "# Scalar Potential-Peak Analysis",
        "",
        "The tortoise-coordinate curvature uses `d/dr_* = f d/dr`; at an optimized peak,",
        "`d2V/dr_*2 = f^2 d2V/dr2` because the first radial derivative vanishes.",
        "",
        "| ell | r_peak/M | peak-height shift | sqrt-height shift | |curvature| shift |",
        "| ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in endpoint_rows:
        lines.append(
            f"| {row.ell} | {row.peak_radius:.8f} | {100*row.fractional_peak_height_shift:.4f}% | "
            f"{100*row.sqrt_peak_height_shift:.4f}% | {100*row.fractional_abs_curvature_shift:.4f}% |"
        )
    lines += [
        "",
        "The `sqrt-height shift` is the leading WKB proxy for the real-frequency shift.",
        "It is used only for qualitative interpretation; the reported QNMs come from the",
        "Chebyshev--Leaver workflow, not from this approximation.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=ROOT_DIR / "outputs" / "results")
    parser.add_argument("--figures-dir", type=Path, default=ROOT_DIR / "outputs" / "figures")
    args = parser.parse_args()
    rows = compute_peak_rows()
    csv_path = args.results_dir / "scalar_potential_peak_analysis.csv"
    report_path = args.results_dir / "scalar_potential_peak_report.md"
    figure_path = args.figures_dir / "scalar_potential_peak_shifts.png"
    write_csv(csv_path, rows)
    write_report(report_path, rows)
    plot_rows(figure_path, rows)
    print(f"Wrote {csv_path}")
    print(f"Wrote {report_path}")
    print(f"Wrote {figure_path}")


if __name__ == "__main__":
    main()
