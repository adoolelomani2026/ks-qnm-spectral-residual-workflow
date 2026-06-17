#!/usr/bin/env python3
"""Catalogue-level physics diagnostics for KS QNM spectra."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from .catalogue import CatalogueRow


@dataclass(frozen=True)
class PhysicsTrendRow:
    perturbation_type: str
    ell: int
    overtone: int
    mode: str
    a: float
    omega: complex
    omega_schwarzschild: complex
    damping: float
    schwarzschild_damping: float
    delta_real: float
    fractional_real_shift: float
    delta_damping: float
    fractional_damping_shift: float
    quality_factor: float
    fractional_quality_shift: float
    complex_shift_abs: float
    fractional_complex_shift: float
    spectral_leaver_relative_difference: float


@dataclass(frozen=True)
class BranchSummary:
    perturbation_type: str
    ell: int
    overtone: int
    mode: str
    endpoint_a: float
    endpoint_real: float
    endpoint_imag: float
    endpoint_fractional_real_shift: float
    endpoint_fractional_damping_shift: float
    endpoint_fractional_quality_shift: float
    endpoint_fractional_complex_shift: float
    real_trend: str
    damping_trend: str
    max_spectral_leaver_relative_difference: float
    max_continued_fraction_abs: float
    sensitivity_rank: int = 0


@dataclass(frozen=True)
class SpectroscopicRatioRow:
    perturbation_type: str
    ell: int
    a: float
    ratio_name: str
    numerator_overtone: int
    denominator_overtone: int | None
    value: complex
    schwarzschild_value: complex
    fractional_shift_abs: float
    fractional_real_shift: float
    validation_max: float


def _branch_key(row: CatalogueRow) -> tuple[str, int, int]:
    return (row.perturbation_type, row.ell, row.overtone)


def _trend(values: list[float], tolerance: float = 1.0e-12) -> str:
    diffs = np.diff(np.array(values, dtype=float))
    if np.all(diffs <= tolerance):
        return "decreasing"
    if np.all(diffs >= -tolerance):
        return "increasing"
    return "mixed"


def compute_trend_rows(rows: list[CatalogueRow]) -> list[PhysicsTrendRow]:
    grouped: dict[tuple[str, int, int], list[CatalogueRow]] = {}
    for row in rows:
        grouped.setdefault(_branch_key(row), []).append(row)

    trend_rows: list[PhysicsTrendRow] = []
    for key_rows in grouped.values():
        ordered = sorted(key_rows, key=lambda row: row.a)
        reference = ordered[0]
        omega0 = reference.omega_leaver
        damping0 = -omega0.imag
        quality0 = omega0.real / (2.0 * damping0)
        for row in ordered:
            omega = row.omega_leaver
            damping = -omega.imag
            quality = omega.real / (2.0 * damping)
            delta_real = omega.real - omega0.real
            delta_damping = damping - damping0
            complex_shift_abs = abs(omega - omega0)
            trend_rows.append(
                PhysicsTrendRow(
                    perturbation_type=row.perturbation_type,
                    ell=row.ell,
                    overtone=row.overtone,
                    mode=row.mode,
                    a=row.a,
                    omega=omega,
                    omega_schwarzschild=omega0,
                    damping=damping,
                    schwarzschild_damping=damping0,
                    delta_real=delta_real,
                    fractional_real_shift=delta_real / omega0.real,
                    delta_damping=delta_damping,
                    fractional_damping_shift=delta_damping / damping0,
                    quality_factor=quality,
                    fractional_quality_shift=(quality - quality0) / quality0,
                    complex_shift_abs=complex_shift_abs,
                    fractional_complex_shift=complex_shift_abs / abs(omega0),
                    spectral_leaver_relative_difference=row.spectral_leaver_relative_difference,
                )
            )
    return sorted(trend_rows, key=lambda row: (row.perturbation_type, row.ell, row.overtone, row.a))


def compute_branch_summaries(rows: list[CatalogueRow]) -> list[BranchSummary]:
    grouped: dict[tuple[str, int, int], list[CatalogueRow]] = {}
    for row in rows:
        grouped.setdefault(_branch_key(row), []).append(row)

    summaries: list[BranchSummary] = []
    for key_rows in grouped.values():
        ordered = sorted(key_rows, key=lambda row: row.a)
        reference = ordered[0].omega_leaver
        endpoint = ordered[-1]
        endpoint_omega = endpoint.omega_leaver
        damping0 = -reference.imag
        endpoint_damping = -endpoint_omega.imag
        quality0 = reference.real / (2.0 * damping0)
        endpoint_quality = endpoint_omega.real / (2.0 * endpoint_damping)
        summaries.append(
            BranchSummary(
                perturbation_type=endpoint.perturbation_type,
                ell=endpoint.ell,
                overtone=endpoint.overtone,
                mode=endpoint.mode,
                endpoint_a=endpoint.a,
                endpoint_real=endpoint_omega.real,
                endpoint_imag=endpoint_omega.imag,
                endpoint_fractional_real_shift=(endpoint_omega.real - reference.real) / reference.real,
                endpoint_fractional_damping_shift=(endpoint_damping - damping0) / damping0,
                endpoint_fractional_quality_shift=(endpoint_quality - quality0) / quality0,
                endpoint_fractional_complex_shift=abs(endpoint_omega - reference) / abs(reference),
                real_trend=_trend([row.omega_leaver.real for row in ordered]),
                damping_trend=_trend([-row.omega_leaver.imag for row in ordered]),
                max_spectral_leaver_relative_difference=max(
                    row.spectral_leaver_relative_difference for row in ordered
                ),
                max_continued_fraction_abs=max(row.continued_fraction_abs for row in ordered),
            )
        )

    ranked = sorted(
        summaries,
        key=lambda row: row.endpoint_fractional_complex_shift,
        reverse=True,
    )
    with_ranks = [
        BranchSummary(
            perturbation_type=row.perturbation_type,
            ell=row.ell,
            overtone=row.overtone,
            mode=row.mode,
            endpoint_a=row.endpoint_a,
            endpoint_real=row.endpoint_real,
            endpoint_imag=row.endpoint_imag,
            endpoint_fractional_real_shift=row.endpoint_fractional_real_shift,
            endpoint_fractional_damping_shift=row.endpoint_fractional_damping_shift,
            endpoint_fractional_quality_shift=row.endpoint_fractional_quality_shift,
            endpoint_fractional_complex_shift=row.endpoint_fractional_complex_shift,
            real_trend=row.real_trend,
            damping_trend=row.damping_trend,
            max_spectral_leaver_relative_difference=row.max_spectral_leaver_relative_difference,
            max_continued_fraction_abs=row.max_continued_fraction_abs,
            sensitivity_rank=rank,
        )
        for rank, row in enumerate(ranked, start=1)
    ]
    return sorted(with_ranks, key=lambda row: (row.perturbation_type, row.ell, row.overtone))


def compute_spectroscopic_ratios(rows: list[CatalogueRow]) -> list[SpectroscopicRatioRow]:
    grouped: dict[tuple[str, int, float], dict[int, CatalogueRow]] = {}
    for row in rows:
        grouped.setdefault((row.perturbation_type, row.ell, row.a), {})[row.overtone] = row

    reference_by_branch: dict[tuple[str, int], dict[int, CatalogueRow]] = {}
    for (perturbation_type, ell, a_value), overtone_rows in grouped.items():
        if abs(a_value) <= 1.0e-14:
            reference_by_branch[(perturbation_type, ell)] = overtone_rows

    ratio_rows: list[SpectroscopicRatioRow] = []
    for (perturbation_type, ell, a_value), overtone_rows in sorted(grouped.items()):
        reference_rows = reference_by_branch.get((perturbation_type, ell))
        if reference_rows is None:
            continue

        if 0 in overtone_rows and 0 in reference_rows:
            omega0 = overtone_rows[0].omega_leaver
            reference_omega0 = reference_rows[0].omega_leaver
            for overtone in (1, 2):
                if overtone not in overtone_rows or overtone not in reference_rows:
                    continue
                value = omega0 / overtone_rows[overtone].omega_leaver
                reference_value = reference_omega0 / reference_rows[overtone].omega_leaver
                validation_max = max(
                    overtone_rows[0].spectral_leaver_relative_difference,
                    overtone_rows[overtone].spectral_leaver_relative_difference,
                )
                ratio_rows.append(
                    SpectroscopicRatioRow(
                        perturbation_type=perturbation_type,
                        ell=ell,
                        a=a_value,
                        ratio_name=f"omega0_over_omega{overtone}",
                        numerator_overtone=0,
                        denominator_overtone=overtone,
                        value=value,
                        schwarzschild_value=reference_value,
                        fractional_shift_abs=abs(value - reference_value) / abs(reference_value),
                        fractional_real_shift=(value.real - reference_value.real) / reference_value.real,
                        validation_max=validation_max,
                    )
                )

        for overtone, row in sorted(overtone_rows.items()):
            reference = reference_rows.get(overtone)
            if reference is None:
                continue
            value = complex(row.omega_leaver.real / (-row.omega_leaver.imag), 0.0)
            reference_value = complex(reference.omega_leaver.real / (-reference.omega_leaver.imag), 0.0)
            ratio_rows.append(
                SpectroscopicRatioRow(
                    perturbation_type=perturbation_type,
                    ell=ell,
                    a=a_value,
                    ratio_name=f"real_to_damping_n{overtone}",
                    numerator_overtone=overtone,
                    denominator_overtone=None,
                    value=value,
                    schwarzschild_value=reference_value,
                    fractional_shift_abs=abs(value - reference_value) / abs(reference_value),
                    fractional_real_shift=(value.real - reference_value.real) / reference_value.real,
                    validation_max=row.spectral_leaver_relative_difference,
                )
            )

    return sorted(
        ratio_rows,
        key=lambda row: (row.perturbation_type, row.ell, row.ratio_name, row.a),
    )


def write_trend_csv(output: Path, trend_rows: list[PhysicsTrendRow]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "perturbation_type",
                "ell",
                "overtone",
                "mode",
                "a_over_M",
                "leaver_real",
                "leaver_imag",
                "schwarzschild_real",
                "schwarzschild_imag",
                "delta_real",
                "fractional_real_shift",
                "damping",
                "schwarzschild_damping",
                "delta_damping",
                "fractional_damping_shift",
                "quality_factor",
                "fractional_quality_shift",
                "complex_shift_abs",
                "fractional_complex_shift",
                "spectral_leaver_relative_difference",
            ]
        )
        for row in trend_rows:
            writer.writerow(
                [
                    row.perturbation_type,
                    row.ell,
                    row.overtone,
                    row.mode,
                    row.a,
                    row.omega.real,
                    row.omega.imag,
                    row.omega_schwarzschild.real,
                    row.omega_schwarzschild.imag,
                    row.delta_real,
                    row.fractional_real_shift,
                    row.damping,
                    row.schwarzschild_damping,
                    row.delta_damping,
                    row.fractional_damping_shift,
                    row.quality_factor,
                    row.fractional_quality_shift,
                    row.complex_shift_abs,
                    row.fractional_complex_shift,
                    row.spectral_leaver_relative_difference,
                ]
            )


def write_summary_csv(output: Path, summaries: list[BranchSummary]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "sensitivity_rank",
                "perturbation_type",
                "ell",
                "overtone",
                "mode",
                "endpoint_a_over_M",
                "endpoint_real",
                "endpoint_imag",
                "endpoint_fractional_real_shift",
                "endpoint_fractional_damping_shift",
                "endpoint_fractional_quality_shift",
                "endpoint_fractional_complex_shift",
                "real_trend",
                "damping_trend",
                "max_spectral_leaver_relative_difference",
                "max_continued_fraction_abs",
            ]
        )
        for row in sorted(summaries, key=lambda item: item.sensitivity_rank):
            writer.writerow(
                [
                    row.sensitivity_rank,
                    row.perturbation_type,
                    row.ell,
                    row.overtone,
                    row.mode,
                    row.endpoint_a,
                    row.endpoint_real,
                    row.endpoint_imag,
                    row.endpoint_fractional_real_shift,
                    row.endpoint_fractional_damping_shift,
                    row.endpoint_fractional_quality_shift,
                    row.endpoint_fractional_complex_shift,
                    row.real_trend,
                    row.damping_trend,
                    row.max_spectral_leaver_relative_difference,
                    row.max_continued_fraction_abs,
                ]
            )


def write_spectroscopic_ratio_csv(output: Path, ratio_rows: list[SpectroscopicRatioRow]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "perturbation_type",
                "ell",
                "a_over_M",
                "ratio_name",
                "numerator_overtone",
                "denominator_overtone",
                "value_real",
                "value_imag",
                "schwarzschild_value_real",
                "schwarzschild_value_imag",
                "fractional_shift_abs",
                "fractional_real_shift",
                "validation_max",
            ]
        )
        for row in ratio_rows:
            writer.writerow(
                [
                    row.perturbation_type,
                    row.ell,
                    row.a,
                    row.ratio_name,
                    row.numerator_overtone,
                    "" if row.denominator_overtone is None else row.denominator_overtone,
                    row.value.real,
                    row.value.imag,
                    row.schwarzschild_value.real,
                    row.schwarzschild_value.imag,
                    row.fractional_shift_abs,
                    row.fractional_real_shift,
                    row.validation_max,
                ]
            )


def plot_l2_spectroscopic_ratios(ratio_rows: list[SpectroscopicRatioRow], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    selected_names = ["omega0_over_omega1", "omega0_over_omega2", "real_to_damping_n0"]
    labels = {
        "omega0_over_omega1": r"$|\Delta(\omega_0/\omega_1)|$",
        "omega0_over_omega2": r"$|\Delta(\omega_0/\omega_2)|$",
        "real_to_damping_n0": r"$\Delta[\mathrm{Re}(\omega_0)/-\mathrm{Im}(\omega_0)]$",
    }
    colors = {
        "omega0_over_omega1": "#1f5f8b",
        "omega0_over_omega2": "#d1495b",
        "real_to_damping_n0": "#2a9d8f",
    }

    fig, axis = plt.subplots(figsize=(7.8, 4.2))
    for ratio_name in selected_names:
        rows = sorted(
            [
                row
                for row in ratio_rows
                if row.perturbation_type == "scalar" and row.ell == 2 and row.ratio_name == ratio_name
            ],
            key=lambda row: row.a,
        )
        if not rows:
            continue
        if ratio_name.startswith("real_to_damping"):
            values = [100.0 * row.fractional_real_shift for row in rows]
        else:
            values = [100.0 * row.fractional_shift_abs for row in rows]
        axis.plot(
            [row.a for row in rows],
            values,
            marker="o",
            linewidth=2.0,
            color=colors[ratio_name],
            label=labels[ratio_name],
        )

    axis.axhline(0.0, color="black", linewidth=0.8, alpha=0.35)
    axis.set_xlabel(r"$a/M$")
    axis.set_ylabel("Schwarzschild-relative ratio shift (%)")
    axis.grid(alpha=0.25)
    axis.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_l2_fractional_shifts(trend_rows: list[PhysicsTrendRow], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    subset = [row for row in trend_rows if row.ell == 2]
    colors = {0: "#1f5f8b", 1: "#2a9d8f", 2: "#d1495b"}
    linestyles = {"scalar": "-", "gravitational": "--"}
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.2))
    for perturbation_type in sorted({row.perturbation_type for row in subset}):
        for overtone in sorted({row.overtone for row in subset}):
            rows = sorted(
                [
                    row
                    for row in subset
                    if row.perturbation_type == perturbation_type and row.overtone == overtone
                ],
                key=lambda row: row.a,
            )
            if not rows:
                continue
            label = f"{perturbation_type}, n={overtone}"
            a_values = np.array([row.a for row in rows])
            axes[0].plot(
                a_values,
                [100.0 * row.fractional_real_shift for row in rows],
                marker="o",
                color=colors[overtone],
                linestyle=linestyles[perturbation_type],
                label=label,
            )
            axes[1].plot(
                a_values,
                [100.0 * row.fractional_damping_shift for row in rows],
                marker="o",
                color=colors[overtone],
                linestyle=linestyles[perturbation_type],
                label=label,
            )
    axes[0].set_xlabel(r"$a/M$")
    axes[0].set_ylabel(r"$\Delta \mathrm{Re}(\omega) / \mathrm{Re}(\omega_0)$ (%)")
    axes[1].set_xlabel(r"$a/M$")
    axes[1].set_ylabel(r"$\Delta[-\mathrm{Im}(\omega)] / [-\mathrm{Im}(\omega_0)]$ (%)")
    for axis in axes:
        axis.axhline(0.0, color="black", linewidth=0.8, alpha=0.35)
        axis.grid(alpha=0.25)
        axis.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def _find_ratio(
    ratio_rows: list[SpectroscopicRatioRow],
    perturbation_type: str,
    ell: int,
    ratio_name: str,
    a_value: float,
) -> SpectroscopicRatioRow:
    return next(
        row
        for row in ratio_rows
        if row.perturbation_type == perturbation_type
        and row.ell == ell
        and row.ratio_name == ratio_name
        and abs(row.a - a_value) <= 1.0e-12
    )


def write_physics_report(
    output: Path,
    summaries: list[BranchSummary],
    ratio_rows: list[SpectroscopicRatioRow],
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    ranked = sorted(summaries, key=lambda row: row.sensitivity_rank)
    top = ranked[0]
    real_trends = {row.real_trend for row in summaries}
    damping_trends = {row.damping_trend for row in summaries}
    worst_validation = max(summaries, key=lambda row: row.max_spectral_leaver_relative_difference)
    scalar_l2_fund = next(
        row for row in summaries if row.perturbation_type == "scalar" and row.ell == 2 and row.overtone == 0
    )
    scalar_l2_second = next(
        row for row in summaries if row.perturbation_type == "scalar" and row.ell == 2 and row.overtone == 2
    )
    grav_l2_fund = next(
        row
        for row in summaries
        if row.perturbation_type == "gravitational" and row.ell == 2 and row.overtone == 0
    )
    scalar_l2_omega01 = _find_ratio(ratio_rows, "scalar", 2, "omega0_over_omega1", 1.0)
    scalar_l2_omega02 = _find_ratio(ratio_rows, "scalar", 2, "omega0_over_omega2", 1.0)
    scalar_l2_re_damp_0 = _find_ratio(ratio_rows, "scalar", 2, "real_to_damping_n0", 1.0)
    scalar_l2_re_damp_1 = _find_ratio(ratio_rows, "scalar", 2, "real_to_damping_n1", 1.0)
    scalar_l2_re_damp_2 = _find_ratio(ratio_rows, "scalar", 2, "real_to_damping_n2", 1.0)

    lines = [
        "# Catalogue Physics Analysis",
        "",
        "This report analyzes the Leaver-validated catalogue by comparing each KS",
        "branch against its Schwarzschild endpoint at a/M=0.",
        "",
        "## Main Observations",
        "",
        f"- Largest endpoint fractional frequency shift: `{100.0 * top.endpoint_fractional_complex_shift:.2f}%` "
        f"for `{top.perturbation_type} ell={top.ell} n={top.overtone}` at a/M={top.endpoint_a:g}.",
        f"- Scalar ell=2 fundamental endpoint shift: real part `{100.0 * scalar_l2_fund.endpoint_fractional_real_shift:.2f}%`, "
        f"damping magnitude `{100.0 * scalar_l2_fund.endpoint_fractional_damping_shift:.2f}%`.",
        f"- Scalar ell=2 fundamental quality-factor shift: `{100.0 * scalar_l2_fund.endpoint_fractional_quality_shift:.2f}%`.",
        f"- Scalar ell=2 spectroscopic ratio shift `omega0/omega1`: "
        f"`{100.0 * scalar_l2_omega01.fractional_shift_abs:.2f}%`.",
        f"- Scalar ell=2 spectroscopic ratio shift `omega0/omega2`: "
        f"`{100.0 * scalar_l2_omega02.fractional_shift_abs:.2f}%`.",
        f"- Scalar ell=2 `Re(omega)/[-Im(omega)]` shifts: n=0 "
        f"`{100.0 * scalar_l2_re_damp_0.fractional_real_shift:.2f}%`, "
        f"n=1 `{100.0 * scalar_l2_re_damp_1.fractional_real_shift:.2f}%`, "
        f"n=2 `{100.0 * scalar_l2_re_damp_2.fractional_real_shift:.2f}%`.",
        f"- Scalar ell=2 second overtone endpoint shift: real part `{100.0 * scalar_l2_second.endpoint_fractional_real_shift:.2f}%`, "
        f"damping magnitude `{100.0 * scalar_l2_second.endpoint_fractional_damping_shift:.2f}%`.",
        f"- Axial gravitational ell=2 fundamental endpoint shift: real part `{100.0 * grav_l2_fund.endpoint_fractional_real_shift:.2f}%`, "
        f"damping magnitude `{100.0 * grav_l2_fund.endpoint_fractional_damping_shift:.2f}%`.",
        f"- Real-part trends across the grid: `{', '.join(sorted(real_trends))}`.",
        f"- Damping-magnitude trends across the grid: `{', '.join(sorted(damping_trends))}`.",
        f"- Worst spectral/Leaver mismatch remains `{worst_validation.max_spectral_leaver_relative_difference:.3e}` "
        f"for `{worst_validation.perturbation_type} ell={worst_validation.ell} n={worst_validation.overtone}`.",
        "",
        "## Dimensionless Spectroscopic Ratios",
        "",
        "The ratio table records mode ratios such as `omega0/omega1`,",
        "`omega0/omega2`, and `Re(omega)/[-Im(omega)]`. These diagnostics are",
        "less directly absorbable into a simple mass rescaling than individual",
        "frequencies.",
        "",
        "| branch | ratio | endpoint value | Schwarzschild value | shift |",
        "|---|---|---:|---:|---:|",
        f"| scalar ell=2 | omega0/omega1 | "
        f"{scalar_l2_omega01.value.real:.6f}{scalar_l2_omega01.value.imag:+.6f}i | "
        f"{scalar_l2_omega01.schwarzschild_value.real:.6f}{scalar_l2_omega01.schwarzschild_value.imag:+.6f}i | "
        f"{100.0 * scalar_l2_omega01.fractional_shift_abs:.2f}% |",
        f"| scalar ell=2 | omega0/omega2 | "
        f"{scalar_l2_omega02.value.real:.6f}{scalar_l2_omega02.value.imag:+.6f}i | "
        f"{scalar_l2_omega02.schwarzschild_value.real:.6f}{scalar_l2_omega02.schwarzschild_value.imag:+.6f}i | "
        f"{100.0 * scalar_l2_omega02.fractional_shift_abs:.2f}% |",
        f"| scalar ell=2 | Re/[-Im], n=0 | {scalar_l2_re_damp_0.value.real:.6f} | "
        f"{scalar_l2_re_damp_0.schwarzschild_value.real:.6f} | "
        f"{100.0 * scalar_l2_re_damp_0.fractional_real_shift:.2f}% |",
        "",
        "## Sensitivity Ranking",
        "",
        "| rank | branch | endpoint shift | real shift | damping shift | validation max |",
        "|---:|---|---:|---:|---:|---:|",
    ]
    for row in ranked[:8]:
        lines.append(
            f"| {row.sensitivity_rank} | {row.perturbation_type} ell={row.ell} n={row.overtone} | "
            f"{100.0 * row.endpoint_fractional_complex_shift:.2f}% | "
            f"{100.0 * row.endpoint_fractional_real_shift:.2f}% | "
            f"{100.0 * row.endpoint_fractional_damping_shift:.2f}% | "
            f"{row.max_spectral_leaver_relative_difference:.3e} |"
        )

    lines += [
        "",
        "## Interpretation Guardrails",
        "",
        "- The scalar sector is the cleanest physics target.",
        "- The axial gravitational rows use a KS-lapse-deformed Regge-Wheeler model;",
        "  they are useful phenomenological diagnostics, not a full gauge-invariant KS",
        "  gravitational perturbation derivation.",
        "- Second overtones carry the largest validation and branch-selection uncertainty.",
        "- The endpoint shifts are catalogue diagnostics, not claims about detectability.",
    ]
    output.write_text("\n".join(lines) + "\n")


def read_catalogue_csv(path: Path) -> list[CatalogueRow]:
    rows: list[CatalogueRow] = []
    with path.open(newline="") as handle:
        for item in csv.DictReader(handle):
            literature_real = item["schwarzschild_literature_real"]
            literature_imag = item["schwarzschild_literature_imag"]
            literature = (
                None
                if literature_real == "" or literature_imag == ""
                else complex(float(literature_real), float(literature_imag))
            )
            literature_error = (
                None if item["literature_relative_error"] == "" else float(item["literature_relative_error"])
            )
            literature_tolerance = None if item["literature_tolerance"] == "" else float(item["literature_tolerance"])
            literature_passed = (
                None
                if item["literature_validation_passed"] == ""
                else item["literature_validation_passed"].lower() == "true"
            )
            rows.append(
                CatalogueRow(
                    perturbation_type=item["perturbation_type"],
                    ell=int(item["ell"]),
                    overtone=int(item["overtone"]),
                    mode=item["mode"],
                    a=float(item["a_over_M"]),
                    spectral_n=int(item["spectral_N"]),
                    omega_spectral=complex(float(item["spectral_real"]), float(item["spectral_imag"])),
                    omega_leaver=complex(float(item["leaver_real"]), float(item["leaver_imag"])),
                    spectral_leaver_relative_difference=float(item["spectral_leaver_relative_difference"]),
                    continued_fraction_abs=float(item["continued_fraction_abs"]),
                    schwarzschild_literature=literature,
                    literature_relative_error=literature_error,
                    literature_tolerance=literature_tolerance,
                    spectral_leaver_validation_passed=item["spectral_leaver_validation_passed"].lower() == "true",
                    literature_validation_passed=literature_passed,
                )
            )
    return rows


def write_physics_analysis(
    results_dir: Path,
    figures_dir: Path,
    catalogue_rows: list[CatalogueRow],
) -> dict[str, Path]:
    trend_rows = compute_trend_rows(catalogue_rows)
    summaries = compute_branch_summaries(catalogue_rows)
    ratio_rows = compute_spectroscopic_ratios(catalogue_rows)
    outputs = {
        "trend_csv": results_dir / "catalogue_physics_trends.csv",
        "summary_csv": results_dir / "catalogue_physics_summary.csv",
        "ratio_csv": results_dir / "catalogue_spectroscopic_ratios.csv",
        "report": results_dir / "catalogue_physics_report.md",
        "l2_shifts": figures_dir / "catalogue_l2_fractional_shifts.png",
        "l2_ratios": figures_dir / "catalogue_l2_spectroscopic_ratios.png",
    }
    write_trend_csv(outputs["trend_csv"], trend_rows)
    write_summary_csv(outputs["summary_csv"], summaries)
    write_spectroscopic_ratio_csv(outputs["ratio_csv"], ratio_rows)
    write_physics_report(outputs["report"], summaries, ratio_rows)
    plot_l2_fractional_shifts(trend_rows, outputs["l2_shifts"])
    plot_l2_spectroscopic_ratios(ratio_rows, outputs["l2_ratios"])
    return outputs
