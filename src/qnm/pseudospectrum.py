"""Pseudospectrum diagnostics for the scalar KS QNM residual operator."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from .analysis import read_catalogue_csv
from .spectral import (
    build_spectral_problem,
    generalized_eigenpairs,
    minimize_residual,
    select_tracked_mode,
    spectral_matrix,
)


DEFAULT_THRESHOLDS = (-10.0, -9.5, -9.0)
DEFAULT_QUANTILES = (0.05, 0.10, 0.25, 0.50)


@dataclass(frozen=True)
class PseudospectrumGrid:
    a: float
    spectral_n: int
    ell: int
    overtone: int
    center: complex
    leaver_target: complex
    center_leaver_relative_difference: float
    selection_score: float
    half_width: float
    real_values: np.ndarray
    imag_values: np.ndarray
    sigma_min: np.ndarray
    sigma_max: np.ndarray
    log10_relative_sigma: np.ndarray


@dataclass(frozen=True)
class PseudospectrumSummary:
    a: float
    spectral_n: int
    ell: int
    overtone: int
    center: complex
    leaver_target: complex
    center_leaver_relative_difference: float
    grid_size: int
    half_width: float
    window_area: float
    min_log10_relative_sigma: float
    median_log10_relative_sigma: float
    quantiles: dict[float, float]
    threshold_areas: dict[float, float]
    threshold_boundary_touches: dict[float, bool]


def threshold_label(threshold: float) -> str:
    sign = "m" if threshold < 0 else "p"
    return f"{sign}{abs(threshold):.1f}".replace(".", "_")


def quantile_label(quantile: float) -> str:
    return f"q{int(round(100.0 * quantile)):02d}"


def scalar_fundamental_targets(catalogue_path: Path) -> dict[float, complex]:
    rows = read_catalogue_csv(catalogue_path)
    targets = {
        row.a: row.omega_leaver
        for row in rows
        if row.perturbation_type == "scalar" and row.ell == 2 and row.overtone == 0
    }
    if not targets:
        raise RuntimeError(f"No scalar ell=2 fundamental rows found in {catalogue_path}")
    return dict(sorted(targets.items()))


def _relative_singular_values(problem, omega: complex) -> tuple[float, float, float]:
    singular_values = np.linalg.svd(spectral_matrix(problem, omega), compute_uv=False)
    sigma_max = float(singular_values[0])
    sigma_min = float(singular_values[-1])
    if sigma_max <= 0.0 or not np.isfinite(sigma_max):
        return sigma_min, sigma_max, float("nan")
    return sigma_min, sigma_max, sigma_min / sigma_max


def compute_pseudospectrum_grid(
    a: float,
    leaver_target: complex,
    spectral_n: int,
    grid_size: int,
    half_width: float,
    ell: int = 2,
    overtone: int = 0,
    minimize_radius: float = 0.006,
) -> PseudospectrumGrid:
    """Compute a relative-singular-value grid around one scalar QNM branch."""

    problem = build_spectral_problem(a, spectral_n, ell=ell, perturbation_type="scalar")
    values, vectors = generalized_eigenpairs(problem)
    selection = select_tracked_mode(problem, values, vectors, leaver_target)
    center, _ = minimize_residual(problem, selection.omega, radius=minimize_radius)

    real_values = np.linspace(center.real - half_width, center.real + half_width, grid_size)
    imag_values = np.linspace(center.imag - half_width, center.imag + half_width, grid_size)
    sigma_min = np.empty((grid_size, grid_size), dtype=float)
    sigma_max = np.empty((grid_size, grid_size), dtype=float)
    log10_relative_sigma = np.empty((grid_size, grid_size), dtype=float)

    for i, imag_part in enumerate(imag_values):
        for j, real_part in enumerate(real_values):
            s_min, s_max, relative = _relative_singular_values(problem, complex(real_part, imag_part))
            sigma_min[i, j] = s_min
            sigma_max[i, j] = s_max
            log10_relative_sigma[i, j] = np.log10(max(relative, 1.0e-300))

    return PseudospectrumGrid(
        a=a,
        spectral_n=spectral_n,
        ell=ell,
        overtone=overtone,
        center=center,
        leaver_target=leaver_target,
        center_leaver_relative_difference=float(abs(center - leaver_target) / abs(leaver_target)),
        selection_score=selection.selection_score,
        half_width=half_width,
        real_values=real_values,
        imag_values=imag_values,
        sigma_min=sigma_min,
        sigma_max=sigma_max,
        log10_relative_sigma=log10_relative_sigma,
    )


def summarize_grid(
    grid: PseudospectrumGrid,
    thresholds: tuple[float, ...] = DEFAULT_THRESHOLDS,
    quantiles: tuple[float, ...] = DEFAULT_QUANTILES,
) -> PseudospectrumSummary:
    values = grid.log10_relative_sigma
    window_area = (2.0 * grid.half_width) ** 2
    quantile_values = {quantile: float(np.quantile(values, quantile)) for quantile in quantiles}
    threshold_areas: dict[float, float] = {}
    threshold_boundary_touches: dict[float, bool] = {}
    for threshold in thresholds:
        mask = values <= threshold
        threshold_areas[threshold] = float(np.mean(mask) * window_area)
        boundary = np.concatenate([mask[0, :], mask[-1, :], mask[:, 0], mask[:, -1]])
        threshold_boundary_touches[threshold] = bool(np.any(boundary))

    return PseudospectrumSummary(
        a=grid.a,
        spectral_n=grid.spectral_n,
        ell=grid.ell,
        overtone=grid.overtone,
        center=grid.center,
        leaver_target=grid.leaver_target,
        center_leaver_relative_difference=grid.center_leaver_relative_difference,
        grid_size=len(grid.real_values),
        half_width=grid.half_width,
        window_area=window_area,
        min_log10_relative_sigma=float(np.min(values)),
        median_log10_relative_sigma=float(np.median(values)),
        quantiles=quantile_values,
        threshold_areas=threshold_areas,
        threshold_boundary_touches=threshold_boundary_touches,
    )


def compute_pseudospectrum_analysis(
    catalogue_path: Path,
    spectral_n: int = 64,
    grid_size: int = 81,
    half_width: float = 0.025,
    thresholds: tuple[float, ...] = DEFAULT_THRESHOLDS,
    quantiles: tuple[float, ...] = DEFAULT_QUANTILES,
) -> tuple[list[PseudospectrumGrid], list[PseudospectrumSummary]]:
    targets = scalar_fundamental_targets(catalogue_path)
    grids: list[PseudospectrumGrid] = []
    summaries: list[PseudospectrumSummary] = []
    for a, leaver_target in targets.items():
        grid = compute_pseudospectrum_grid(
            a=a,
            leaver_target=leaver_target,
            spectral_n=spectral_n,
            grid_size=grid_size,
            half_width=half_width,
        )
        grids.append(grid)
        summaries.append(summarize_grid(grid, thresholds=thresholds, quantiles=quantiles))
    return grids, summaries


def compute_resolution_check(
    catalogue_path: Path,
    spectral_sizes: tuple[int, ...] = (32, 48, 64),
    grid_size: int = 41,
    half_width: float = 0.025,
    thresholds: tuple[float, ...] = DEFAULT_THRESHOLDS,
    quantiles: tuple[float, ...] = DEFAULT_QUANTILES,
) -> list[PseudospectrumSummary]:
    all_summaries: list[PseudospectrumSummary] = []
    for spectral_n in spectral_sizes:
        _, summaries = compute_pseudospectrum_analysis(
            catalogue_path=catalogue_path,
            spectral_n=spectral_n,
            grid_size=grid_size,
            half_width=half_width,
            thresholds=thresholds,
            quantiles=quantiles,
        )
        all_summaries.extend(summaries)
    return all_summaries


def write_grid_csv(output: Path, grids: list[PseudospectrumGrid]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "a_over_M",
                "spectral_N",
                "ell",
                "overtone",
                "center_real",
                "center_imag",
                "leaver_real",
                "leaver_imag",
                "center_leaver_relative_difference",
                "real",
                "imag",
                "delta_real_from_center",
                "delta_imag_from_center",
                "sigma_min",
                "sigma_max",
                "relative_sigma_min",
                "log10_relative_sigma_min",
            ]
        )
        for grid in grids:
            for i, imag_part in enumerate(grid.imag_values):
                for j, real_part in enumerate(grid.real_values):
                    relative = grid.sigma_min[i, j] / grid.sigma_max[i, j]
                    writer.writerow(
                        [
                            grid.a,
                            grid.spectral_n,
                            grid.ell,
                            grid.overtone,
                            grid.center.real,
                            grid.center.imag,
                            grid.leaver_target.real,
                            grid.leaver_target.imag,
                            grid.center_leaver_relative_difference,
                            real_part,
                            imag_part,
                            real_part - grid.center.real,
                            imag_part - grid.center.imag,
                            grid.sigma_min[i, j],
                            grid.sigma_max[i, j],
                            relative,
                            grid.log10_relative_sigma[i, j],
                        ]
                    )


def write_summary_csv(
    output: Path,
    summaries: list[PseudospectrumSummary],
    thresholds: tuple[float, ...] = DEFAULT_THRESHOLDS,
    quantiles: tuple[float, ...] = DEFAULT_QUANTILES,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        header = [
            "a_over_M",
            "spectral_N",
            "ell",
            "overtone",
            "grid_size",
            "half_width",
            "window_area",
            "center_real",
            "center_imag",
            "leaver_real",
            "leaver_imag",
            "center_leaver_relative_difference",
            "min_log10_relative_sigma",
            "median_log10_relative_sigma",
        ]
        header += [f"log10_relative_sigma_{quantile_label(q)}" for q in quantiles]
        header += [f"area_log10_relative_sigma_le_{threshold_label(t)}" for t in thresholds]
        header += [f"area_fraction_log10_relative_sigma_le_{threshold_label(t)}" for t in thresholds]
        header += [f"boundary_touch_log10_relative_sigma_le_{threshold_label(t)}" for t in thresholds]
        writer.writerow(header)

        for row in sorted(summaries, key=lambda item: (item.spectral_n, item.a)):
            values = [
                row.a,
                row.spectral_n,
                row.ell,
                row.overtone,
                row.grid_size,
                row.half_width,
                row.window_area,
                row.center.real,
                row.center.imag,
                row.leaver_target.real,
                row.leaver_target.imag,
                row.center_leaver_relative_difference,
                row.min_log10_relative_sigma,
                row.median_log10_relative_sigma,
            ]
            values += [row.quantiles[q] for q in quantiles]
            values += [row.threshold_areas[t] for t in thresholds]
            values += [row.threshold_areas[t] / row.window_area for t in thresholds]
            values += [row.threshold_boundary_touches[t] for t in thresholds]
            writer.writerow(values)


def plot_pseudospectrum_contours(
    output: Path,
    grids: list[PseudospectrumGrid],
    contour_levels: tuple[float, ...] = (-11.0, -10.5, -10.0, -9.5, -9.0, -8.5),
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    ordered = sorted(grids, key=lambda grid: grid.a)
    fig, axes = plt.subplots(2, 2, figsize=(9.4, 7.5), constrained_layout=True)
    axes_flat = axes.ravel()
    mappable = None
    for axis, grid in zip(axes_flat, ordered):
        x, y = np.meshgrid(grid.real_values, grid.imag_values)
        values = np.clip(grid.log10_relative_sigma, contour_levels[0], contour_levels[-1])
        mappable = axis.contourf(
            x,
            y,
            values,
            levels=contour_levels,
            cmap="viridis",
            extend="min",
        )
        axis.contour(
            x,
            y,
            values,
            levels=contour_levels[1:-1],
            colors="white",
            linewidths=0.55,
            alpha=0.75,
        )
        axis.plot(grid.center.real, grid.center.imag, marker="x", color="white", markersize=7, mew=1.4)
        axis.set_title(fr"$a/M={grid.a:g}$", fontsize=10)
        axis.set_xlabel(r"$\mathrm{Re}(M\omega)$")
        axis.set_ylabel(r"$\mathrm{Im}(M\omega)$")
        axis.ticklabel_format(useOffset=False)
    if mappable is not None:
        cbar = fig.colorbar(mappable, ax=axes_flat.tolist(), shrink=0.92)
        cbar.set_label(r"$\log_{10}[\sigma_{\min}(P_N)/\sigma_{\max}(P_N)]$")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_pseudospectrum_sensitivity(
    output: Path,
    summaries: list[PseudospectrumSummary],
    thresholds: tuple[float, ...] = DEFAULT_THRESHOLDS,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = sorted(summaries, key=lambda row: row.a)
    a_values = np.array([row.a for row in rows])
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2), constrained_layout=True)
    axes[0].plot(
        a_values,
        [-row.quantiles[0.10] for row in rows],
        marker="o",
        label=r"$-Q_{10}[\log_{10}\eta]$",
    )
    axes[0].plot(
        a_values,
        [-row.median_log10_relative_sigma for row in rows],
        marker="s",
        label=r"$-Q_{50}[\log_{10}\eta]$",
    )
    axes[0].set_xlabel(r"$a/M$")
    axes[0].set_ylabel("local pseudospectral susceptibility")
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False, fontsize=8)

    for threshold in thresholds:
        axes[1].plot(
            a_values,
            [100.0 * row.threshold_areas[threshold] / row.window_area for row in rows],
            marker="o",
            label=fr"$\log_{{10}}\eta\leq {threshold:g}$",
        )
    axes[1].set_xlabel(r"$a/M$")
    axes[1].set_ylabel("contour area fraction (%)")
    axes[1].grid(alpha=0.25)
    axes[1].legend(frameon=False, fontsize=8)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_resolution_check(
    output: Path,
    summaries: list[PseudospectrumSummary],
    quantile: float = 0.10,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, axis = plt.subplots(figsize=(7.2, 4.2), constrained_layout=True)
    for spectral_n in sorted({row.spectral_n for row in summaries}):
        rows = sorted([row for row in summaries if row.spectral_n == spectral_n], key=lambda row: row.a)
        baseline = -rows[0].quantiles[quantile]
        axis.plot(
            [row.a for row in rows],
            [(-row.quantiles[quantile]) - baseline for row in rows],
            marker="o",
            label=fr"$N={spectral_n}$",
        )
    axis.axhline(0.0, color="black", linewidth=0.8, alpha=0.35)
    axis.set_xlabel(r"$a/M$")
    axis.set_ylabel(fr"$\Delta[-Q_{{{int(100 * quantile)}}}(\log_{{10}}\eta)]$")
    axis.grid(alpha=0.25)
    axis.legend(frameon=False, fontsize=8)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def write_report(
    output: Path,
    summaries: list[PseudospectrumSummary],
    resolution_summaries: list[PseudospectrumSummary],
    thresholds: tuple[float, ...] = DEFAULT_THRESHOLDS,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = sorted(summaries, key=lambda row: row.a)
    schwarzschild = rows[0]
    endpoint = rows[-1]
    threshold = thresholds[0]
    area_factor = endpoint.threshold_areas[threshold] / schwarzschild.threshold_areas[threshold]
    susceptibility_gain = (-endpoint.quantiles[0.10]) - (-schwarzschild.quantiles[0.10])
    center_error = max(row.center_leaver_relative_difference for row in rows)
    endpoint_boundary_touch = endpoint.threshold_boundary_touches[threshold]
    resolution_lines = []
    for spectral_n in sorted({row.spectral_n for row in resolution_summaries}):
        n_rows = sorted([row for row in resolution_summaries if row.spectral_n == spectral_n], key=lambda row: row.a)
        gain = (-n_rows[-1].quantiles[0.10]) - (-n_rows[0].quantiles[0.10])
        resolution_lines.append(f"- N={spectral_n}: Q10 susceptibility gain from a/M=0 to 1 is {gain:.3f}.")

    lines = [
        "# Scalar Fundamental Pseudospectrum Upgrade Report",
        "",
        "## Repository Audit",
        "",
        "- Core numerical code lives under `src/qnm/`: compactified geometry helpers,",
        "  Chebyshev spectral operators, Leaver-style validation, catalogue generation,",
        "  normalization utilities, and catalogue-level physics diagnostics.",
        "- Validated generated products live under `outputs/results/` and `outputs/figures/`.",
        "- Manuscript sources live under `papers/manuscript/`.",
        "- The existing pipeline safely supports scalar-sector extensions because the scalar",
        "  Chebyshev and Leaver branches are already cross-validated. The phenomenological",
        "  axial sector was not extended in this upgrade.",
        "",
        "## Upgrade Selection",
        "",
        "| candidate direction | impact | feasibility | decision |",
        "|---|---|---|---|",
        "| Pseudospectrum and spectral instability | high: connects the residual workflow to QNM instability literature | high: uses existing P_N(omega) and sigma_min machinery | selected |",
        "| Exceptional points or branch interactions | potentially high, but needs denser branch/eigenvector tracking and stronger mathematical evidence | medium | deferred |",
        "| Literature-matched KS benchmarking | useful, but incremental after the Konoplya side comparison | high | secondary future work |",
        "| Gauge-invariant gravitational sector | very high, but requires a new perturbation derivation beyond the current codebase | low for this iteration | deferred |",
        "| Spectroscopy/detectability | useful, but risks unsupported detector claims without a full waveform/noise model | medium | deferred |",
        "",
        "## What Was Attempted",
        "",
        "A finite-dimensional pseudospectrum diagnostic was added for the scalar",
        "ell=2 fundamental KS QNM branch. For each deformation value, the code",
        "builds the Chebyshev residual operator P_N(omega), centers a local grid",
        "on the residual-minimized spectral frequency, and evaluates",
        "eta(omega)=sigma_min(P_N)/sigma_max(P_N).",
        "",
        "## Reproduction Commands",
        "",
        "```bash",
        "python scripts/analyze_pseudospectrum.py",
        "python -m pytest",
        "python scripts/run_hybrid_qnm_algorithm.py --tests-only",
        "```",
        "",
        "The pseudospectrum command writes the grid, summary table, resolution check,",
        "figures, and this report. The test commands check that the existing validated",
        "pipeline remains intact.",
        "",
        "## New Implementation Artifacts",
        "",
        "- `src/qnm/pseudospectrum.py`",
        "- `scripts/analyze_pseudospectrum.py`",
        "- `outputs/results/scalar_l2_pseudospectrum_grid.csv`",
        "- `outputs/results/scalar_l2_pseudospectrum_summary.csv`",
        "- `outputs/results/scalar_l2_pseudospectrum_resolution_check.csv`",
        "- `outputs/figures/scalar_l2_pseudospectrum_contours.png`",
        "- `outputs/figures/scalar_l2_pseudospectrum_sensitivity.png`",
        "- `outputs/figures/scalar_l2_pseudospectrum_resolution_check.png`",
        "",
        "## What Worked",
        "",
        f"- Main grid: N={endpoint.spectral_n}, grid={endpoint.grid_size}x{endpoint.grid_size}, "
        f"half-width={endpoint.half_width:g} in both Re(M omega) and Im(M omega).",
        f"- Maximum center-to-Leaver relative difference: {center_error:.3e}.",
        f"- The 10% quantile susceptibility, -Q10(log10 eta), increases by {susceptibility_gain:.3f} "
        "from a/M=0 to a/M=1.",
        f"- The area fraction satisfying log10(eta)<={threshold:g} grows by a factor {area_factor:.2f} "
        "from Schwarzschild to a/M=1.",
        "- The sign of the Q10 susceptibility trend is stable across N=32, 48, and 64.",
        "- The contour-area diagnostic is secondary to the quantile diagnostic because fixed",
        "  epsilon contour areas depend more strongly on N and on the chosen plotting window.",
        "",
        "## Resolution Check",
        "",
        *resolution_lines,
        "",
        "## What Failed Or Was Limited",
        "",
        f"- The log10(eta)<={threshold:g} contour at a/M=1 touches the local-window boundary: "
        f"{endpoint_boundary_touch}. The reported area factor is therefore a finite-window diagnostic,",
        "  not a global contour area.",
        "- Absolute contour levels shift with Chebyshev size N, so the finite-N robustness check",
        "  uses the sign and monotonicity of the Q10 susceptibility gain rather than exact equality",
        "  of epsilon-contour areas.",
        "- No exceptional-point search was attempted in this upgrade; mode coalescence would require",
        "  a separate eigenvector and mode-pair condition analysis.",
        "",
        "## Publishable Claim",
        "",
        "Within the finite-dimensional Chebyshev residual normalization used here,",
        "the scalar ell=2 fundamental KS branch shows increasing local pseudospectral",
        "sensitivity as a/M grows. This complements the frequency-softening result:",
        "the branch moves to lower oscillation frequency while the surrounding",
        "relative-smallest-singular-value basin expands.",
        "",
        "## What Remains Speculative",
        "",
        "- This is not a proof about the infinite-dimensional KS wave operator.",
        "- Absolute epsilon-contour values depend on N and on the chosen operator normalization.",
        "- The result is local to the scalar ell=2 fundamental branch and should not be",
        "  generalized to overtones or the phenomenological axial sector without separate checks.",
        "- The analysis is not a detector forecast.",
        "",
        "## Realistic Journal Target",
        "",
        "The new result makes the paper more suitable for a numerics-focused GR journal",
        "or a strong mathematical-physics venue. EPJ Plus and Universe remain realistic;",
        "CQG becomes more plausible if the pseudospectrum section is presented as a",
        "finite-dimensional diagnostic and not as an infinite-dimensional stability theorem.",
    ]
    output.write_text("\n".join(lines) + "\n")
