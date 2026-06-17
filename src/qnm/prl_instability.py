"""PRL-level stress test for KS scalar QNM pseudospectral instability claims.

This module deliberately separates a broad discovery-claim scan from the
publication catalogue.  The goal is not to certify new modes, but to test
whether the scalar KS branches support a simple Letter-level statement about
ringdown softening and increased non-Hermitian sensitivity.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eig
from matplotlib.lines import Line2D

from .common import SCHWARZSCHILD_REFERENCES, f_ks, horizon_radius, scalar_potential
from .leaver import solve_leaver_mode
from .spectral import (
    build_spectral_problem,
    generalized_eigenpairs,
    minimize_residual,
    mode_overlap,
    scale_generalized_pencil,
    select_tracked_mode,
    spectral_matrix,
)


PRL_A_VALUES = (0.0, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0)
PRL_ELL_VALUES = (0, 1, 2, 3, 4)
PRL_OVERTONES = (0, 1, 2)
PRL_SPECTRAL_SIZES = (32, 48, 64)
PRL_N96_BRANCHES = tuple((ell, 0) for ell in PRL_ELL_VALUES)
PRL_THRESHOLDS = (-10.0, -9.5, -9.0)
PRL_QUANTILES = (0.10, 0.50)

# The ell=0 and ell=1 estimates are used only to initialize branch tracking.
# Endpoint Leaver checks are written separately; unreliable branches are not
# allowed into a PRL-level claim.
SCALAR_INITIAL_TARGETS: dict[tuple[int, int], complex | None] = {
    (0, 0): 0.11045493636871656 - 0.10489571229328226j,
    (0, 1): 0.08615890818056095 - 0.3480584541826301j,
    # The next ell=0 branch is a purely imaginary Leaver root in this setup and
    # is not tracked by the oscillatory-mode selection filter used below.
    (0, 2): None,
    (1, 0): 0.2929361332672828 - 0.09765998891357816j,
    (1, 1): 0.2644486505967102 - 0.30625739154784654j,
    (1, 2): 0.22953931525220106 - 0.5401334449274291j,
}
for ell in (2, 3, 4):
    for overtone in PRL_OVERTONES:
        SCALAR_INITIAL_TARGETS[(ell, overtone)] = SCHWARZSCHILD_REFERENCES[("scalar", ell, overtone)]


@dataclass
class BranchScanRow:
    ell: int
    overtone: int
    a: float
    spectral_n: int
    center: complex
    initial_eigenvalue: complex
    residual_norm: float
    selection_score: float
    eigenvector_overlap_previous_a: float | None
    half_width: float
    grid_size: int
    min_log10_relative_sigma: float
    q10_log10_relative_sigma: float
    q50_log10_relative_sigma: float
    q10_log10_sigma_min: float
    q50_log10_sigma_min: float
    area_fraction_m10: float
    area_fraction_m9_5: float
    area_fraction_m9: float
    boundary_touch_m10: bool
    curvature_trace_indicator: float
    eigen_condition_indicator: float
    nearest_eigenvalue_distance: float
    hit_minimizer_boundary: bool
    shape: np.ndarray


@dataclass
class BranchVerdict:
    ell: int
    overtone: int
    q10_gain_n64: float
    q50_gain_n64: float
    absolute_q10_gain_n64: float
    real_frequency_shift_n64: float
    damping_shift_n64: float
    monotonic_q10_n64: bool
    positive_gain_all_tested_n: bool
    max_center_spread_endpoint: float
    reliability: str
    prl_support: bool
    comment: str


@dataclass
class BarrierMetric:
    ell: int
    a: float
    peak_height: float
    r_peak: float
    rstar_width_halfmax: float
    curvature_rstar_at_peak: float


@dataclass
class ModePairDiagnostic:
    ell: int
    a: float
    spectral_n: int
    overtone_i: int
    overtone_j: int
    frequency_distance: float
    normalized_frequency_distance: float
    mode_shape_overlap: float
    possible_coalescence: bool


def active_branches() -> list[tuple[int, int]]:
    return [
        (ell, overtone)
        for ell in PRL_ELL_VALUES
        for overtone in PRL_OVERTONES
        if SCALAR_INITIAL_TARGETS.get((ell, overtone)) is not None
    ]


def _relative_and_absolute_singular_values(problem, omega: complex) -> tuple[float, float, float]:
    singular_values = np.linalg.svd(spectral_matrix(problem, omega), compute_uv=False)
    sigma_max = float(singular_values[0])
    sigma_min = float(singular_values[-1])
    relative = sigma_min / sigma_max if sigma_max > 0.0 else float("nan")
    return sigma_min, sigma_max, relative


def _condition_indicator(problem, target: complex) -> tuple[float, float]:
    left_scaled, right_scaled, _ = scale_generalized_pencil(problem.left, problem.right)
    values, left_vectors, right_vectors = eig(left_scaled, right_scaled, left=True, right=True)
    finite_indices = [
        index
        for index, value in enumerate(values)
        if np.isfinite(value)
        and value.real > 0.03
        and value.imag < -0.02
        and value.real < 2.0
        and value.imag > -3.0
    ]
    if not finite_indices:
        return float("nan"), float("nan")
    index = min(finite_indices, key=lambda item: abs(values[item] - target))
    right = right_vectors[:, index]
    left = left_vectors[:, index]
    denominator = abs(np.vdot(left, right_scaled @ right))
    if denominator == 0.0 or not np.isfinite(denominator):
        condition = float("inf")
    else:
        condition = float(np.linalg.norm(left) * np.linalg.norm(right) / denominator)
    distances = [abs(values[index] - values[other]) for other in finite_indices if other != index]
    nearest = float(min(distances)) if distances else float("nan")
    return condition, nearest


def _half_width(center: complex) -> float:
    return max(0.012, 0.05 * abs(center))


def _grid_metrics(problem, center: complex, grid_size: int, half_width: float) -> dict[str, float | bool]:
    real_values = np.linspace(center.real - half_width, center.real + half_width, grid_size)
    imag_values = np.linspace(center.imag - half_width, center.imag + half_width, grid_size)
    relative_log = np.empty((grid_size, grid_size), dtype=float)
    absolute_log = np.empty((grid_size, grid_size), dtype=float)
    for i, imag_part in enumerate(imag_values):
        for j, real_part in enumerate(real_values):
            sigma_min, _, relative = _relative_and_absolute_singular_values(
                problem,
                complex(real_part, imag_part),
            )
            relative_log[i, j] = np.log10(max(relative, 1.0e-300))
            absolute_log[i, j] = np.log10(max(sigma_min, 1.0e-300))

    metrics: dict[str, float | bool] = {
        "min_log10_relative_sigma": float(np.min(relative_log)),
        "q10_log10_relative_sigma": float(np.quantile(relative_log, 0.10)),
        "q50_log10_relative_sigma": float(np.quantile(relative_log, 0.50)),
        "q10_log10_sigma_min": float(np.quantile(absolute_log, 0.10)),
        "q50_log10_sigma_min": float(np.quantile(absolute_log, 0.50)),
    }
    for threshold in PRL_THRESHOLDS:
        mask = relative_log <= threshold
        suffix = str(abs(threshold)).replace(".", "_")
        metrics[f"area_fraction_m{suffix}"] = float(np.mean(mask))
        boundary = np.concatenate([mask[0, :], mask[-1, :], mask[:, 0], mask[:, -1]])
        metrics[f"boundary_touch_m{suffix}"] = bool(np.any(boundary))

    center_index = grid_size // 2
    susceptibility = -relative_log
    dr = real_values[1] - real_values[0]
    di = imag_values[1] - imag_values[0]
    if 0 < center_index < grid_size - 1:
        d2_re = (
            susceptibility[center_index, center_index + 1]
            - 2.0 * susceptibility[center_index, center_index]
            + susceptibility[center_index, center_index - 1]
        ) / (dr * dr)
        d2_im = (
            susceptibility[center_index + 1, center_index]
            - 2.0 * susceptibility[center_index, center_index]
            + susceptibility[center_index - 1, center_index]
        ) / (di * di)
        metrics["curvature_trace_indicator"] = float(-(d2_re + d2_im))
    else:
        metrics["curvature_trace_indicator"] = float("nan")
    return metrics


def scan_branch(
    ell: int,
    overtone: int,
    spectral_n: int,
    a_values: tuple[float, ...],
    grid_size: int,
    initial_target: complex,
) -> list[BranchScanRow]:
    rows: list[BranchScanRow] = []
    target = initial_target
    previous_omega: complex | None = None
    previous_shape: np.ndarray | None = None
    previous_s_nodes: np.ndarray | None = None
    for a in a_values:
        problem = build_spectral_problem(a, spectral_n, ell=ell, perturbation_type="scalar")
        values, vectors = generalized_eigenpairs(problem)
        selection = select_tracked_mode(
            problem,
            values,
            vectors,
            target,
            previous_omega=previous_omega,
            previous_shape=previous_shape,
            previous_s_nodes=previous_s_nodes,
        )
        radius = 0.010 if overtone == 0 else 0.016
        center, residual = minimize_residual(problem, selection.omega, radius=radius)
        half_width = _half_width(center)
        metrics = _grid_metrics(problem, center, grid_size=grid_size, half_width=half_width)
        condition, nearest = _condition_indicator(problem, center)
        hit_boundary = (
            abs(center.real - selection.omega.real) > 0.96 * radius
            or abs(center.imag - selection.omega.imag) > 0.96 * radius
        )
        rows.append(
            BranchScanRow(
                ell=ell,
                overtone=overtone,
                a=a,
                spectral_n=spectral_n,
                center=center,
                initial_eigenvalue=selection.omega,
                residual_norm=residual,
                selection_score=selection.selection_score,
                eigenvector_overlap_previous_a=selection.eigenvector_overlap,
                half_width=half_width,
                grid_size=grid_size,
                min_log10_relative_sigma=float(metrics["min_log10_relative_sigma"]),
                q10_log10_relative_sigma=float(metrics["q10_log10_relative_sigma"]),
                q50_log10_relative_sigma=float(metrics["q50_log10_relative_sigma"]),
                q10_log10_sigma_min=float(metrics["q10_log10_sigma_min"]),
                q50_log10_sigma_min=float(metrics["q50_log10_sigma_min"]),
                area_fraction_m10=float(metrics["area_fraction_m10_0"]),
                area_fraction_m9_5=float(metrics["area_fraction_m9_5"]),
                area_fraction_m9=float(metrics["area_fraction_m9_0"]),
                boundary_touch_m10=bool(metrics["boundary_touch_m10_0"]),
                curvature_trace_indicator=float(metrics["curvature_trace_indicator"]),
                eigen_condition_indicator=condition,
                nearest_eigenvalue_distance=nearest,
                hit_minimizer_boundary=hit_boundary,
                shape=selection.shape,
            )
        )
        target = center
        previous_omega = center
        previous_shape = selection.shape
        previous_s_nodes = problem.s_nodes
    return rows


def run_scan(
    a_values: tuple[float, ...] = PRL_A_VALUES,
    spectral_sizes: tuple[int, ...] = PRL_SPECTRAL_SIZES,
    main_grid_size: int = 31,
    resolution_grid_size: int = 21,
    include_n96_fundamentals: bool = True,
) -> list[BranchScanRow]:
    rows: list[BranchScanRow] = []
    branches = active_branches()
    for spectral_n in spectral_sizes:
        grid_size = main_grid_size if spectral_n == 64 else resolution_grid_size
        for ell, overtone in branches:
            target = SCALAR_INITIAL_TARGETS[(ell, overtone)]
            if target is None:
                continue
            rows.extend(scan_branch(ell, overtone, spectral_n, a_values, grid_size, target))
    if include_n96_fundamentals:
        for ell, overtone in PRL_N96_BRANCHES:
            target = SCALAR_INITIAL_TARGETS[(ell, overtone)]
            if target is None:
                continue
            rows.extend(scan_branch(ell, overtone, 96, a_values, resolution_grid_size, target))
    return rows


def summarize_branches(rows: list[BranchScanRow]) -> list[BranchVerdict]:
    verdicts: list[BranchVerdict] = []
    for ell, overtone in active_branches():
        branch_rows = [row for row in rows if row.ell == ell and row.overtone == overtone]
        if not branch_rows:
            continue
        n64 = sorted([row for row in branch_rows if row.spectral_n == 64], key=lambda row: row.a)
        if len(n64) < 2:
            continue
        susceptibilities = [-row.q10_log10_relative_sigma for row in n64]
        medians = [-row.q50_log10_relative_sigma for row in n64]
        q10_gain = susceptibilities[-1] - susceptibilities[0]
        q50_gain = medians[-1] - medians[0]
        absolute_gain = (-n64[-1].q10_log10_sigma_min) - (-n64[0].q10_log10_sigma_min)
        real_shift = n64[-1].center.real / n64[0].center.real - 1.0
        damping_shift = (-n64[-1].center.imag) / (-n64[0].center.imag) - 1.0
        monotonic = all(b >= a - 1.0e-3 for a, b in zip(susceptibilities, susceptibilities[1:]))
        gains_by_n: list[float] = []
        for spectral_n in sorted({row.spectral_n for row in branch_rows}):
            n_rows = sorted([row for row in branch_rows if row.spectral_n == spectral_n], key=lambda row: row.a)
            if len(n_rows) >= 2:
                gains_by_n.append((-n_rows[-1].q10_log10_relative_sigma) - (-n_rows[0].q10_log10_relative_sigma))
        positive_all = bool(gains_by_n) and all(gain > 0.0 for gain in gains_by_n)
        endpoint_rows = [row for row in branch_rows if abs(row.a - 1.0) < 1.0e-12 and row.spectral_n in (32, 48, 64)]
        centers = [row.center for row in endpoint_rows]
        max_spread = 0.0
        for i, left in enumerate(centers):
            for right in centers[i + 1 :]:
                max_spread = max(max_spread, abs(left - right) / max(abs(left), 0.1))
        any_boundary = any(row.hit_minimizer_boundary for row in branch_rows)
        any_contour_boundary = any(row.boundary_touch_m10 for row in n64)
        if max_spread > 5.0e-2 or any_boundary:
            reliability = "exploratory"
        elif overtone == 2 or ell == 0:
            reliability = "caution"
        else:
            reliability = "usable"
        prl_support = reliability == "usable" and monotonic and positive_all and q10_gain > 0.0
        comments: list[str] = []
        if any_boundary:
            comments.append("minimizer hit search boundary")
        if any_contour_boundary:
            comments.append("deep contour touches local window")
        if not monotonic:
            comments.append("N=64 susceptibility not monotonic in a/M")
        if not positive_all:
            comments.append("endpoint gain not positive for every tested N")
        if max_spread > 5.0e-2:
            comments.append("endpoint center spread across N is large")
        if not comments:
            comments.append("trend passes finite-N scalar checks")
        verdicts.append(
            BranchVerdict(
                ell=ell,
                overtone=overtone,
                q10_gain_n64=q10_gain,
                q50_gain_n64=q50_gain,
                absolute_q10_gain_n64=absolute_gain,
                real_frequency_shift_n64=real_shift,
                damping_shift_n64=damping_shift,
                monotonic_q10_n64=monotonic,
                positive_gain_all_tested_n=positive_all,
                max_center_spread_endpoint=max_spread,
                reliability=reliability,
                prl_support=prl_support,
                comment="; ".join(comments),
            )
        )
    return verdicts


def compute_barrier_metrics(
    ell_values: tuple[int, ...] = PRL_ELL_VALUES,
    a_values: tuple[float, ...] = PRL_A_VALUES,
) -> list[BarrierMetric]:
    metrics: list[BarrierMetric] = []
    for ell in ell_values:
        for a in a_values:
            rh = horizon_radius(a)
            r = np.geomspace(rh * (1.0 + 1.0e-5), 120.0, 8000)
            potential = scalar_potential(r, ell, a)
            peak_index = int(np.argmax(potential))
            rstar = np.zeros_like(r)
            integrand = 1.0 / np.asarray(f_ks(r, a), dtype=float)
            dr = np.diff(r)
            rstar[1:] = np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * dr)
            peak = float(potential[peak_index])
            half = 0.5 * peak
            left = peak_index
            while left > 0 and potential[left] > half:
                left -= 1
            right = peak_index
            while right < len(r) - 1 and potential[right] > half:
                right += 1
            width = float(rstar[right] - rstar[left]) if right > left else float("nan")
            if 2 <= peak_index < len(r) - 2:
                local_x = rstar[peak_index - 2 : peak_index + 3]
                local_y = potential[peak_index - 2 : peak_index + 3]
                coeffs = np.polyfit(local_x - local_x[2], local_y, 2)
                curvature = float(2.0 * coeffs[0])
            else:
                curvature = float("nan")
            metrics.append(
                BarrierMetric(
                    ell=ell,
                    a=a,
                    peak_height=peak,
                    r_peak=float(r[peak_index]),
                    rstar_width_halfmax=width,
                    curvature_rstar_at_peak=curvature,
                )
            )
    return metrics


def compute_pair_diagnostics(rows: list[BranchScanRow], spectral_n: int = 64) -> list[ModePairDiagnostic]:
    diagnostics: list[ModePairDiagnostic] = []
    for ell in PRL_ELL_VALUES:
        for a in PRL_A_VALUES:
            local = sorted(
                [row for row in rows if row.ell == ell and row.a == a and row.spectral_n == spectral_n],
                key=lambda row: row.overtone,
            )
            for i, left in enumerate(local):
                for right in local[i + 1 :]:
                    distance = abs(left.center - right.center)
                    normalized = distance / max(abs(left.center), abs(right.center), 0.1)
                    overlap = mode_overlap(left.shape, right.shape)
                    diagnostics.append(
                        ModePairDiagnostic(
                            ell=ell,
                            a=a,
                            spectral_n=spectral_n,
                            overtone_i=left.overtone,
                            overtone_j=right.overtone,
                            frequency_distance=float(distance),
                            normalized_frequency_distance=float(normalized),
                            mode_shape_overlap=overlap,
                            possible_coalescence=bool(normalized < 0.05 and overlap > 0.95),
                        )
                    )
    return diagnostics


def compute_window_sensitivity(
    rows: list[BranchScanRow],
    scales: tuple[float, ...] = (0.75, 1.25),
    grid_size: int = 21,
) -> list[dict[str, float | int]]:
    output: list[dict[str, float | int]] = []
    endpoints = [
        row
        for row in rows
        if row.spectral_n == 64 and (abs(row.a - 0.0) < 1.0e-12 or abs(row.a - 1.0) < 1.0e-12)
    ]
    for row in endpoints:
        for scale in scales:
            problem = build_spectral_problem(row.a, row.spectral_n, ell=row.ell, perturbation_type="scalar")
            metrics = _grid_metrics(problem, row.center, grid_size=grid_size, half_width=row.half_width * scale)
            output.append(
                {
                    "ell": row.ell,
                    "overtone": row.overtone,
                    "a_over_M": row.a,
                    "spectral_N": row.spectral_n,
                    "window_scale": scale,
                    "q10_susceptibility": -float(metrics["q10_log10_relative_sigma"]),
                    "q50_susceptibility": -float(metrics["q50_log10_relative_sigma"]),
                    "area_fraction_m10": float(metrics["area_fraction_m10_0"]),
                }
            )
    return output


def run_endpoint_leaver_checks(rows: list[BranchScanRow]) -> list[dict[str, float | int | str]]:
    checks: list[dict[str, float | int | str]] = []
    endpoints = [
        row
        for row in rows
        if row.spectral_n == 64 and (abs(row.a - 0.0) < 1.0e-12 or abs(row.a - 1.0) < 1.0e-12)
    ]
    for row in endpoints:
        try:
            omega_leaver, cf_abs = solve_leaver_mode(
                row.a,
                row.center,
                ell=row.ell,
                perturbation_type="scalar",
                depth=180,
                order=80,
                residual_tolerance=1.0e-3,
            )
            relative_difference = abs(omega_leaver - row.center) / max(abs(row.center), 0.1)
            status = "ok"
        except Exception as exc:  # pragma: no cover - diagnostic path
            omega_leaver = complex(float("nan"), float("nan"))
            cf_abs = float("nan")
            relative_difference = float("nan")
            status = f"failed: {type(exc).__name__}"
        checks.append(
            {
                "ell": row.ell,
                "overtone": row.overtone,
                "a_over_M": row.a,
                "spectral_N": row.spectral_n,
                "center_real": row.center.real,
                "center_imag": row.center.imag,
                "leaver_real": omega_leaver.real,
                "leaver_imag": omega_leaver.imag,
                "relative_difference": float(relative_difference),
                "continued_fraction_abs": float(cf_abs),
                "status": status,
            }
        )
    return checks


def _write_dict_csv(output: Path, rows: list[dict]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        output.write_text("")
        return
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_scan_csv(output: Path, rows: list[BranchScanRow]) -> None:
    dict_rows = []
    for row in rows:
        dict_rows.append(
            {
                "ell": row.ell,
                "overtone": row.overtone,
                "a_over_M": row.a,
                "spectral_N": row.spectral_n,
                "center_real": row.center.real,
                "center_imag": row.center.imag,
                "initial_eigenvalue_real": row.initial_eigenvalue.real,
                "initial_eigenvalue_imag": row.initial_eigenvalue.imag,
                "residual_norm": row.residual_norm,
                "selection_score": row.selection_score,
                "eigenvector_overlap_previous_a": row.eigenvector_overlap_previous_a,
                "half_width": row.half_width,
                "grid_size": row.grid_size,
                "min_log10_relative_sigma": row.min_log10_relative_sigma,
                "q10_susceptibility": -row.q10_log10_relative_sigma,
                "q50_susceptibility": -row.q50_log10_relative_sigma,
                "q10_log10_sigma_min": row.q10_log10_sigma_min,
                "q50_log10_sigma_min": row.q50_log10_sigma_min,
                "area_fraction_m10": row.area_fraction_m10,
                "area_fraction_m9_5": row.area_fraction_m9_5,
                "area_fraction_m9": row.area_fraction_m9,
                "boundary_touch_m10": row.boundary_touch_m10,
                "curvature_trace_indicator": row.curvature_trace_indicator,
                "eigen_condition_indicator": row.eigen_condition_indicator,
                "nearest_eigenvalue_distance": row.nearest_eigenvalue_distance,
                "hit_minimizer_boundary": row.hit_minimizer_boundary,
            }
        )
    _write_dict_csv(output, dict_rows)


def write_verdict_csv(output: Path, verdicts: list[BranchVerdict]) -> None:
    _write_dict_csv(output, [verdict.__dict__ for verdict in verdicts])


def write_barrier_csv(output: Path, metrics: list[BarrierMetric]) -> None:
    _write_dict_csv(output, [metric.__dict__ for metric in metrics])


def write_pair_csv(output: Path, diagnostics: list[ModePairDiagnostic]) -> None:
    _write_dict_csv(output, [diagnostic.__dict__ for diagnostic in diagnostics])


def _endpoint_gain_matrix(verdicts: list[BranchVerdict]) -> np.ndarray:
    matrix = np.full((len(PRL_OVERTONES), len(PRL_ELL_VALUES)), np.nan)
    for verdict in verdicts:
        row = PRL_OVERTONES.index(verdict.overtone)
        col = PRL_ELL_VALUES.index(verdict.ell)
        matrix[row, col] = verdict.q10_gain_n64
    return matrix


def plot_central_heatmap(output: Path, verdicts: list[BranchVerdict]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    matrix = _endpoint_gain_matrix(verdicts)
    fig, ax = plt.subplots(figsize=(7.0, 3.2), constrained_layout=True)
    vmax = max(0.25, float(np.nanmax(np.abs(matrix))))
    image = ax.imshow(matrix, cmap="coolwarm", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(PRL_ELL_VALUES)), [str(ell) for ell in PRL_ELL_VALUES])
    ax.set_yticks(range(len(PRL_OVERTONES)), [str(n) for n in PRL_OVERTONES])
    ax.set_xlabel(r"scalar multipole $\ell$")
    ax.set_ylabel("overtone n")
    for verdict in verdicts:
        row = PRL_OVERTONES.index(verdict.overtone)
        col = PRL_ELL_VALUES.index(verdict.ell)
        marker = "" if verdict.prl_support else "*"
        ax.text(col, row, f"{verdict.q10_gain_n64:+.2f}{marker}", ha="center", va="center", fontsize=8)
    ax.text(0.01, 1.04, r"$\Delta[-Q_{10}(\log_{10}\eta_N)]$, $N=64$, $a/M:0\to1$",
            transform=ax.transAxes, fontsize=9)
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label("susceptibility gain")
    ax.text(
        0.99,
        -0.18,
        "* caution/exploratory branch",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
    )
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_robustness(output: Path, rows: list[BranchScanRow]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.0), constrained_layout=True)
    for ell in PRL_ELL_VALUES:
        gains = []
        ns = []
        for spectral_n in (32, 48, 64, 96):
            n_rows = sorted(
                [
                    row
                    for row in rows
                    if row.ell == ell
                    and row.overtone == 0
                    and row.spectral_n == spectral_n
                ],
                key=lambda row: row.a,
            )
            if len(n_rows) >= 2:
                gains.append((-n_rows[-1].q10_log10_relative_sigma) - (-n_rows[0].q10_log10_relative_sigma))
                ns.append(spectral_n)
        if gains:
            ax.plot(ns, gains, marker="o", label=fr"$\ell={ell}$")
    ax.axhline(0.0, color="black", linewidth=0.8, alpha=0.4)
    ax.set_xlabel("Chebyshev size N")
    ax.set_ylabel(r"endpoint gain in $-Q_{10}(\log_{10}\eta_N)$")
    ax.text(
        0.01,
        1.03,
        "fundamental-branch finite-N robustness",
        transform=ax.transAxes,
        fontsize=9,
    )
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, ncol=3, fontsize=8)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_softening_scatter(output: Path, verdicts: list[BranchVerdict]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.2, 4.0), constrained_layout=True)
    markers = {0: "o", 1: "s", 2: "^"}
    colors = dict(zip(PRL_ELL_VALUES, plt.cm.viridis(np.linspace(0.1, 0.9, len(PRL_ELL_VALUES)))))
    for verdict in verdicts:
        ax.scatter(
            -100.0 * verdict.real_frequency_shift_n64,
            verdict.q10_gain_n64,
            marker=markers.get(verdict.overtone, "o"),
            color=colors[verdict.ell],
            s=55,
            alpha=0.85,
        )
        ax.text(
            -100.0 * verdict.real_frequency_shift_n64,
            verdict.q10_gain_n64,
            f" {verdict.ell}",
            fontsize=7,
            va="center",
        )
    ax.axhline(0.0, color="black", linewidth=0.8, alpha=0.4)
    ax.set_xlabel(r"frequency softening $-\Delta\mathrm{Re}(\omega)$ (%)")
    ax.set_ylabel(r"$\Delta[-Q_{10}(\log_{10}\eta_N)]$")
    ax.grid(alpha=0.25)
    marker_handles = [
        Line2D([0], [0], marker=markers[n], color="black", linestyle="", label=fr"$n={n}$")
        for n in PRL_OVERTONES
    ]
    color_handles = [
        Line2D([0], [0], marker="o", color=colors[ell], linestyle="", label=fr"$\ell={ell}$")
        for ell in PRL_ELL_VALUES
    ]
    first = ax.legend(handles=marker_handles, frameon=False, fontsize=8, loc="upper left")
    ax.add_artist(first)
    ax.legend(handles=color_handles, frameon=False, fontsize=8, loc="upper right", ncol=2)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def write_assessment_report(
    output: Path,
    verdicts: list[BranchVerdict],
    pair_diagnostics: list[ModePairDiagnostic],
    leaver_checks: list[dict],
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    usable = [verdict for verdict in verdicts if verdict.reliability == "usable"]
    support = [verdict for verdict in verdicts if verdict.prl_support]
    negative_gains = [verdict for verdict in verdicts if verdict.q10_gain_n64 <= 0.0]
    nonmonotonic = [verdict for verdict in verdicts if not verdict.monotonic_q10_n64]
    possible_pairs = [diag for diag in pair_diagnostics if diag.possible_coalescence]
    failed_leaver = [check for check in leaver_checks if check["status"] != "ok"]
    max_leaver_diff = max(
        [float(check["relative_difference"]) for check in leaver_checks if check["status"] == "ok"],
        default=float("nan"),
    )

    verdict_line = "not PRL-level yet"
    if usable and len(support) == len(usable) and not negative_gains and not nonmonotonic and not failed_leaver:
        verdict_line = "interesting but still not PRL-level without a continuum/operator mechanism and comparator metric"

    lines = [
        "# PRL Instability Assessment",
        "",
        f"## Short Verdict: {verdict_line}",
        "",
        "The scan tests the proposed Letter claim that KS quantum deformation softens",
        "scalar ringdown while universally amplifying finite-dimensional pseudospectral",
        "sensitivity.  The result is interesting, but the evidence is not yet strong",
        "enough for a Physical Review Letters claim.",
        "",
        "## Strongest Supported One-Sentence Claim",
        "",
        "In the current Chebyshev residual normalization, several scalar KS branches",
        "show simultaneous frequency softening and increased local finite-N",
        "pseudospectral susceptibility, but the trend is branch- and reliability-dependent",
        "rather than a demonstrated universal instability law.",
        "",
        "## Branch-Level Outcome",
        "",
        f"- Active oscillatory scalar branches scanned: {len(verdicts)}.",
        f"- Usable branches by finite-N criteria: {len(usable)}.",
        f"- Usable branches supporting positive monotonic gain: {len(support)}.",
        f"- Branches with non-positive N=64 endpoint gain: {len(negative_gains)}.",
        f"- Branches with nonmonotonic N=64 susceptibility in a/M: {len(nonmonotonic)}.",
        f"- Endpoint Leaver checks that failed to converge: {len(failed_leaver)}.",
        f"- Largest endpoint Leaver relative difference among successful checks: {max_leaver_diff:.3e}.",
        f"- Exceptional-point-like pair candidates found: {len(possible_pairs)}.",
        "",
        "## Why This Is Not A PRL Yet",
        "",
        "- The claim is still based on finite Chebyshev matrices, not a continuum",
        "  pseudospectrum theorem for the KS wave operator.",
        "- The low-ell and second-overtone sectors require caution; ell=0,n=2 is not",
        "  included as an oscillatory branch in the current selector.",
        "- Endpoint Leaver checks expose overtone risk: some checks fail to converge,",
        "  and the largest successful exploratory discrepancy is too large for a",
        "  Letter-level universal claim.",
        "- No independent quantum-corrected or regular black-hole comparator has been",
        "  implemented, so the effect cannot yet be called universal across quantum",
        "  black-hole models.",
        "- No analytical mechanism has been derived; barrier metrics are diagnostic",
        "  rather than explanatory.",
        "- No exceptional point is supported: frequency distances and mode-shape",
        "  overlaps do not demonstrate eigenvalue/eigenvector coalescence.",
        "",
        "## Publishable Direction",
        "",
        "This is strong PRD/CQG material if presented as a scalar finite-N",
        "pseudospectrum audit of KS ringdown softening.  It should not be submitted",
        "as a PRL unless a stronger mechanism, continuum robustness argument, and at",
        "least one comparator metric are added.",
        "",
        "## Branch Table",
        "",
        "| ell | n | gain | monotonic | all-N positive | reliability | comment |",
        "|---:|---:|---:|:---:|:---:|---|---|",
    ]
    for verdict in verdicts:
        lines.append(
            f"| {verdict.ell} | {verdict.overtone} | {verdict.q10_gain_n64:+.3f} | "
            f"{verdict.monotonic_q10_n64} | {verdict.positive_gain_all_tested_n} | "
            f"{verdict.reliability} | {verdict.comment} |"
        )
    output.write_text("\n".join(lines) + "\n")
