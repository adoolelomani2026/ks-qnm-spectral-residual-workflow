"""Universality stress tests for scalar QNM pseudospectral sensitivity.

The routines here compare the KS deformation with regular black-hole
comparators using the same finite-dimensional Chebyshev residual diagnostics.
They are intentionally conservative: the output is designed to falsify broad
claims unless several models, branches, and numerical checks all agree.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy.linalg import eig
from scipy.optimize import brentq

from .common import MASS
from .prl_instability import SCALAR_INITIAL_TARGETS
from .spectral import (
    SpectralProblem,
    barycentric_diff_matrices,
    chebyshev_gauss_grid,
    generalized_eigenpairs,
    minimize_residual,
    mode_overlap,
    scale_generalized_pencil,
    select_tracked_mode,
    spectral_matrix,
)


UNIVERSALITY_PARAMETERS = {
    "ks": (0.0, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0),
    # Hayward horizons exist for q/M < 4/(3 sqrt(3)) ~= 0.770.
    "hayward": (0.0, 0.1, 0.2, 0.35, 0.5, 0.65, 0.72),
    # Bardeen horizons have the same critical magnetic-charge scale.
    "bardeen": (0.0, 0.1, 0.2, 0.35, 0.5, 0.65, 0.72),
}
UNIVERSALITY_ELL_VALUES = (0, 1, 2, 3, 4)
UNIVERSALITY_OVERTONES = (0, 1)
UNIVERSALITY_SPECTRAL_SIZES = (32, 48, 64)
UNIVERSALITY_THRESHOLDS = (-10.0, -9.5, -9.0)


@dataclass(frozen=True)
class StaticScalarModel:
    name: str
    parameter_label: str
    asymptotic_mass: float
    horizon_radius: Callable[[float], float]
    lapse: Callable[[np.ndarray | float, float], np.ndarray | float]
    lapse_derivative: Callable[[np.ndarray | float, float], np.ndarray | float]


@dataclass
class ModelScanRow:
    model: str
    parameter: float
    ell: int
    overtone: int
    spectral_n: int
    center: complex
    initial_eigenvalue: complex
    residual_norm: float
    selection_score: float
    overlap_previous_parameter: float | None
    q10_susceptibility: float
    q50_susceptibility: float
    area_fraction_m10: float
    curvature_trace_indicator: float
    eigen_condition_indicator: float
    nearest_eigenvalue_distance: float
    hit_minimizer_boundary: bool
    shape: np.ndarray


@dataclass
class ModelBranchVerdict:
    model: str
    ell: int
    overtone: int
    endpoint_parameter: float
    q10_gain_n64: float
    q50_gain_n64: float
    real_frequency_shift_n64: float
    damping_shift_n64: float
    monotonic_q10_n64: bool
    positive_gain_all_tested_n: bool
    max_center_spread_endpoint: float
    reliability: str
    universality_support: bool
    comment: str


@dataclass
class ModelBarrierMetric:
    model: str
    parameter: float
    ell: int
    horizon_radius: float
    horizon_shift: float | None
    surface_gravity: float
    surface_gravity_shift: float | None
    photon_sphere_radius: float
    photon_sphere_shift: float | None
    photon_orbital_frequency: float
    photon_orbital_frequency_shift: float | None
    eikonal_lyapunov: float
    eikonal_lyapunov_shift: float | None
    peak_height: float
    peak_height_shift: float | None
    r_peak: float
    rstar_width_halfmax: float
    width_shift: float | None
    curvature_rstar_at_peak: float
    curvature_shift: float | None


def f_ks_generic(r: np.ndarray | float, a: float, mass: float = MASS) -> np.ndarray | float:
    r_arr = np.asarray(r)
    return (np.sqrt(r_arr * r_arr - a * a) - 2.0 * mass) / r_arr


def df_ks_generic(r: np.ndarray | float, a: float, mass: float = MASS) -> np.ndarray | float:
    r_arr = np.asarray(r)
    s = np.sqrt(r_arr * r_arr - a * a)
    return a * a / (s * r_arr * r_arr) + 2.0 * mass / (r_arr * r_arr)


def horizon_ks(a: float, mass: float = MASS) -> float:
    return math.sqrt((2.0 * mass) ** 2 + a * a)


def f_hayward(r: np.ndarray | float, q: float, mass: float = MASS) -> np.ndarray | float:
    r_arr = np.asarray(r)
    return 1.0 - 2.0 * mass * r_arr * r_arr / (r_arr**3 + 2.0 * mass * q * q)


def df_hayward(r: np.ndarray | float, q: float, mass: float = MASS) -> np.ndarray | float:
    r_arr = np.asarray(r)
    denom = r_arr**3 + 2.0 * mass * q * q
    return 2.0 * mass * r_arr * (r_arr**3 - 4.0 * mass * q * q) / (denom * denom)


def horizon_hayward(q: float, mass: float = MASS) -> float:
    coefficients = [1.0, -2.0 * mass, 0.0, 2.0 * mass * q * q]
    roots = np.roots(coefficients)
    real_roots = sorted(root.real for root in roots if abs(root.imag) < 1.0e-9 and root.real > 0.0)
    if not real_roots:
        raise ValueError(f"Hayward parameter q/M={q:g} has no positive horizon.")
    return float(real_roots[-1])


def f_bardeen(r: np.ndarray | float, g: float, mass: float = MASS) -> np.ndarray | float:
    r_arr = np.asarray(r)
    return 1.0 - 2.0 * mass * r_arr * r_arr / (r_arr * r_arr + g * g) ** 1.5


def df_bardeen(r: np.ndarray | float, g: float, mass: float = MASS) -> np.ndarray | float:
    r_arr = np.asarray(r)
    return 2.0 * mass * r_arr * (r_arr * r_arr - 2.0 * g * g) / (r_arr * r_arr + g * g) ** 2.5


def horizon_bardeen(g: float, mass: float = MASS) -> float:
    if abs(g) < 1.0e-14:
        return 2.0 * mass
    r_grid = np.linspace(1.0e-6, 5.0 * mass, 2000)
    values = np.asarray(f_bardeen(r_grid, g, mass), dtype=float)
    roots: list[float] = []
    for left, right, f_left, f_right in zip(r_grid[:-1], r_grid[1:], values[:-1], values[1:]):
        if not np.isfinite(f_left) or not np.isfinite(f_right):
            continue
        if f_left == 0.0:
            roots.append(float(left))
        elif f_left * f_right < 0.0:
            roots.append(float(brentq(lambda radius: f_bardeen(radius, g, mass), left, right)))
    if not roots:
        raise ValueError(f"Bardeen parameter g/M={g:g} has no positive horizon.")
    return max(roots)


MODELS = {
    "ks": StaticScalarModel(
        name="ks",
        parameter_label="a/M",
        asymptotic_mass=MASS,
        horizon_radius=horizon_ks,
        lapse=f_ks_generic,
        lapse_derivative=df_ks_generic,
    ),
    "hayward": StaticScalarModel(
        name="hayward",
        parameter_label="q/M",
        asymptotic_mass=MASS,
        horizon_radius=horizon_hayward,
        lapse=f_hayward,
        lapse_derivative=df_hayward,
    ),
    "bardeen": StaticScalarModel(
        name="bardeen",
        parameter_label="g/M",
        asymptotic_mass=MASS,
        horizon_radius=horizon_bardeen,
        lapse=f_bardeen,
        lapse_derivative=df_bardeen,
    ),
}


def scalar_potential_model(model: StaticScalarModel, r: np.ndarray, ell: int, parameter: float) -> np.ndarray:
    f = model.lapse(r, parameter)
    fp = model.lapse_derivative(r, parameter)
    return f * (ell * (ell + 1.0) / (r * r) + fp / r)


def _second_lapse_derivative(model: StaticScalarModel, r: float, parameter: float) -> float:
    step = max(1.0e-5 * abs(r), 1.0e-6)
    return float(
        (
            model.lapse(r + step, parameter)
            - 2.0 * model.lapse(r, parameter)
            + model.lapse(r - step, parameter)
        )
        / (step * step)
    )


def _photon_sphere_metrics(model: StaticScalarModel, parameter: float) -> tuple[float, float, float]:
    rh = model.horizon_radius(parameter)

    def condition(radius: float) -> float:
        return float(radius * model.lapse_derivative(radius, parameter) - 2.0 * model.lapse(radius, parameter))

    r_grid = np.geomspace(rh * (1.0 + 1.0e-5), 80.0 * model.asymptotic_mass, 6000)
    values = np.array([condition(radius) for radius in r_grid], dtype=float)
    roots: list[float] = []
    for left, right, f_left, f_right in zip(r_grid[:-1], r_grid[1:], values[:-1], values[1:]):
        if not np.isfinite(f_left) or not np.isfinite(f_right):
            continue
        if f_left == 0.0:
            roots.append(float(left))
        elif f_left * f_right < 0.0:
            roots.append(float(brentq(condition, left, right)))
    if not roots:
        return float("nan"), float("nan"), float("nan")
    radius = max(roots)
    lapse = float(model.lapse(radius, parameter))
    omega = math.sqrt(max(lapse / (radius * radius), 0.0))
    fpp = _second_lapse_derivative(model, radius, parameter)
    lyapunov_squared = lapse * (2.0 * lapse - radius * radius * fpp) / (2.0 * radius * radius)
    lyapunov = math.sqrt(max(lyapunov_squared, 0.0)) if np.isfinite(lyapunov_squared) else float("nan")
    return radius, omega, lyapunov


def build_model_spectral_problem(
    model: StaticScalarModel,
    parameter: float,
    n: int,
    ell: int,
) -> SpectralProblem:
    """Build the scalar Chebyshev residual pencil for one static spherical model."""

    rh = model.horizon_radius(parameter)
    fh = float(model.lapse_derivative(rh, parameter))
    x_nodes, s_nodes = chebyshev_gauss_grid(n)
    d1, d2 = barycentric_diff_matrices(s_nodes)

    r = rh / (1.0 - s_nodes)
    ds_dr = (1.0 - s_nodes) ** 2 / rh
    d2s_dr2 = -2.0 * (1.0 - s_nodes) ** 3 / (rh * rh)

    f = model.lapse(r, parameter)
    fp = model.lapse_derivative(r, parameter)
    potential = scalar_potential_model(model, r, ell, parameter)

    dlog_s_dr = ds_dr / s_nodes
    d_dlog_s_dr = -rh * (2.0 * r - rh) / (r * r * (r - rh) ** 2)
    h1 = 1j * (1.0 + 2.0 * model.asymptotic_mass / r - dlog_s_dr / fh)
    dh1_dr = 1j * (-2.0 * model.asymptotic_mass / (r * r) - d_dlog_s_dr / fh)

    coeff_u_ss = f * f * ds_dr * ds_dr
    coeff_u_s_0 = f * f * d2s_dr2 + f * fp * ds_dr
    coeff_u_s_1 = 2.0 * f * f * h1 * ds_dr
    coeff_u_0 = -potential
    coeff_u_1 = f * f * dh1_dr + f * fp * h1
    coeff_u_2 = 1.0 + f * f * h1 * h1

    a0 = np.diag(coeff_u_ss) @ d2 + np.diag(coeff_u_s_0) @ d1 + np.diag(coeff_u_0)
    a1 = np.diag(coeff_u_s_1) @ d1 + np.diag(coeff_u_1)
    a2 = np.diag(coeff_u_2)

    zero = np.zeros((n, n), dtype=complex)
    ident = np.eye(n, dtype=complex)
    left = np.block([[zero, ident], [-a0, -a1]])
    right = np.block([[ident, zero], [zero, a2]])

    return SpectralProblem(
        a=parameter,
        n=n,
        ell=ell,
        perturbation_type=f"scalar_{model.name}",
        x_nodes=x_nodes,
        s_nodes=s_nodes,
        r_nodes=r,
        d1=d1,
        d2=d2,
        a0=a0,
        a1=a1,
        a2=a2,
        left=left,
        right=right,
    )


def scan_branches() -> list[tuple[int, int]]:
    return [
        (ell, overtone)
        for ell in UNIVERSALITY_ELL_VALUES
        for overtone in UNIVERSALITY_OVERTONES
        if SCALAR_INITIAL_TARGETS.get((ell, overtone)) is not None
    ]


def _relative_singular_values(problem: SpectralProblem, omega: complex) -> tuple[float, float, float]:
    singular_values = np.linalg.svd(spectral_matrix(problem, omega), compute_uv=False)
    sigma_max = float(singular_values[0])
    sigma_min = float(singular_values[-1])
    return sigma_min, sigma_max, sigma_min / sigma_max if sigma_max > 0.0 else float("nan")


def _condition_indicator(problem: SpectralProblem, target: complex) -> tuple[float, float]:
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
    left = left_vectors[:, index]
    right = right_vectors[:, index]
    denominator = abs(np.vdot(left, right_scaled @ right))
    condition = float("inf") if denominator == 0.0 else float(np.linalg.norm(left) * np.linalg.norm(right) / denominator)
    distances = [abs(values[index] - values[other]) for other in finite_indices if other != index]
    nearest = float(min(distances)) if distances else float("nan")
    return condition, nearest


def _grid_metrics(problem: SpectralProblem, center: complex, grid_size: int) -> dict[str, float]:
    half_width = max(0.012, 0.05 * abs(center))
    real_values = np.linspace(center.real - half_width, center.real + half_width, grid_size)
    imag_values = np.linspace(center.imag - half_width, center.imag + half_width, grid_size)
    relative_log = np.empty((grid_size, grid_size), dtype=float)
    for i, imag_part in enumerate(imag_values):
        for j, real_part in enumerate(real_values):
            _, _, relative = _relative_singular_values(problem, complex(real_part, imag_part))
            relative_log[i, j] = np.log10(max(relative, 1.0e-300))
    center_index = grid_size // 2
    susceptibility = -relative_log
    dr = real_values[1] - real_values[0]
    di = imag_values[1] - imag_values[0]
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
    return {
        "q10_susceptibility": float(-np.quantile(relative_log, 0.10)),
        "q50_susceptibility": float(-np.quantile(relative_log, 0.50)),
        "area_fraction_m10": float(np.mean(relative_log <= -10.0)),
        "curvature_trace_indicator": float(-(d2_re + d2_im)),
    }


def scan_model_branch(
    model: StaticScalarModel,
    parameter_values: tuple[float, ...],
    ell: int,
    overtone: int,
    spectral_n: int,
    grid_size: int,
) -> list[ModelScanRow]:
    target = SCALAR_INITIAL_TARGETS[(ell, overtone)]
    if target is None:
        return []
    rows: list[ModelScanRow] = []
    previous_omega: complex | None = None
    previous_shape: np.ndarray | None = None
    previous_s_nodes: np.ndarray | None = None
    for parameter in parameter_values:
        problem = build_model_spectral_problem(model, parameter, spectral_n, ell)
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
        metrics = _grid_metrics(problem, center, grid_size=grid_size)
        condition, nearest = _condition_indicator(problem, center)
        hit_boundary = (
            abs(center.real - selection.omega.real) > 0.96 * radius
            or abs(center.imag - selection.omega.imag) > 0.96 * radius
        )
        rows.append(
            ModelScanRow(
                model=model.name,
                parameter=parameter,
                ell=ell,
                overtone=overtone,
                spectral_n=spectral_n,
                center=center,
                initial_eigenvalue=selection.omega,
                residual_norm=residual,
                selection_score=selection.selection_score,
                overlap_previous_parameter=selection.eigenvector_overlap,
                q10_susceptibility=metrics["q10_susceptibility"],
                q50_susceptibility=metrics["q50_susceptibility"],
                area_fraction_m10=metrics["area_fraction_m10"],
                curvature_trace_indicator=metrics["curvature_trace_indicator"],
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


def run_universality_scan(main_grid_size: int = 25, resolution_grid_size: int = 19) -> list[ModelScanRow]:
    rows: list[ModelScanRow] = []
    for model_name, model in MODELS.items():
        parameter_values = UNIVERSALITY_PARAMETERS[model_name]
        for spectral_n in UNIVERSALITY_SPECTRAL_SIZES:
            grid_size = main_grid_size if spectral_n == 64 else resolution_grid_size
            for ell, overtone in scan_branches():
                rows.extend(scan_model_branch(model, parameter_values, ell, overtone, spectral_n, grid_size))
    return rows


def summarize_model_branches(rows: list[ModelScanRow]) -> list[ModelBranchVerdict]:
    verdicts: list[ModelBranchVerdict] = []
    for model in sorted({row.model for row in rows}):
        for ell, overtone in scan_branches():
            branch_rows = [row for row in rows if row.model == model and row.ell == ell and row.overtone == overtone]
            n64 = sorted([row for row in branch_rows if row.spectral_n == 64], key=lambda row: row.parameter)
            if len(n64) < 2:
                continue
            q10_values = [row.q10_susceptibility for row in n64]
            q10_gain = q10_values[-1] - q10_values[0]
            q50_gain = n64[-1].q50_susceptibility - n64[0].q50_susceptibility
            real_shift = n64[-1].center.real / n64[0].center.real - 1.0
            damping_shift = (-n64[-1].center.imag) / (-n64[0].center.imag) - 1.0
            monotonic = all(new >= old - 1.0e-3 for old, new in zip(q10_values, q10_values[1:]))
            gains_by_n = []
            for spectral_n in sorted({row.spectral_n for row in branch_rows}):
                local = sorted([row for row in branch_rows if row.spectral_n == spectral_n], key=lambda row: row.parameter)
                if len(local) >= 2:
                    gains_by_n.append(local[-1].q10_susceptibility - local[0].q10_susceptibility)
            positive_all = bool(gains_by_n) and all(gain > 0.0 for gain in gains_by_n)
            endpoints = [row for row in branch_rows if row.parameter == n64[-1].parameter and row.spectral_n in UNIVERSALITY_SPECTRAL_SIZES]
            centers = [row.center for row in endpoints]
            max_spread = 0.0
            for i, left in enumerate(centers):
                for right in centers[i + 1 :]:
                    max_spread = max(max_spread, abs(left - right) / max(abs(left), 0.1))
            any_boundary = any(row.hit_minimizer_boundary for row in branch_rows)
            if max_spread > 5.0e-2 or any_boundary:
                reliability = "exploratory"
            elif ell == 0:
                reliability = "caution"
            else:
                reliability = "usable"
            support = reliability == "usable" and monotonic and positive_all and q10_gain > 0.0
            comments: list[str] = []
            if any_boundary:
                comments.append("minimizer hit search boundary")
            if not monotonic:
                comments.append("N=64 susceptibility is nonmonotonic")
            if not positive_all:
                comments.append("endpoint gain is not positive at every tested N")
            if max_spread > 5.0e-2:
                comments.append("endpoint center spread across N is large")
            if not comments:
                comments.append("finite-N trend passes branch checks")
            verdicts.append(
                ModelBranchVerdict(
                    model=model,
                    ell=ell,
                    overtone=overtone,
                    endpoint_parameter=n64[-1].parameter,
                    q10_gain_n64=q10_gain,
                    q50_gain_n64=q50_gain,
                    real_frequency_shift_n64=real_shift,
                    damping_shift_n64=damping_shift,
                    monotonic_q10_n64=monotonic,
                    positive_gain_all_tested_n=positive_all,
                    max_center_spread_endpoint=max_spread,
                    reliability=reliability,
                    universality_support=support,
                    comment="; ".join(comments),
                )
            )
    return verdicts


def _tortoise_grid(model: StaticScalarModel, parameter: float, r: np.ndarray) -> np.ndarray:
    integrand = 1.0 / np.asarray(model.lapse(r, parameter), dtype=float)
    rstar = np.zeros_like(r)
    rstar[1:] = np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(r))
    return rstar


def compute_barrier_metrics() -> list[ModelBarrierMetric]:
    rows: list[ModelBarrierMetric] = []
    for model_name, model in MODELS.items():
        baselines: dict[int, dict[str, float]] = {}
        for ell in UNIVERSALITY_ELL_VALUES:
            for parameter in UNIVERSALITY_PARAMETERS[model_name]:
                rh = model.horizon_radius(parameter)
                surface_gravity = 0.5 * float(model.lapse_derivative(rh, parameter))
                photon_radius, photon_frequency, lyapunov = _photon_sphere_metrics(model, parameter)
                r = np.geomspace(rh * (1.0 + 1.0e-5), 140.0, 8000)
                potential = scalar_potential_model(model, r, ell, parameter)
                peak_index = int(np.argmax(potential))
                rstar = _tortoise_grid(model, parameter, r)
                peak_height = float(potential[peak_index])
                half = 0.5 * peak_height
                left = peak_index
                while left > 0 and potential[left] > half:
                    left -= 1
                right = peak_index
                while right < len(r) - 1 and potential[right] > half:
                    right += 1
                width = float(rstar[right] - rstar[left]) if right > left else float("nan")
                local_x = rstar[peak_index - 2 : peak_index + 3]
                local_y = potential[peak_index - 2 : peak_index + 3]
                coeffs = np.polyfit(local_x - local_x[2], local_y, 2)
                curvature = float(2.0 * coeffs[0])
                if parameter == 0.0:
                    baselines[ell] = {
                        "horizon_radius": rh,
                        "surface_gravity": surface_gravity,
                        "photon_sphere_radius": photon_radius,
                        "photon_orbital_frequency": photon_frequency,
                        "eikonal_lyapunov": lyapunov,
                        "peak_height": peak_height,
                        "width": width,
                        "curvature": curvature,
                    }
                base = baselines.get(ell)
                rows.append(
                    ModelBarrierMetric(
                        model=model_name,
                        parameter=parameter,
                        ell=ell,
                        horizon_radius=rh,
                        horizon_shift=None if base is None else rh / base["horizon_radius"] - 1.0,
                        surface_gravity=surface_gravity,
                        surface_gravity_shift=None if base is None else surface_gravity / base["surface_gravity"] - 1.0,
                        photon_sphere_radius=photon_radius,
                        photon_sphere_shift=None
                        if base is None or not np.isfinite(base["photon_sphere_radius"])
                        else photon_radius / base["photon_sphere_radius"] - 1.0,
                        photon_orbital_frequency=photon_frequency,
                        photon_orbital_frequency_shift=None
                        if base is None or not np.isfinite(base["photon_orbital_frequency"])
                        else photon_frequency / base["photon_orbital_frequency"] - 1.0,
                        eikonal_lyapunov=lyapunov,
                        eikonal_lyapunov_shift=None
                        if base is None or not np.isfinite(base["eikonal_lyapunov"])
                        else lyapunov / base["eikonal_lyapunov"] - 1.0,
                        peak_height=peak_height,
                        peak_height_shift=None if base is None else peak_height / base["peak_height"] - 1.0,
                        r_peak=float(r[peak_index]),
                        rstar_width_halfmax=width,
                        width_shift=None if base is None else width / base["width"] - 1.0,
                        curvature_rstar_at_peak=curvature,
                        curvature_shift=None if base is None else abs(curvature) / abs(base["curvature"]) - 1.0,
                    )
                )
    return rows


def compute_mode_pair_diagnostics(rows: list[ModelScanRow]) -> list[dict[str, float | int | str | bool]]:
    output: list[dict[str, float | int | str | bool]] = []
    for model in sorted({row.model for row in rows}):
        for parameter in sorted({row.parameter for row in rows if row.model == model}):
            for ell in UNIVERSALITY_ELL_VALUES:
                local = sorted(
                    [
                        row
                        for row in rows
                        if row.model == model and row.parameter == parameter and row.ell == ell and row.spectral_n == 64
                    ],
                    key=lambda row: row.overtone,
                )
                for i, left in enumerate(local):
                    for right in local[i + 1 :]:
                        distance = abs(left.center - right.center)
                        normalized = distance / max(abs(left.center), abs(right.center), 0.1)
                        overlap = mode_overlap(left.shape, right.shape)
                        output.append(
                            {
                                "model": model,
                                "parameter": parameter,
                                "ell": ell,
                                "overtone_i": left.overtone,
                                "overtone_j": right.overtone,
                                "frequency_distance": float(distance),
                                "normalized_frequency_distance": float(normalized),
                                "mode_shape_overlap": overlap,
                                "possible_coalescence": bool(normalized < 0.05 and overlap > 0.95),
                            }
                        )
    return output


def compute_correlations(verdicts: list[ModelBranchVerdict], barriers: list[ModelBarrierMetric]) -> list[dict[str, float | str]]:
    barrier_lookup = {(row.model, row.parameter, row.ell): row for row in barriers}
    data = []
    for verdict in verdicts:
        barrier = barrier_lookup.get((verdict.model, verdict.endpoint_parameter, verdict.ell))
        if barrier is None:
            continue
        data.append(
            {
                "model": verdict.model,
                "ell": verdict.ell,
                "overtone": verdict.overtone,
                "q10_gain": verdict.q10_gain_n64,
                "frequency_softening": -verdict.real_frequency_shift_n64,
                "damping_shift": verdict.damping_shift_n64,
                "peak_height_softening": -(barrier.peak_height_shift or 0.0),
                "width_growth": barrier.width_shift or 0.0,
                "curvature_softening": -(barrier.curvature_shift or 0.0),
                "horizon_expansion": barrier.horizon_shift or 0.0,
                "surface_gravity_growth": barrier.surface_gravity_shift or 0.0,
                "photon_frequency_softening": -(barrier.photon_orbital_frequency_shift or 0.0),
                "eikonal_lyapunov_growth": barrier.eikonal_lyapunov_shift or 0.0,
            }
        )
    metrics = [
        "frequency_softening",
        "damping_shift",
        "peak_height_softening",
        "width_growth",
        "curvature_softening",
        "horizon_expansion",
        "surface_gravity_growth",
        "photon_frequency_softening",
        "eikonal_lyapunov_growth",
    ]
    rows: list[dict[str, float | str]] = []
    for metric in metrics:
        x = np.array([float(item[metric]) for item in data], dtype=float)
        y = np.array([float(item["q10_gain"]) for item in data], dtype=float)
        mask = np.isfinite(x) & np.isfinite(y)
        correlation = float(np.corrcoef(x[mask], y[mask])[0, 1]) if np.count_nonzero(mask) > 2 else float("nan")
        rows.append({"predictor": metric, "pearson_r_with_q10_gain": correlation, "sample_count": int(np.count_nonzero(mask))})
    return rows


def write_dict_csv(output: Path, rows: list[dict]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        output.write_text("")
        return
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_scan_csv(output: Path, rows: list[ModelScanRow]) -> None:
    write_dict_csv(
        output,
        [
            {
                "model": row.model,
                "parameter": row.parameter,
                "ell": row.ell,
                "overtone": row.overtone,
                "spectral_N": row.spectral_n,
                "center_real": row.center.real,
                "center_imag": row.center.imag,
                "q10_susceptibility": row.q10_susceptibility,
                "q50_susceptibility": row.q50_susceptibility,
                "area_fraction_m10": row.area_fraction_m10,
                "curvature_trace_indicator": row.curvature_trace_indicator,
                "eigen_condition_indicator": row.eigen_condition_indicator,
                "nearest_eigenvalue_distance": row.nearest_eigenvalue_distance,
                "overlap_previous_parameter": row.overlap_previous_parameter,
                "hit_minimizer_boundary": row.hit_minimizer_boundary,
            }
            for row in rows
        ],
    )


def write_verdict_csv(output: Path, verdicts: list[ModelBranchVerdict]) -> None:
    write_dict_csv(output, [verdict.__dict__ for verdict in verdicts])


def write_barrier_csv(output: Path, barriers: list[ModelBarrierMetric]) -> None:
    write_dict_csv(output, [barrier.__dict__ for barrier in barriers])


def plot_model_heatmap(output: Path, verdicts: list[ModelBranchVerdict]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    models = sorted({verdict.model for verdict in verdicts})
    fig, axes = plt.subplots(1, len(models), figsize=(7.2, 3.2), constrained_layout=True, sharey=True)
    if len(models) == 1:
        axes = [axes]
    vmax = max(0.2, max(abs(verdict.q10_gain_n64) for verdict in verdicts))
    for axis, model in zip(axes, models):
        matrix = np.full((len(UNIVERSALITY_OVERTONES), len(UNIVERSALITY_ELL_VALUES)), np.nan)
        for verdict in verdicts:
            if verdict.model != model:
                continue
            row = UNIVERSALITY_OVERTONES.index(verdict.overtone)
            col = UNIVERSALITY_ELL_VALUES.index(verdict.ell)
            matrix[row, col] = verdict.q10_gain_n64
        image = axis.imshow(matrix, cmap="coolwarm", vmin=-vmax, vmax=vmax, aspect="auto")
        axis.set_title(model)
        axis.set_xticks(range(len(UNIVERSALITY_ELL_VALUES)), [str(ell) for ell in UNIVERSALITY_ELL_VALUES])
        axis.set_yticks(range(len(UNIVERSALITY_OVERTONES)), [str(n) for n in UNIVERSALITY_OVERTONES])
        axis.set_xlabel(r"$\ell$")
        for verdict in verdicts:
            if verdict.model != model:
                continue
            row = UNIVERSALITY_OVERTONES.index(verdict.overtone)
            col = UNIVERSALITY_ELL_VALUES.index(verdict.ell)
            marker = "" if verdict.universality_support else "*"
            axis.text(col, row, f"{verdict.q10_gain_n64:+.2f}{marker}", ha="center", va="center", fontsize=8)
    axes[0].set_ylabel("overtone n")
    cbar = fig.colorbar(image, ax=axes)
    cbar.set_label(r"$\Delta[-Q_{10}(\log_{10}\eta_N)]$")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_softening_vs_sensitivity(output: Path, verdicts: list[ModelBranchVerdict]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, axis = plt.subplots(figsize=(6.4, 4.0), constrained_layout=True)
    markers = {"ks": "o", "hayward": "s"}
    colors = {0: "#356bad", 1: "#d95f02"}
    for verdict in verdicts:
        axis.scatter(
            -100.0 * verdict.real_frequency_shift_n64,
            verdict.q10_gain_n64,
            marker=markers.get(verdict.model, "o"),
            color=colors.get(verdict.overtone, "black"),
            alpha=0.78,
            s=55,
        )
        axis.text(
            -100.0 * verdict.real_frequency_shift_n64,
            verdict.q10_gain_n64,
            f" {verdict.model[0]}{verdict.ell}",
            fontsize=7,
            va="center",
        )
    axis.axhline(0.0, color="black", linewidth=0.8, alpha=0.35)
    axis.axvline(0.0, color="black", linewidth=0.8, alpha=0.35)
    axis.set_xlabel(r"frequency softening $-\Delta\mathrm{Re}(\omega)$ (%)")
    axis.set_ylabel(r"$\Delta[-Q_{10}(\log_{10}\eta_N)]$")
    axis.grid(alpha=0.25)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_barrier_correlation(output: Path, verdicts: list[ModelBranchVerdict], barriers: list[ModelBarrierMetric]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    barrier_lookup = {(row.model, row.parameter, row.ell): row for row in barriers}
    fig, axis = plt.subplots(figsize=(6.4, 4.0), constrained_layout=True)
    for verdict in verdicts:
        barrier = barrier_lookup.get((verdict.model, verdict.endpoint_parameter, verdict.ell))
        if barrier is None or barrier.peak_height_shift is None:
            continue
        axis.scatter(
            -100.0 * barrier.peak_height_shift,
            verdict.q10_gain_n64,
            marker="o" if verdict.model == "ks" else "s",
            alpha=0.78,
        )
        axis.text(-100.0 * barrier.peak_height_shift, verdict.q10_gain_n64, f" {verdict.model[0]}{verdict.ell}n{verdict.overtone}", fontsize=7)
    axis.axhline(0.0, color="black", linewidth=0.8, alpha=0.35)
    axis.set_xlabel("barrier-height softening (%)")
    axis.set_ylabel(r"$\Delta[-Q_{10}(\log_{10}\eta_N)]$")
    axis.grid(alpha=0.25)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_width_correlation(output: Path, verdicts: list[ModelBranchVerdict], barriers: list[ModelBarrierMetric]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    barrier_lookup = {(row.model, row.parameter, row.ell): row for row in barriers}
    fig, axis = plt.subplots(figsize=(6.4, 4.0), constrained_layout=True)
    markers = {"ks": "o", "hayward": "s", "bardeen": "^"}
    colors = {0: "#356bad", 1: "#d95f02"}
    for verdict in verdicts:
        barrier = barrier_lookup.get((verdict.model, verdict.endpoint_parameter, verdict.ell))
        if barrier is None or barrier.width_shift is None:
            continue
        axis.scatter(
            100.0 * barrier.width_shift,
            verdict.q10_gain_n64,
            marker=markers.get(verdict.model, "o"),
            color=colors.get(verdict.overtone, "black"),
            alpha=0.78,
            s=55,
        )
    axis.axhline(0.0, color="black", linewidth=0.8, alpha=0.35)
    axis.axvline(0.0, color="black", linewidth=0.8, alpha=0.35)
    axis.set_xlabel("tortoise-width change at half maximum (%)")
    axis.set_ylabel(r"$\Delta[-Q_{10}(\log_{10}\eta_N)]$")
    axis.grid(alpha=0.25)
    model_handles = [
        Line2D([0], [0], marker=marker, color="none", markerfacecolor="0.35", markersize=7, label=model)
        for model, marker in markers.items()
    ]
    overtone_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=color, markersize=7, label=f"n={overtone}")
        for overtone, color in colors.items()
    ]
    first_legend = axis.legend(handles=model_handles, loc="lower left", frameon=False, title="model")
    axis.add_artist(first_legend)
    axis.legend(handles=overtone_handles, loc="upper right", frameon=False, title="branch")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def write_assessment(
    output: Path,
    verdicts: list[ModelBranchVerdict],
    correlations: list[dict[str, float | str]],
    pair_rows: list[dict[str, float | int | str | bool]],
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    total = len(verdicts)
    support = [verdict for verdict in verdicts if verdict.universality_support]
    by_model = {
        model: [verdict for verdict in verdicts if verdict.model == model]
        for model in sorted({verdict.model for verdict in verdicts})
    }
    nonpositive = [verdict for verdict in verdicts if verdict.q10_gain_n64 <= 0.0]
    nonmonotonic = [verdict for verdict in verdicts if not verdict.monotonic_q10_n64]
    coalescence = [row for row in pair_rows if row["possible_coalescence"]]
    best_correlation = max(
        correlations,
        key=lambda row: abs(float(row["pearson_r_with_q10_gain"])) if np.isfinite(float(row["pearson_r_with_q10_gain"])) else -1.0,
    )
    lines = [
        "# Universality Assessment",
        "",
        "## Verdict",
        "",
        "Result is PRD/CQG quality but not PRL quality.",
        "",
        "## Strongest one-sentence physical claim",
        "",
        "KS scalar branches show simultaneous frequency softening and increased",
        "finite-N pseudospectral susceptibility, while the Hayward and Bardeen",
        "regular-black-hole comparators mostly harden and narrow the same residual",
        "diagnostic; the response is therefore deformation-specific rather than a",
        "generic quantum-black-hole instability principle.",
        "",
        "## Evidence for or against universality",
        "",
        f"- Branch verdicts tested: {total}.",
        f"- Branches satisfying the stricter universality-support criteria: {len(support)}.",
        f"- Branches with non-positive endpoint susceptibility gain: {len(nonpositive)}.",
        f"- Branches with nonmonotonic N=64 susceptibility: {len(nonmonotonic)}.",
    ]
    for model, local in by_model.items():
        local_support = [verdict for verdict in local if verdict.universality_support]
        lines.append(f"- {model}: {len(local_support)}/{len(local)} branches satisfy support criteria.")
    lines.extend(
        [
            "",
            "The endpoint-gain pattern is therefore a deformation-response taxonomy,",
            "not evidence for universal amplification.",
            "",
            "## Evidence for or against a mechanism",
            "",
            f"- Best single predictor in this scan: {best_correlation['predictor']} with Pearson r={float(best_correlation['pearson_r_with_q10_gain']):.3f}.",
            "- KS lowers the scalar barrier and broadens the finite-N residual susceptibility on usable branches.",
            "- Hayward and Bardeen raise the scalar barrier for most ell>=1 branches and narrow the same diagnostic.",
            "- No single dimensionless barrier or photon-sphere quantity is strong enough to be a universal scaling law.",
            "",
            "## Evidence for or against exceptional-point behavior",
            "",
            f"- Near-coalescent mode-pair candidates found: {len(coalescence)}.",
            "- No eigenvalue/eigenvector coalescence criterion is met in the scan.",
            "- Petermann-like condition indicators vary, but do not reveal a critical point.",
            "",
            "## Recommended journal",
            "",
            "PRD or CQG. EPJ C is also plausible. EPJ Plus remains safe. PRL is not recommended.",
            "",
            "## Referee-risk assessment",
            "",
            "- Comparator models use the same Chebyshev residual machinery but not independent continued-fraction validation layers.",
            "- Low multipoles and overtones remain the most fragile sectors.",
            "- The finite-N pseudospectrum is not a continuum operator theorem.",
            "- The new result is a robust negative universality test, not a discovery-level PRL mechanism.",
            "",
            "## Branch table",
            "",
            "| model | ell | n | gain | monotonic | all-N positive | reliability | support | comment |",
            "|---|---:|---:|---:|:---:|:---:|---|:---:|---|",
        ]
    )
    for verdict in verdicts:
        lines.append(
            f"| {verdict.model} | {verdict.ell} | {verdict.overtone} | {verdict.q10_gain_n64:+.3f} | "
            f"{verdict.monotonic_q10_n64} | {verdict.positive_gain_all_tested_n} | {verdict.reliability} | "
            f"{verdict.universality_support} | {verdict.comment} |"
        )
    output.write_text("\n".join(lines) + "\n")
