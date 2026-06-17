"""Direct Chebyshev pseudospectral QNM residual operator."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable

import numpy as np
from scipy.linalg import eig
from scipy.interpolate import BarycentricInterpolator
from scipy.optimize import minimize

from .common import (
    CATALOGUE_SPECTRAL_N,
    ELL,
    FINAL_SPECTRAL_N,
    MASS,
    SCHWARZSCHILD_SCALAR_L2,
    SCHWARZSCHILD_SCALAR_L2_OVERTONE_ESTIMATE,
    df_ks,
    f_ks,
    horizon_radius,
    perturbation_potential,
    select_physical_mode,
)


OVERTONE_PUBLICATION_N = CATALOGUE_SPECTRAL_N
OVERTONE_EXPLORATORY_N = FINAL_SPECTRAL_N
PENCIL_SCALE_FLOOR = 1.0e-300
PENCIL_SCALE_LIMIT = 1.0e100
CONDITION_WARNING_THRESHOLD = 1.0e14


@dataclass
class SpectralProblem:
    a: float
    n: int
    ell: int
    perturbation_type: str
    x_nodes: np.ndarray
    s_nodes: np.ndarray
    r_nodes: np.ndarray
    d1: np.ndarray
    d2: np.ndarray
    a0: np.ndarray
    a1: np.ndarray
    a2: np.ndarray
    left: np.ndarray
    right: np.ndarray


@dataclass
class ModeResult:
    a: float
    n: int
    mode: str
    omega: complex
    omega_residual: complex
    residual_norm: float
    relative_change: float | None
    hermiticity_error: float
    psd_min_eigenvalue: float
    matrix_dimension: int
    sparsity: float
    condition_number: float
    conditioning_warning: bool
    selection_score: float | None = None
    eigenvector_overlap: float | None = None
    branch_status: str = "tracked"
    ell: int = ELL
    perturbation_type: str = "scalar"


@dataclass
class EigenSelection:
    omega: complex
    vector: np.ndarray
    shape: np.ndarray
    residual_norm: float
    selection_score: float
    eigenvector_overlap: float | None


@dataclass
class TestResult:
    name: str
    value: float
    threshold: float
    passed: bool


def chebyshev_gauss_grid(n: int) -> tuple[np.ndarray, np.ndarray]:
    """Chebyshev-Gauss grid on x in (-1,1) and s in (0,1)."""

    j = np.arange(n)
    x = np.sort(np.cos(np.pi * (j + 0.5) / n))
    s = 0.5 * (x + 1.0)
    return x, s


def barycentric_diff_matrices(nodes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(nodes, dtype=float)
    n = len(x)
    weights = np.ones(n)
    for j in range(n):
        weights[j] = 1.0 / np.prod(x[j] - np.delete(x, j))

    d1 = np.empty((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            if i != j:
                d1[i, j] = weights[j] / weights[i] / (x[i] - x[j])
    d1[np.diag_indices(n)] = 0.0
    d1[np.diag_indices(n)] = -d1.sum(axis=1)
    return d1, d1 @ d1


def build_spectral_problem(
    a: float,
    n: int,
    ell: int = ELL,
    mass: float = MASS,
    perturbation_type: str = "scalar",
) -> SpectralProblem:
    """Build P_N(omega)=A0+omega A1+omega^2 A2 from the perturbation equation.

    The compactification is x in [-1,1], s=(x+1)/2, r=r_h/(1-s).
    QNM boundary behavior is imposed by factoring

        Psi = exp(i omega r) r^(2 i M omega) s^(-i omega/f'(r_h)) u(s).

    Regularity of u on Chebyshev-Gauss points supplies the ingoing/outgoing
    boundary behavior without placing collocation points directly on either
    singular endpoint.
    """

    rh = horizon_radius(a, mass)
    fh = float(df_ks(rh, a, mass))
    x_nodes, s_nodes = chebyshev_gauss_grid(n)
    d1, d2 = barycentric_diff_matrices(s_nodes)

    r = rh / (1.0 - s_nodes)
    ds_dr = (1.0 - s_nodes) ** 2 / rh
    d2s_dr2 = -2.0 * (1.0 - s_nodes) ** 3 / (rh * rh)

    f = f_ks(r, a, mass)
    fp = df_ks(r, a, mass)
    potential = perturbation_potential(r, ell, a, perturbation_type, mass)

    dlog_s_dr = ds_dr / s_nodes
    d_dlog_s_dr = -rh * (2.0 * r - rh) / (r * r * (r - rh) ** 2)
    h1 = 1j * (1.0 + 2.0 * mass / r - dlog_s_dr / fh)
    dh1_dr = 1j * (-2.0 * mass / (r * r) - d_dlog_s_dr / fh)

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
        a=a,
        n=n,
        ell=ell,
        perturbation_type=perturbation_type,
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


def spectral_matrix(problem: SpectralProblem, omega: complex) -> np.ndarray:
    return problem.a0 + omega * problem.a1 + omega * omega * problem.a2


def scale_generalized_pencil(left: np.ndarray, right: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Apply common row/column equilibration to a generalized pencil.

    If A x = lambda B x, then D_l A D_r y = lambda D_l B D_r y has the
    same eigenvalues and original eigenvectors x = D_r y. The scaling is
    deliberately simple and bounded; it reduces avoidable norm imbalance
    without changing the discretized problem.
    """

    row_norm = np.maximum(np.max(np.abs(left), axis=1), np.max(np.abs(right), axis=1))
    row_scale = np.ones_like(row_norm, dtype=float)
    mask = row_norm > PENCIL_SCALE_FLOOR
    row_scale[mask] = 1.0 / row_norm[mask]
    row_scale = np.clip(row_scale, 1.0 / PENCIL_SCALE_LIMIT, PENCIL_SCALE_LIMIT)

    left_scaled = row_scale[:, None] * left
    right_scaled = row_scale[:, None] * right

    col_norm = np.maximum(np.max(np.abs(left_scaled), axis=0), np.max(np.abs(right_scaled), axis=0))
    col_scale = np.ones_like(col_norm, dtype=float)
    mask = col_norm > PENCIL_SCALE_FLOOR
    col_scale[mask] = 1.0 / col_norm[mask]
    col_scale = np.clip(col_scale, 1.0 / PENCIL_SCALE_LIMIT, PENCIL_SCALE_LIMIT)

    return left_scaled * col_scale[None, :], right_scaled * col_scale[None, :], col_scale


def generalized_eigenpairs(problem: SpectralProblem, scale: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """Return generalized eigenvalues and right eigenvectors.

    The default common pencil scaling preserves the eigenvalues while making
    branch selection less vulnerable to raw matrix norm imbalance.
    """

    if not scale:
        values, vectors = eig(problem.left, problem.right, left=False, right=True)
        return values, vectors

    left_scaled, right_scaled, col_scale = scale_generalized_pencil(problem.left, problem.right)
    values, scaled_vectors = eig(left_scaled, right_scaled, left=False, right=True)
    vectors = col_scale[:, None] * scaled_vectors
    return values, vectors


def generalized_eigenvalues(problem: SpectralProblem, scale: bool = True) -> np.ndarray:
    values, _ = generalized_eigenpairs(problem, scale=scale)
    return values


def normalized_mode_shape(problem: SpectralProblem, vector: np.ndarray) -> np.ndarray:
    """Extract, normalize, and phase-fix the collocated regular mode shape."""

    shape = np.asarray(vector[: problem.n], dtype=complex)
    norm = np.linalg.norm(shape)
    if norm == 0.0 or not np.isfinite(norm):
        return shape
    shape = shape / norm
    pivot = int(np.argmax(np.abs(shape)))
    if abs(shape[pivot]) > 0.0:
        shape = shape * np.exp(-1.0j * np.angle(shape[pivot]))
    return shape


def interpolate_mode_shape(
    old_s_nodes: np.ndarray,
    old_shape: np.ndarray,
    new_s_nodes: np.ndarray,
) -> np.ndarray:
    interpolator = BarycentricInterpolator(old_s_nodes, old_shape)
    shape = np.asarray(interpolator(new_s_nodes), dtype=complex)
    norm = np.linalg.norm(shape)
    if norm == 0.0 or not np.isfinite(norm):
        return shape
    shape = shape / norm
    pivot = int(np.argmax(np.abs(shape)))
    if abs(shape[pivot]) > 0.0:
        shape = shape * np.exp(-1.0j * np.angle(shape[pivot]))
    return shape


def mode_overlap(left_shape: np.ndarray, right_shape: np.ndarray) -> float:
    denominator = np.linalg.norm(left_shape) * np.linalg.norm(right_shape)
    if denominator == 0.0 or not np.isfinite(denominator):
        return 0.0
    return float(abs(np.vdot(left_shape, right_shape)) / denominator)


def physical_eigenpair_indices(values: Iterable[complex], exclude: list[complex] | None = None) -> list[int]:
    exclude = exclude or []
    indices: list[int] = []
    for index, value in enumerate(values):
        if (
            np.isfinite(value)
            and value.real > 0.05
            and value.imag < -0.02
            and value.real < 2.0
            and value.imag > -3.0
            and all(abs(value - old) > 1.0e-7 for old in exclude)
        ):
            indices.append(index)
    return indices


def select_tracked_mode(
    problem: SpectralProblem,
    values: np.ndarray,
    vectors: np.ndarray,
    target: complex,
    previous_omega: complex | None = None,
    previous_shape: np.ndarray | None = None,
    previous_s_nodes: np.ndarray | None = None,
    exclude: list[complex] | None = None,
    candidate_limit: int = 20,
) -> EigenSelection:
    """Select a physical eigenbranch using frequency, residual, and shape continuity."""

    exclude = exclude or []
    indices = physical_eigenpair_indices(values, exclude=exclude)
    if not indices:
        indices = [index for index, value in enumerate(values) if np.isfinite(value)]
    if not indices:
        raise RuntimeError("No finite eigenvalues found.")

    target_scale = max(abs(target), 0.1)
    if previous_omega is None:
        indices.sort(key=lambda index: abs(values[index] - target) / target_scale)
    else:
        previous_scale = max(abs(previous_omega), 0.1)
        indices.sort(
            key=lambda index: min(
                abs(values[index] - target) / target_scale,
                abs(values[index] - previous_omega) / previous_scale,
            )
        )
    indices = indices[:candidate_limit]

    interpolated_previous: np.ndarray | None = None
    if previous_shape is not None and previous_s_nodes is not None:
        interpolated_previous = interpolate_mode_shape(previous_s_nodes, previous_shape, problem.s_nodes)

    candidates = []
    for index in indices:
        omega = complex(values[index])
        shape = normalized_mode_shape(problem, vectors[:, index])
        overlap = None if interpolated_previous is None else mode_overlap(interpolated_previous, shape)
        freq_term = abs(omega - target) / target_scale
        continuity_term = 0.0
        if previous_omega is not None:
            continuity_term = abs(omega - previous_omega) / max(abs(previous_omega), 0.1)
        residual_at_eigenvalue = residual_norm(problem, omega)
        candidates.append(
            {
                "omega": omega,
                "vector": vectors[:, index],
                "shape": shape,
                "overlap": overlap,
                "freq_term": freq_term,
                "continuity_term": continuity_term,
                "residual_norm": residual_at_eigenvalue,
            }
        )

    residual_scale = min(candidate["residual_norm"] for candidate in candidates) + 1.0e-30
    best = None
    best_score = math.inf
    for candidate in candidates:
        overlap = candidate["overlap"]
        overlap_term = 0.0 if overlap is None else (1.0 - overlap)
        residual_term = min(candidate["residual_norm"] / residual_scale, 100.0) * 1.0e-3
        score = (
            float(candidate["freq_term"])
            + 0.45 * float(candidate["continuity_term"])
            + 0.35 * overlap_term
            + residual_term
        )
        if score < best_score:
            best_score = score
            best = candidate

    if best is None:
        raise RuntimeError("No tracked eigenmode candidate selected.")

    return EigenSelection(
        omega=complex(best["omega"]),
        vector=np.asarray(best["vector"], dtype=complex),
        shape=np.asarray(best["shape"], dtype=complex),
        residual_norm=float(best["residual_norm"]),
        selection_score=float(best_score),
        eigenvector_overlap=None if best["overlap"] is None else float(best["overlap"]),
    )


def select_fundamental_and_overtone(
    values: Iterable[complex],
    fundamental_target: complex,
    overtone_target: complex | None = None,
) -> dict[str, complex]:
    fundamental = select_physical_mode(values, fundamental_target)
    overtone_target = overtone_target or SCHWARZSCHILD_SCALAR_L2_OVERTONE_ESTIMATE
    overtone = select_physical_mode(values, overtone_target, exclude=[fundamental])
    return {"fundamental": fundamental, "first_overtone": overtone}


def residual_operator(problem: SpectralProblem, omega: complex) -> np.ndarray:
    p = spectral_matrix(problem, omega)
    residual = p.conj().T @ p
    return 0.5 * (residual + residual.conj().T)


def residual_norm(problem: SpectralProblem, omega: complex) -> float:
    p = spectral_matrix(problem, omega)
    singular_values = np.linalg.svd(p, compute_uv=False)
    return float(singular_values[-1])


def minimize_residual(problem: SpectralProblem, omega_initial: complex, radius: float = 0.02) -> tuple[complex, float]:
    bounds = [
        (omega_initial.real - radius, omega_initial.real + radius),
        (omega_initial.imag - radius, omega_initial.imag + radius),
    ]

    def objective(params: np.ndarray) -> float:
        omega = complex(params[0], params[1])
        return residual_norm(problem, omega) ** 2

    result = minimize(
        objective,
        np.array([omega_initial.real, omega_initial.imag]),
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": 500, "ftol": 1.0e-18, "gtol": 1.0e-12},
    )
    omega = complex(float(result.x[0]), float(result.x[1]))
    return omega, math.sqrt(float(result.fun))


def residual_diagnostics(problem: SpectralProblem, omega: complex, tol: float = 1.0e-9) -> dict[str, float | int]:
    residual = residual_operator(problem, omega)
    p = spectral_matrix(problem, omega)
    nnz = int(np.count_nonzero(np.abs(residual) > tol))
    total = residual.size
    condition_number = float(np.linalg.cond(p))
    return {
        "matrix_dimension": problem.n,
        "sparsity": 1.0 - nnz / total,
        "condition_number": condition_number,
        "conditioning_warning": condition_number >= CONDITION_WARNING_THRESHOLD,
        "hermiticity_error": float(np.linalg.norm(residual - residual.conj().T)),
        "psd_min_eigenvalue": float(np.linalg.eigvalsh(residual)[0].real),
    }


def run_spectral_study(a_values: list[float], sizes: list[int], baseline_targets: dict[float, complex]) -> list[ModeResult]:
    rows: list[ModeResult] = []
    previous_final_overtone = SCHWARZSCHILD_SCALAR_L2_OVERTONE_ESTIMATE

    for a in a_values:
        fundamental_target = baseline_targets.get(a, SCHWARZSCHILD_SCALAR_L2)
        overtone_target = previous_final_overtone
        previous_by_mode: dict[str, complex] = {}
        previous_shapes_by_mode: dict[str, np.ndarray] = {}
        previous_s_nodes_by_mode: dict[str, np.ndarray] = {}

        for n in sizes:
            problem = build_spectral_problem(a, n)
            values, vectors = generalized_eigenpairs(problem)
            selections: dict[str, EigenSelection] = {}

            fundamental_selection = select_tracked_mode(
                problem,
                values,
                vectors,
                fundamental_target,
                previous_omega=previous_by_mode.get("fundamental"),
                previous_shape=previous_shapes_by_mode.get("fundamental"),
                previous_s_nodes=previous_s_nodes_by_mode.get("fundamental"),
            )
            selections["fundamental"] = fundamental_selection

            overtone_selection = select_tracked_mode(
                problem,
                values,
                vectors,
                overtone_target,
                previous_omega=previous_by_mode.get("first_overtone"),
                previous_shape=previous_shapes_by_mode.get("first_overtone"),
                previous_s_nodes=previous_s_nodes_by_mode.get("first_overtone"),
                exclude=[fundamental_selection.omega],
            )
            selections["first_overtone"] = overtone_selection

            for mode_name, selection in selections.items():
                omega_residual, residual_at_min = minimize_residual(problem, selection.omega)
                diagnostics = residual_diagnostics(problem, omega_residual)
                previous = previous_by_mode.get(mode_name)
                relative_change = None if previous is None else float(abs(selection.omega - previous) / abs(previous))
                if mode_name == "fundamental":
                    branch_status = "publication_high_n" if n == FINAL_SPECTRAL_N else "tracked"
                elif n == OVERTONE_PUBLICATION_N:
                    branch_status = "leaver_validated_reference_grid"
                elif n > OVERTONE_PUBLICATION_N:
                    branch_status = "exploratory_high_n_overtone"
                else:
                    branch_status = "tracked"

                rows.append(
                    ModeResult(
                        a=a,
                        n=n,
                        mode=mode_name,
                        omega=selection.omega,
                        omega_residual=omega_residual,
                        residual_norm=residual_at_min,
                        relative_change=relative_change,
                        hermiticity_error=float(diagnostics["hermiticity_error"]),
                        psd_min_eigenvalue=float(diagnostics["psd_min_eigenvalue"]),
                        matrix_dimension=int(diagnostics["matrix_dimension"]),
                        sparsity=float(diagnostics["sparsity"]),
                        condition_number=float(diagnostics["condition_number"]),
                        conditioning_warning=bool(diagnostics["conditioning_warning"]),
                        selection_score=selection.selection_score,
                        eigenvector_overlap=selection.eigenvector_overlap,
                        branch_status=branch_status,
                    )
                )
                previous_by_mode[mode_name] = selection.omega
                previous_shapes_by_mode[mode_name] = selection.shape
                previous_s_nodes_by_mode[mode_name] = problem.s_nodes

            fundamental_target = selections["fundamental"].omega
            overtone_target = selections["first_overtone"].omega

        previous_final_overtone = overtone_target

    return rows


def run_self_tests() -> list[TestResult]:
    tests: list[TestResult] = []
    problem = build_spectral_problem(0.0, 32)
    values = generalized_eigenvalues(problem)
    omega = select_physical_mode(values, SCHWARZSCHILD_SCALAR_L2)
    residual = residual_operator(problem, omega)
    hermiticity_error = float(np.linalg.norm(residual - residual.conj().T))
    psd_min = float(np.linalg.eigvalsh(residual)[0].real)
    schwarzschild_error = float(abs(omega - SCHWARZSCHILD_SCALAR_L2) / abs(SCHWARZSCHILD_SCALAR_L2))

    tests.append(TestResult("R_N Hermiticity", hermiticity_error, 1.0e-9, hermiticity_error < 1.0e-9))
    tests.append(TestResult("R_N positive semidefinite", max(0.0, -psd_min), 1.0e-8, psd_min >= -1.0e-8))
    tests.append(TestResult("Schwarzschild spectral validation", schwarzschild_error, 1.0e-8, schwarzschild_error < 1.0e-8))

    singular_residual = residual_norm(problem, omega)
    residual_eigenvalue = max(0.0, psd_min)
    consistency_error = abs(singular_residual * singular_residual - residual_eigenvalue)
    tests.append(
        TestResult(
            "Residual singular-value consistency",
            consistency_error,
            1.0e-8,
            consistency_error < 1.0e-8,
        )
    )
    return tests
