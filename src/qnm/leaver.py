"""Leaver-style Frobenius continued-fraction validation for KS scalar QNMs.

This module is intentionally separate from the Chebyshev spectral residual
implementation. It constructs a horizon Frobenius series for the same factored
scalar perturbation equation used by the spectral solver, reduces the resulting
multi-term recurrence to a three-term recurrence, and solves the associated
continued-fraction condition.

The scalar wave equation is

    d^2 Psi / dr_*^2 + [omega^2 - V_l(r; a)] Psi = 0,

with the QNM behavior factored as

    Psi = exp(i omega r) r^(2 i M omega) s^(-i omega / f'(r_h)) u(s),
    s = 1 - r_h/r.

After substitution, the regular part obeys

    B_2(s) u'' + B_1(s) u' + B_0(s) u = 0,

where B_j are analytic at the horizon. The Frobenius expansion
u(s)=sum_n a_n s^n gives a finite Taylor-truncated multi-term recurrence. We
use Gaussian elimination to reduce it to

    alpha_n a_{n+1} + beta_n a_n + gamma_n a_{n-1} = 0,

and impose Leaver's minimal-solution condition through the continued fraction

    beta_0 - alpha_0 gamma_1 /
        (beta_1 - alpha_1 gamma_2 / (beta_2 - ...)) = 0.

This is an independent validation layer in the numerical sense: it avoids
Chebyshev collocation, waveform data, and residual minimization, while still
sharing the same perturbation equation, compact coordinate, endpoint
factorization, and potential model as the spectral solver.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import root

from .common import (
    A_VALUES,
    ELL,
    MASS,
    SCHWARZSCHILD_SCALAR_L2,
    SCHWARZSCHILD_SCALAR_L2_OVERTONE_ESTIMATE,
    df_ks,
    horizon_radius,
    select_physical_mode,
)
from .spectral import build_spectral_problem, generalized_eigenvalues


DEFAULT_SPECTRAL_VALIDATION_N = 32
DEFAULT_TAYLOR_ORDER = 96
DEFAULT_CF_DEPTH = 240
DEFAULT_VALIDATION_THRESHOLD = 1.0e-6


@dataclass
class LeaverModeResult:
    a: float
    mode: str
    omega_leaver: complex
    omega_spectral: complex
    spectral_n: int
    relative_difference: float
    schwarzschild_reference_relative_error: float | None
    continued_fraction_abs: float
    taylor_order: int
    cf_depth: int
    validation_threshold: float
    validation_passed: bool


Series = dict[int, complex]


def _trim(series: Series, max_order: int, min_order: int = -8) -> Series:
    return {
        power: coeff
        for power, coeff in series.items()
        if min_order <= power <= max_order and abs(coeff) > 0.0
    }


def _add(left: Series, right: Series, max_order: int) -> Series:
    result = dict(left)
    for power, coeff in right.items():
        result[power] = result.get(power, 0.0j) + coeff
    return _trim(result, max_order)


def _scale(series: Series, factor: complex, max_order: int) -> Series:
    return _trim({power: factor * coeff for power, coeff in series.items()}, max_order)


def _mul(left: Series, right: Series, max_order: int) -> Series:
    result: Series = {}
    for left_power, left_coeff in left.items():
        for right_power, right_coeff in right.items():
            power = left_power + right_power
            if -8 <= power <= max_order:
                result[power] = result.get(power, 0.0j) + left_coeff * right_coeff
    return _trim(result, max_order)


def _shift(series: Series, offset: int, max_order: int) -> Series:
    return _trim({power + offset: coeff for power, coeff in series.items()}, max_order)


def _sqrt_series(series: Series, max_order: int) -> Series:
    coeffs = [0.0j] * (max_order + 1)
    coeffs[0] = complex(math.sqrt(series.get(0, 0.0).real))
    for n in range(1, max_order + 1):
        convolution = sum(coeffs[k] * coeffs[n - k] for k in range(1, n))
        coeffs[n] = (series.get(n, 0.0j) - convolution) / (2.0 * coeffs[0])
    return {n: coeff for n, coeff in enumerate(coeffs) if abs(coeff) > 0.0}


def _inverse_series(series: Series, max_order: int) -> Series:
    coeffs = [0.0j] * (max_order + 1)
    coeffs[0] = 1.0 / series.get(0, 0.0j)
    for n in range(1, max_order + 1):
        convolution = sum(series.get(k, 0.0j) * coeffs[n - k] for k in range(1, n + 1))
        coeffs[n] = -convolution / series.get(0, 0.0j)
    return {n: coeff for n, coeff in enumerate(coeffs) if abs(coeff) > 0.0}


def frobenius_coefficient_series(
    a: float,
    omega: complex,
    order: int = DEFAULT_TAYLOR_ORDER,
    ell: int = ELL,
    mass: float = MASS,
    perturbation_type: str = "scalar",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return Taylor coefficients of B_2, B_1, B_0 in the factored equation.

    The KS square root is expanded formally through

        sqrt(r^2-a^2) = sqrt(r_h^2 - a^2(1-s)^2)/(1-s),

    so the coefficient construction remains independent of the Chebyshev
    collocation machinery.
    """

    rh = horizon_radius(a, mass)
    fh = float(df_ks(rh, a, mass))

    one: Series = {0: 1.0 + 0.0j}
    s_power: Series = {1: 1.0 + 0.0j}
    one_minus_s: Series = {0: 1.0 + 0.0j, 1: -1.0 + 0.0j}
    one_minus_s2 = _mul(one_minus_s, one_minus_s, order)
    one_minus_s3 = _mul(one_minus_s2, one_minus_s, order)

    sqrt_argument: Series = {
        0: (2.0 * mass) ** 2 + 0.0j,
        1: 2.0 * a * a + 0.0j,
        2: -a * a + 0.0j,
    }
    q_series = _sqrt_series(sqrt_argument, order)
    inv_q_series = _inverse_series(q_series, order)

    f_series = _scale(
        _add(q_series, _scale(one_minus_s, -2.0 * mass, order), order),
        1.0 / rh,
        order,
    )
    fp_series = _add(
        _scale(_mul(one_minus_s3, inv_q_series, order), a * a / (rh * rh), order),
        _scale(one_minus_s2, 2.0 * mass / (rh * rh), order),
        order,
    )

    ds_dr = _scale(one_minus_s2, 1.0 / rh, order)
    d2s_dr2 = _scale(one_minus_s3, -2.0 / (rh * rh), order)
    inv_r = _scale(one_minus_s, 1.0 / rh, order)
    inv_r2 = _scale(one_minus_s2, 1.0 / (rh * rh), order)

    if perturbation_type == "scalar":
        potential_inner = _add(
            _scale(inv_r2, ell * (ell + 1.0), order),
            _mul(fp_series, inv_r, order),
            order,
        )
    elif perturbation_type == "gravitational":
        inv_r3 = _mul(inv_r2, inv_r, order)
        potential_inner = _add(
            _scale(inv_r2, ell * (ell + 1.0), order),
            _scale(inv_r3, -6.0 * mass, order),
            order,
        )
    else:
        raise ValueError(f"Unknown perturbation type: {perturbation_type}")

    potential = _mul(f_series, potential_inner, order)

    dlog_s_dr = _shift(ds_dr, -1, order)
    d_dlog_s_dr = _shift(
        _scale(_mul(_add(one, s_power, order), one_minus_s3, order), -1.0 / (rh * rh), order),
        -2,
        order,
    )

    h1 = _scale(
        _add(
            _add(one, _scale(inv_r, 2.0 * mass, order), order),
            _scale(dlog_s_dr, -1.0 / fh, order),
            order,
        ),
        1.0j,
        order,
    )
    dh1_dr = _scale(
        _add(
            _scale(inv_r2, -2.0 * mass, order),
            _scale(d_dlog_s_dr, -1.0 / fh, order),
            order,
        ),
        1.0j,
        order,
    )

    f2 = _mul(f_series, f_series, order)
    b2 = _shift(_mul(_mul(f2, ds_dr, order), ds_dr, order), -1, order)
    b1 = _shift(
        _add(
            _add(_mul(f2, d2s_dr2, order), _mul(_mul(f_series, fp_series, order), ds_dr, order), order),
            _scale(_mul(_mul(f2, h1, order), ds_dr, order), 2.0 * omega, order),
            order,
        ),
        -1,
        order,
    )
    b0 = _shift(
        _add(
            _add(
                _scale(potential, -1.0, order),
                _scale(_add(_mul(f2, dh1_dr, order), _mul(_mul(f_series, fp_series, order), h1, order), order), omega, order),
                order,
            ),
            _scale(_add(one, _mul(_mul(f2, h1, order), h1, order), order), omega * omega, order),
            order,
        ),
        -1,
        order,
    )

    return (
        np.array([b2.get(n, 0.0j) for n in range(order + 1)], dtype=complex),
        np.array([b1.get(n, 0.0j) for n in range(order + 1)], dtype=complex),
        np.array([b0.get(n, 0.0j) for n in range(order + 1)], dtype=complex),
    )


def recurrence_rows(
    a: float,
    omega: complex,
    depth: int = DEFAULT_CF_DEPTH,
    order: int = DEFAULT_TAYLOR_ORDER,
    ell: int = ELL,
    perturbation_type: str = "scalar",
) -> np.ndarray:
    """Build the Taylor-truncated Frobenius recurrence rows.

    Row n stores coefficients of a_{n+1}, a_n, a_{n-1}, ... by offset.
    """

    b2, b1, b0 = frobenius_coefficient_series(
        a,
        omega,
        order=order,
        ell=ell,
        perturbation_type=perturbation_type,
    )
    rows = np.zeros((depth + 1, order + 3), dtype=complex)

    for n in range(depth + 1):
        for k, coeff in enumerate(b0):
            index = n - k
            offset = n + 1 - index
            if index >= 0 and offset < rows.shape[1]:
                rows[n, offset] += coeff
        for k, coeff in enumerate(b1):
            index = n - k + 1
            offset = n + 1 - index
            if index >= 0 and offset < rows.shape[1]:
                rows[n, offset] += coeff * index
        for k, coeff in enumerate(b2):
            index = n - k + 2
            offset = n + 1 - index
            if index >= 0 and offset < rows.shape[1]:
                rows[n, offset] += coeff * index * (index - 1)

    while rows.shape[1] > 3 and np.all(np.abs(rows[:, -1]) < 1.0e-13):
        rows = rows[:, :-1]
    return rows


def reduce_recurrence_to_three_terms(rows: np.ndarray) -> np.ndarray:
    """Gaussian-eliminate a multi-term recurrence to three terms."""

    reduced = rows.astype(complex, copy=True)
    depth = reduced.shape[0] - 1
    max_offset = reduced.shape[1] - 1

    for offset in range(max_offset, 2, -1):
        for n in range(offset - 1, depth + 1):
            denominator = reduced[n - 1, offset - 1]
            if abs(denominator) < 1.0e-24:
                continue
            factor = reduced[n, offset] / denominator
            reduced[n, 1 : offset + 1] -= factor * reduced[n - 1, 0:offset]
            reduced[n, offset] = 0.0

    return reduced[:, :3]


def continued_fraction_residual(
    a: float,
    omega: complex,
    depth: int = DEFAULT_CF_DEPTH,
    order: int = DEFAULT_TAYLOR_ORDER,
    ell: int = ELL,
    perturbation_type: str = "scalar",
) -> complex:
    """Evaluate the finite-depth Leaver continued-fraction equation."""

    three_term = reduce_recurrence_to_three_terms(
        recurrence_rows(
            a,
            omega,
            depth=depth,
            order=order,
            ell=ell,
            perturbation_type=perturbation_type,
        )
    )
    alpha = three_term[:, 0]
    beta = three_term[:, 1]
    gamma = three_term[:, 2]

    tail = beta[depth]
    for n in range(depth - 1, 0, -1):
        tail = beta[n] - alpha[n] * gamma[n + 1] / tail
    return beta[0] - alpha[0] * gamma[1] / tail


def solve_leaver_mode(
    a: float,
    initial_guess: complex,
    depth: int = DEFAULT_CF_DEPTH,
    order: int = DEFAULT_TAYLOR_ORDER,
    ell: int = ELL,
    perturbation_type: str = "scalar",
    residual_tolerance: float = 1.0e-7,
) -> tuple[complex, float]:
    """Solve the continued-fraction root near an initial QNM guess."""

    def residual_components(params: np.ndarray) -> list[float]:
        omega = complex(float(params[0]), float(params[1]))
        residual = continued_fraction_residual(
            a,
            omega,
            depth=depth,
            order=order,
            ell=ell,
            perturbation_type=perturbation_type,
        )
        return [float(residual.real), float(residual.imag)]

    solution = root(
        residual_components,
        np.array([initial_guess.real, initial_guess.imag], dtype=float),
        method="hybr",
        tol=1.0e-11,
    )
    omega = complex(float(solution.x[0]), float(solution.x[1]))
    residual_abs = abs(
        continued_fraction_residual(
            a,
            omega,
            depth=depth,
            order=order,
            ell=ell,
            perturbation_type=perturbation_type,
        )
    )
    if residual_abs > residual_tolerance:
        raise RuntimeError(
            f"Leaver root solve failed for a={a:g}, guess={initial_guess}: "
            f"{solution.message}; |CF|={residual_abs:.3e}"
        )
    return omega, float(residual_abs)


def spectral_mode_targets(
    a_values: list[float],
    spectral_n: int,
    ell: int = ELL,
    perturbation_type: str = "scalar",
) -> dict[float, dict[str, complex]]:
    """Compute spectral target modes without changing the spectral module."""

    targets: dict[float, dict[str, complex]] = {}
    fundamental_target = SCHWARZSCHILD_SCALAR_L2
    overtone_target = SCHWARZSCHILD_SCALAR_L2_OVERTONE_ESTIMATE

    for a in a_values:
        problem = build_spectral_problem(a, spectral_n, ell=ell, perturbation_type=perturbation_type)
        values = generalized_eigenvalues(problem)
        fundamental = select_physical_mode(values, fundamental_target)
        overtone = select_physical_mode(values, overtone_target, exclude=[fundamental])
        targets[a] = {"fundamental": fundamental, "first_overtone": overtone}
        fundamental_target = fundamental
        overtone_target = overtone

    return targets


def run_leaver_validation(
    a_values: list[float] | None = None,
    spectral_n: int = DEFAULT_SPECTRAL_VALIDATION_N,
    depth: int = DEFAULT_CF_DEPTH,
    order: int = DEFAULT_TAYLOR_ORDER,
    threshold: float = DEFAULT_VALIDATION_THRESHOLD,
    ell: int = ELL,
    perturbation_type: str = "scalar",
) -> list[LeaverModeResult]:
    """Solve Leaver modes and compare them against the spectral implementation."""

    a_values = A_VALUES if a_values is None else a_values
    targets = spectral_mode_targets(a_values, spectral_n, ell=ell, perturbation_type=perturbation_type)
    rows: list[LeaverModeResult] = []

    for a in a_values:
        for mode_name in ("fundamental", "first_overtone"):
            spectral_omega = targets[a][mode_name]
            omega_leaver, cf_abs = solve_leaver_mode(
                a,
                spectral_omega,
                depth=depth,
                order=order,
                ell=ell,
                perturbation_type=perturbation_type,
            )
            relative_difference = float(abs(omega_leaver - spectral_omega) / abs(spectral_omega))
            reference_error: float | None = None
            if a == 0.0 and mode_name == "fundamental":
                reference_error = float(abs(omega_leaver - SCHWARZSCHILD_SCALAR_L2) / abs(SCHWARZSCHILD_SCALAR_L2))

            rows.append(
                LeaverModeResult(
                    a=a,
                    mode=mode_name,
                    omega_leaver=omega_leaver,
                    omega_spectral=spectral_omega,
                    spectral_n=spectral_n,
                    relative_difference=relative_difference,
                    schwarzschild_reference_relative_error=reference_error,
                    continued_fraction_abs=cf_abs,
                    taylor_order=order,
                    cf_depth=depth,
                    validation_threshold=threshold,
                    validation_passed=relative_difference < threshold,
                )
            )

    return rows


def write_leaver_comparison(output: Path, rows: list[LeaverModeResult]) -> None:
    output.parent.mkdir(exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "a_over_M",
                "mode",
                "spectral_N",
                "spectral_real",
                "spectral_imag",
                "leaver_real",
                "leaver_imag",
                "relative_difference_spectral_leaver",
                "schwarzschild_reference_relative_error",
                "continued_fraction_abs",
                "taylor_order",
                "continued_fraction_depth",
                "validation_threshold",
                "validation_passed",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.a,
                    row.mode,
                    row.spectral_n,
                    row.omega_spectral.real,
                    row.omega_spectral.imag,
                    row.omega_leaver.real,
                    row.omega_leaver.imag,
                    row.relative_difference,
                    "" if row.schwarzschild_reference_relative_error is None else row.schwarzschild_reference_relative_error,
                    row.continued_fraction_abs,
                    row.taylor_order,
                    row.cf_depth,
                    row.validation_threshold,
                    row.validation_passed,
                ]
            )


def assert_leaver_validation(rows: list[LeaverModeResult]) -> None:
    failed = [row for row in rows if not row.validation_passed]
    if failed:
        details = ", ".join(
            f"a/M={row.a:g} {row.mode}: {row.relative_difference:.3e}" for row in failed
        )
        raise AssertionError(f"Leaver/spectral disagreement exceeded tolerance: {details}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("outputs/results/leaver_comparison.csv"))
    parser.add_argument("--spectral-n", type=int, default=DEFAULT_SPECTRAL_VALIDATION_N)
    parser.add_argument("--order", type=int, default=DEFAULT_TAYLOR_ORDER)
    parser.add_argument("--depth", type=int, default=DEFAULT_CF_DEPTH)
    parser.add_argument("--threshold", type=float, default=DEFAULT_VALIDATION_THRESHOLD)
    args = parser.parse_args()

    rows = run_leaver_validation(
        spectral_n=args.spectral_n,
        depth=args.depth,
        order=args.order,
        threshold=args.threshold,
    )
    write_leaver_comparison(args.output, rows)
    assert_leaver_validation(rows)
    for row in rows:
        print(
            f"a/M={row.a:g}, {row.mode}: "
            f"Leaver={row.omega_leaver.real:.12f}{row.omega_leaver.imag:+.12f}i, "
            f"spectral(N={row.spectral_n})={row.omega_spectral.real:.12f}{row.omega_spectral.imag:+.12f}i, "
            f"rel_diff={row.relative_difference:.3e}"
        )


if __name__ == "__main__":
    main()
