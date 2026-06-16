"""Baseline time-domain, ringdown-fit, and matrix-pencil QNM workflow.

This module intentionally preserves the earlier proof-of-concept pipeline:

    time-domain evolution -> waveform fitting -> matrix-pencil diagnostic.

The matrix-pencil operator is not treated as the final spectral operator.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.linalg import eig, logm, svd
from scipy.optimize import curve_fit, minimize_scalar

from .common import ELL, MASS, WINDOW, f_ks, horizon_radius, scalar_potential, select_physical_mode


@dataclass
class Background:
    a: float
    x: np.ndarray
    r: np.ndarray
    potential: np.ndarray
    horizon: float
    r_peak: float
    x_table_min: float
    x_table_max: float


@dataclass
class BaselineResult:
    a: float
    omega_fit: complex
    omega_pencil: complex
    fit_rms_relative: float
    r_horizon: float
    r_peak: float


def build_background(
    a: float,
    ell: int,
    x: np.ndarray,
    mass: float = MASS,
    r_table_max: float = 1200.0,
    n_table: int = 220_000,
) -> Background:
    horizon = horizon_radius(a, mass)
    eps = 1.0e-10
    offsets = np.geomspace(eps, r_table_max - horizon, n_table)
    r_table = horizon + offsets
    rstar_table = cumulative_trapezoid(1.0 / f_ks(r_table, a, mass), r_table, initial=0.0)

    peak = minimize_scalar(
        lambda rr: -float(scalar_potential(np.array([rr]), ell, a, mass)[0]),
        bounds=(horizon + 1.0e-7, 20.0 * mass),
        method="bounded",
    )
    r_peak = float(peak.x)
    rstar_peak = float(np.interp(r_peak, r_table, rstar_table))
    rstar_table = rstar_table - rstar_peak

    r_of_x = PchipInterpolator(rstar_table, r_table, extrapolate=False)
    mask = (x >= rstar_table[0]) & (x <= rstar_table[-1])
    r = np.full_like(x, np.nan, dtype=float)
    potential = np.zeros_like(x, dtype=float)
    r[mask] = r_of_x(x[mask])
    potential[mask] = scalar_potential(r[mask], ell, a, mass)

    return Background(
        a=a,
        x=x,
        r=r,
        potential=potential,
        horizon=horizon,
        r_peak=r_peak,
        x_table_min=float(rstar_table[0]),
        x_table_max=float(rstar_table[-1]),
    )


def evolve_waveform(
    background: Background,
    dx: float,
    dt: float,
    tmax: float,
    obs_x: float = 20.0,
    source_x: float = 0.0,
    sigma: float = 4.0,
) -> tuple[np.ndarray, np.ndarray]:
    x = background.x
    potential = background.potential
    psi0 = np.exp(-((x - source_x) / sigma) ** 2)
    cfl2 = (dt / dx) ** 2
    lap = np.zeros_like(psi0)
    lap[1:-1] = psi0[2:] - 2.0 * psi0[1:-1] + psi0[:-2]

    previous = psi0.copy()
    current = psi0 + 0.5 * cfl2 * lap - 0.5 * dt * dt * potential * psi0
    obs_index = int(np.argmin(np.abs(x - obs_x)))
    boundary_q = (dt - dx) / (dt + dx)

    times: list[float] = []
    signal: list[float] = []
    for n in range(1, int(tmax / dt)):
        nxt = np.zeros_like(current)
        nxt[1:-1] = (
            -previous[1:-1]
            + 2.0 * current[1:-1]
            + cfl2 * (current[2:] - 2.0 * current[1:-1] + current[:-2])
            - dt * dt * potential[1:-1] * current[1:-1]
        )
        nxt[0] = current[1] + boundary_q * (nxt[1] - current[0])
        nxt[-1] = current[-2] + boundary_q * (nxt[-2] - current[-1])
        previous, current = current, nxt

        if n % 2 == 0:
            times.append((n + 1) * dt)
            signal.append(float(current[obs_index]))

    return np.array(times), np.array(signal)


def fit_ringdown(
    times: np.ndarray,
    signal: np.ndarray,
    window: tuple[float, float] = WINDOW,
) -> tuple[complex, float, np.ndarray, np.ndarray, np.ndarray]:
    mask = (times >= window[0]) & (times <= window[1])
    t = times[mask]
    y = signal[mask]
    t0 = window[0]

    def model(tvals: np.ndarray, amp: float, alpha: float, omega: float, phi: float, offset: float) -> np.ndarray:
        return amp * np.exp(-alpha * (tvals - t0)) * np.cos(omega * tvals + phi) + offset

    p0 = [0.5 * float(np.ptp(y)), 0.095, 0.48, 0.0, 0.0]
    bounds = ([-np.inf, 0.01, 0.1, -20.0, -np.inf], [np.inf, 0.5, 1.5, 20.0, np.inf])
    params, _ = curve_fit(model, t, y, p0=p0, bounds=bounds, maxfev=50_000)
    fitted = model(t, *params)
    rms_relative = float(np.sqrt(np.mean((fitted - y) ** 2)) / np.sqrt(np.mean(y**2)))
    _, alpha, omega, _, _ = params
    return complex(omega, -alpha), rms_relative, t, y, fitted


def matrix_pencil_operator(
    times: np.ndarray,
    signal: np.ndarray,
    rank: int = 4,
    window: tuple[float, float] = WINDOW,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the old waveform-derived diagnostic operator.

    The current spectral method does not use this matrix as the final operator.
    It is preserved only for comparison with the earlier proof-of-concept run.
    """

    mask = (times >= window[0]) & (times <= window[1])
    y = signal[mask].astype(complex)
    y = y - np.mean(y)
    dt = float(times[1] - times[0])
    m = len(y)
    rows = m // 2
    cols = m - rows - 1

    h0 = np.column_stack([y[i : i + rows] for i in range(cols)])
    h1 = np.column_stack([y[i + 1 : i + rows + 1] for i in range(cols)])
    u, s, vh = svd(h0, full_matrices=False)
    rank = min(rank, len(s))
    reduced_shift = u[:, :rank].conj().T @ h1 @ vh[:rank, :].conj().T @ np.diag(1.0 / s[:rank])
    c_operator = 1j * logm(reduced_shift) / dt
    omega_values = eig(c_operator, left=False, right=False)
    return c_operator, omega_values


def plot_waveform(
    output: Path,
    a: float,
    times: np.ndarray,
    signal: np.ndarray,
    fit_t: np.ndarray,
    fit_y: np.ndarray,
) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    ax.plot(times, signal, color="#1f5f8b", lw=1.0, label="time-domain waveform")
    ax.plot(fit_t, fit_y, color="#d1495b", lw=2.0, label="ringdown fit")
    ax.axvspan(WINDOW[0], WINDOW[1], color="#f2cc8f", alpha=0.22, label="fit window")
    ax.set_xlim(0, 220)
    ax.set_xlabel(r"$t/M$")
    ax.set_ylabel(r"$\Psi(t,r_*^{obs})$")
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right", frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def run_baseline(a_values: list[float], figures_dir: Path) -> dict[float, BaselineResult]:
    dx = 0.1
    dt = 0.05
    tmax = 320.0
    x = np.arange(-250.0, 500.0 + dx / 2.0, dx)
    results: dict[float, BaselineResult] = {}

    for a in a_values:
        background = build_background(a, ELL, x)
        times, signal = evolve_waveform(background, dx=dx, dt=dt, tmax=tmax)
        omega_fit, fit_rms, fit_t, _, fit_y = fit_ringdown(times, signal)
        _, omega_values = matrix_pencil_operator(times, signal, rank=4)
        omega_pencil = select_physical_mode(omega_values, omega_fit)
        plot_waveform(figures_dir / f"waveform_fit_a_{a:g}.png", a, times, signal, fit_t, fit_y)
        results[a] = BaselineResult(
            a=a,
            omega_fit=omega_fit,
            omega_pencil=omega_pencil,
            fit_rms_relative=fit_rms,
            r_horizon=background.horizon,
            r_peak=background.r_peak,
        )
    return results
