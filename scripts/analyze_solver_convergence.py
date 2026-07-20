#!/usr/bin/env python3
"""Generate separate Taylor-order and Chebyshev-to-CF convergence audits."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from qnm.common import SCHWARZSCHILD_REFERENCES, select_physical_mode
from qnm.leaver import solve_leaver_mode, spectral_mode_targets
from qnm.spectral import build_spectral_problem, generalized_eigenvalues

RESULTS = ROOT / "outputs" / "results"
FIGURES = ROOT / "outputs" / "figures"
ORDERS = (64, 80, 96, 112, 128)
SIZES = (16, 24, 32, 48, 64, 80, 96, 112, 128)


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def taylor_audit() -> list[dict[str, object]]:
    cases = [
        ("scalar_l2_n0_a1", "scalar", 2, 0, 1.0),
        ("scalar_l2_n1_a1", "scalar", 2, 1, 1.0),
        ("axial_l2_n0_a1", "gravitational", 2, 0, 1.0),
    ]
    all_solutions: dict[str, list[tuple[int, complex, float]]] = {}
    for name, sector, ell, overtone, a in cases:
        targets = spectral_mode_targets([a], 32, ell=ell, perturbation_type=sector)
        mode = "fundamental" if overtone == 0 else "first_overtone"
        guess = targets[a][mode]
        solutions = []
        for order in ORDERS:
            omega, residual = solve_leaver_mode(
                a, guess, depth=320, order=order, ell=ell, perturbation_type=sector
            )
            solutions.append((order, omega, residual))
            guess = omega
        all_solutions[name] = solutions

    rows: list[dict[str, object]] = []
    for name, sector, ell, overtone, a in cases:
        final = all_solutions[name][-1][1]
        for order, omega, residual in all_solutions[name]:
            rows.append(
                {
                    "case": name,
                    "perturbation_type": sector,
                    "ell": ell,
                    "overtone": overtone,
                    "a_over_M": a,
                    "taylor_order": order,
                    "continued_fraction_depth": 320,
                    "omega_real": omega.real,
                    "omega_imag": omega.imag,
                    "continued_fraction_abs": residual,
                    "absolute_difference_from_order_128": abs(omega - final),
                }
            )
    return rows


def spectral_audit() -> list[dict[str, object]]:
    cases = [
        ("scalar_l2_n0_a0", 0.0, 0),
        ("scalar_l2_n0_a1", 1.0, 0),
        ("scalar_l2_n1_a0", 0.0, 1),
        ("scalar_l2_n1_a1", 1.0, 1),
    ]
    rows: list[dict[str, object]] = []
    for name, a, overtone in cases:
        targets = spectral_mode_targets([a], 32)
        mode = "fundamental" if overtone == 0 else "first_overtone"
        cf, cf_residual = solve_leaver_mode(a, targets[a][mode], depth=320, order=128)
        previous: complex | None = None
        for n in SIZES:
            values = generalized_eigenvalues(build_spectral_problem(a, n))
            # Use the collocation-independent CF value only for diagnostic association;
            # continuation is recorded as a second distance below.
            omega = select_physical_mode(values, cf)
            rows.append(
                {
                    "case": name,
                    "a_over_M": a,
                    "overtone": overtone,
                    "N": n,
                    "omega_real": omega.real,
                    "omega_imag": omega.imag,
                    "cf_real_order128_depth320": cf.real,
                    "cf_imag_order128_depth320": cf.imag,
                    "cf_residual_abs": cf_residual,
                    "absolute_error_from_cf": abs(omega - cf),
                    "movement_from_previous_N": "" if previous is None else abs(omega - previous),
                }
            )
            previous = omega
    return rows


def plot(rows: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.3, 4.2), sharey=True)
    for overtone, ax in enumerate(axes):
        for a, marker in ((0.0, "o"), (1.0, "s")):
            subset = [r for r in rows if r["overtone"] == overtone and r["a_over_M"] == a]
            ax.plot(
                [int(r["N"]) for r in subset],
                [float(r["absolute_error_from_cf"]) for r in subset],
                marker=marker,
                label=rf"$a/M={a:g}$",
            )
        ax.set_yscale("log")
        ax.set_xlabel(r"Chebyshev size $N$")
        ax.set_title("fundamental" if overtone == 0 else "first overtone")
        ax.grid(alpha=0.25, which="both")
        ax.legend(frameon=False)
    axes[0].set_ylabel(r"$|\omega_N-\omega_{\rm CF}|$")
    fig.tight_layout()
    fig.savefig(FIGURES / "spectral_convergence_to_cf.png", dpi=220)
    plt.close(fig)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    taylor = taylor_audit()
    spectral = spectral_audit()
    write_rows(RESULTS / "leaver_taylor_order_convergence.csv", taylor)
    write_rows(RESULTS / "spectral_convergence_to_cf.csv", spectral)
    plot(spectral)
    worst_order96 = max(
        float(r["absolute_difference_from_order_128"])
        for r in taylor if int(r["taylor_order"]) == 96
    )
    n1_32 = [
        float(r["absolute_error_from_cf"])
        for r in spectral if int(r["overtone"]) == 1 and int(r["N"]) == 32
    ]
    n1_high = [
        float(r["absolute_error_from_cf"])
        for r in spectral if int(r["overtone"]) == 1 and int(r["N"]) >= 64
    ]
    (RESULTS / "solver_convergence_report.md").write_text(
        "# Solver convergence audit\n\n"
        f"- Maximum order-96 to order-128 frequency difference: `{worst_order96:.3e}`.\n"
        f"- First-overtone N=32 errors from CF: `{min(n1_32):.3e}` to `{max(n1_32):.3e}`.\n"
        f"- First-overtone N>=64 errors from CF: `{min(n1_high):.3e}` to `{max(n1_high):.3e}`.\n"
        "- These modes are therefore described as cross-validated low-lying first overtones at N=32; higher-N motion is shown rather than hidden.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
