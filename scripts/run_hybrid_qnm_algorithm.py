#!/usr/bin/env python3
"""Run baseline and direct spectral-residual KS scalar QNM pipelines.

The matrix-pencil workflow is preserved as a baseline diagnostic. The upgraded
operator path is:

    KS perturbation equation
    -> Chebyshev pseudospectral discretization
    -> P_N(omega) u = 0
    -> R_N(omega)=P_N(omega)^dagger P_N(omega)
    -> residual minimization
    -> convergence and Leaver-style validation
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import scipy

from qnm.analysis import write_physics_analysis
from qnm.baseline import BaselineResult, run_baseline
from qnm.catalogue import (
    assert_catalogue_validation,
    plot_mode_trajectories,
    run_catalogue,
    write_catalogue,
    write_catalogue_report,
)
from qnm.common import A_VALUES, FINAL_SPECTRAL_N, SCHWARZSCHILD_SCALAR_L2, SPECTRAL_SIZES
from qnm.leaver import assert_leaver_validation, run_leaver_validation, write_leaver_comparison
from qnm.spectral import (
    CONDITION_WARNING_THRESHOLD,
    ModeResult,
    OVERTONE_PUBLICATION_N,
    run_self_tests,
    run_spectral_study,
)


def final_fundamental_rows(rows: list[ModeResult]) -> list[ModeResult]:
    return [row for row in rows if row.n == FINAL_SPECTRAL_N and row.mode == "fundamental"]


def publication_mode_rows(rows: list[ModeResult]) -> list[ModeResult]:
    """Rows safe for headline tables.

    Fundamentals use the high-resolution tracked branch. The first overtone is
    frozen at the Leaver-validated catalogue grid until high-N overtone branch
    tracking has stronger independent validation.
    """

    selected = [
        row
        for row in rows
        if (row.mode == "fundamental" and row.n == FINAL_SPECTRAL_N)
        or (row.mode == "first_overtone" and row.n == OVERTONE_PUBLICATION_N)
    ]
    return sorted(selected, key=lambda row: (row.a, row.mode != "fundamental", row.n))


def write_spectral_results(
    output: Path,
    baselines: dict[float, BaselineResult],
    rows: list[ModeResult],
    publication_only: bool = True,
) -> None:
    rows_to_write = publication_mode_rows(rows) if publication_only else rows
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "a_over_M",
                "N",
                "mode",
                "baseline_time_domain_real",
                "baseline_time_domain_imag",
                "baseline_matrix_pencil_real",
                "baseline_matrix_pencil_imag",
                "spectral_real",
                "spectral_imag",
                "spectral_residual_real",
                "spectral_residual_imag",
                "schwarzschild_reference_relative_error",
                "matrix_dimension",
                "effective_qubits",
                "pauli_terms",
                "sparsity",
                "condition_number",
                "conditioning_warning",
                "residual_norm",
                "hermiticity_error",
                "psd_min_eigenvalue",
                "branch_status",
                "selection_score",
                "eigenvector_overlap",
            ]
        )
        for row in rows_to_write:
            baseline = baselines[row.a]
            if row.a == 0.0 and row.mode == "fundamental":
                ref_error = abs(row.omega - SCHWARZSCHILD_SCALAR_L2) / abs(SCHWARZSCHILD_SCALAR_L2)
            else:
                ref_error = ""

            writer.writerow(
                [
                    row.a,
                    row.n,
                    row.mode,
                    baseline.omega_fit.real if row.mode == "fundamental" else "",
                    baseline.omega_fit.imag if row.mode == "fundamental" else "",
                    baseline.omega_pencil.real if row.mode == "fundamental" else "",
                    baseline.omega_pencil.imag if row.mode == "fundamental" else "",
                    row.omega.real,
                    row.omega.imag,
                    row.omega_residual.real,
                    row.omega_residual.imag,
                    ref_error,
                    row.matrix_dimension,
                    row.effective_qubits,
                    row.pauli_terms,
                    row.sparsity,
                    row.condition_number,
                    row.conditioning_warning,
                    row.residual_norm,
                    row.hermiticity_error,
                    row.psd_min_eigenvalue,
                    row.branch_status,
                    "" if row.selection_score is None else row.selection_score,
                    "" if row.eigenvector_overlap is None else row.eigenvector_overlap,
                ]
            )


def write_convergence_table(output: Path, rows: list[ModeResult]) -> None:
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "a_over_M",
                "N",
                "mode",
                "omega_real",
                "omega_imag",
                "relative_change_from_previous_N",
                "residual_norm",
                "matrix_dimension",
                "effective_qubits",
                "pauli_terms",
                "sparsity",
                "condition_number",
                "conditioning_warning",
                "branch_status",
                "selection_score",
                "eigenvector_overlap",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.a,
                    row.n,
                    row.mode,
                    row.omega.real,
                    row.omega.imag,
                    "" if row.relative_change is None else row.relative_change,
                    row.residual_norm,
                    row.matrix_dimension,
                    row.effective_qubits,
                    row.pauli_terms,
                    row.sparsity,
                    row.condition_number,
                    row.conditioning_warning,
                    row.branch_status,
                    "" if row.selection_score is None else row.selection_score,
                    "" if row.eigenvector_overlap is None else row.eigenvector_overlap,
                ]
            )


def write_resource_estimates(output: Path, rows: list[ModeResult]) -> None:
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "a_over_M",
                "N",
                "mode",
                "matrix_dimension",
                "effective_qubits",
                "pauli_terms",
                "sparsity",
                "condition_number",
                "conditioning_warning",
                "residual_norm",
                "hermiticity_error",
                "psd_min_eigenvalue",
                "branch_status",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.a,
                    row.n,
                    row.mode,
                    row.matrix_dimension,
                    row.effective_qubits,
                    row.pauli_terms,
                    row.sparsity,
                    row.condition_number,
                    row.conditioning_warning,
                    row.residual_norm,
                    row.hermiticity_error,
                    row.psd_min_eigenvalue,
                    row.branch_status,
                ]
            )


def write_run_metadata(output: Path) -> None:
    metadata = {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
        "spectral_sizes": SPECTRAL_SIZES,
        "final_spectral_n_for_fundamentals": FINAL_SPECTRAL_N,
        "publication_overtone_n": OVERTONE_PUBLICATION_N,
        "condition_warning_threshold": CONDITION_WARNING_THRESHOLD,
        "note": (
            "Publication-facing spectral_results.csv freezes first overtones at "
            "the Leaver-validated reference grid. exploratory_spectral_results.csv "
            "contains tracked high-N overtone rows for diagnostics."
        ),
    }
    output.write_text(json.dumps(metadata, indent=2) + "\n")


def write_residual_operator_report(
    output: Path,
    baselines: dict[float, BaselineResult],
    rows: list[ModeResult],
) -> None:
    publication_rows = publication_mode_rows(rows)
    fundamental_rows = final_fundamental_rows(rows)
    fundamental_zero = next(row for row in fundamental_rows if row.a == 0.0)
    rel_error = abs(fundamental_zero.omega - SCHWARZSCHILD_SCALAR_L2) / abs(SCHWARZSCHILD_SCALAR_L2)

    lines = [
        "# Spectral Residual Operator Report",
        "",
        "This update preserves the baseline time-domain and matrix-pencil workflow, but no longer treats",
        "the waveform-derived matrix-pencil operator as the final QNM operator.",
        "",
        "## Baseline Workflow",
        "",
        "`time-domain evolution -> waveform fitting -> matrix-pencil reduction`",
        "",
        "This path remains useful for diagnostics and for checking the expected deformation trend.",
        "It is not used as the final residual operator.",
        "",
        "## Direct Spectral Residual Workflow",
        "",
        "`KS perturbation equation -> Chebyshev pseudospectral discretization -> P_N(omega) -> R_N(omega)`",
        "",
        "The scalar perturbation equation is discretized directly after compactifying the horizon-to-infinity",
        "domain and factoring the QNM asymptotic behavior. The residual operator",
        "`R_N(omega)=P_N(omega)^dagger P_N(omega)` is Hermitian and positive semidefinite by construction.",
        "Candidate modes are accepted only after residual, convergence, branch-status, and Leaver-style validation checks.",
        "",
        "## Schwarzschild Check",
        "",
        f"- Reference scalar l=2 mode: `{SCHWARZSCHILD_SCALAR_L2.real:.12f} {SCHWARZSCHILD_SCALAR_L2.imag:+.12f}i`.",
        f"- Direct spectral N={FINAL_SPECTRAL_N}: `{fundamental_zero.omega.real:.12f} {fundamental_zero.omega.imag:+.12f}i`.",
        f"- Relative error: `{rel_error:.3e}`.",
        "",
        "## Publication-Safe Spectral Results",
        "",
        "The fundamental branch is reported at N=96. The first-overtone branch is",
        f"reported at the Leaver-validated reference grid N={OVERTONE_PUBLICATION_N};",
        "tracked high-N overtone rows are written separately as exploratory diagnostics.",
        "",
        "| a/M | N | mode | status | baseline fit | matrix pencil | spectral | spectral residual |",
        "|---:|---:|---|---|---:|---:|---:|---:|",
    ]

    for row in publication_rows:
        baseline = baselines[row.a]
        if row.mode == "fundamental":
            fit = f"`{baseline.omega_fit.real:.9f}{baseline.omega_fit.imag:+.9f}i`"
            pencil = f"`{baseline.omega_pencil.real:.9f}{baseline.omega_pencil.imag:+.9f}i`"
        else:
            fit = ""
            pencil = ""
        lines.append(
            f"| {row.a:.1f} | {row.n} | {row.mode} | {row.branch_status} | {fit} | {pencil} | "
            f"`{row.omega.real:.9f}{row.omega.imag:+.9f}i` | "
            f"`{row.omega_residual.real:.9f}{row.omega_residual.imag:+.9f}i` |"
        )

    lines += [
        "",
        "## Conditioning Diagnostics",
        "",
        f"Rows with `cond(P_N) >= {CONDITION_WARNING_THRESHOLD:.1e}` are flagged because",
        "ill-conditioning limits how many digits should be interpreted as physically stable.",
        "",
        "| a/M | N | mode | dim | cond(P_N) | warning | residual norm |",
        "|---:|---:|---|---:|---:|---|---:|",
    ]
    for row in publication_rows:
        lines.append(
            f"| {row.a:.1f} | {row.n} | {row.mode} | {row.matrix_dimension} | "
            f"{row.condition_number:.3e} | {row.conditioning_warning} | {row.residual_norm:.3e} |"
        )

    lines += [
        "",
        "## Verification",
        "",
        "The code verifies Hermiticity, positive semidefiniteness, Schwarzschild recovery,",
        "Leaver agreement, and catalogue consistency through pytest-compatible tests,",
        "`tests/test_qnm_algorithm.py --full`, and `scripts/run_hybrid_qnm_algorithm.py --tests-only`.",
        "",
        "Small negative PSD eigenvalues can appear at the 1e-13 level because of dense Hermitian eigensolver",
        "roundoff. The reported residual norm uses the smallest singular value of `P_N(omega)`, which is",
        "nonnegative by construction.",
        "",
        "High-N first-overtone rows are treated as exploratory until branch-by-branch Leaver or",
        "higher-precision validation is added.",
    ]
    output.write_text("\n".join(lines) + "\n")


def plot_convergence(rows: list[ModeResult], output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4))
    for a in A_VALUES:
        subset = [row for row in rows if row.a == a and row.mode == "fundamental"]
        ns = np.array([row.n for row in subset])
        omegas = np.array([row.omega for row in subset])
        axes[0].plot(ns, omegas.real, marker="o", label=f"a/M={a:g}")
        axes[1].plot(ns, -omegas.imag, marker="o", label=f"a/M={a:g}")
    axes[0].set_xlabel("Chebyshev size N")
    axes[0].set_ylabel(r"$\mathrm{Re}(M\omega)$")
    axes[1].set_xlabel("Chebyshev size N")
    axes[1].set_ylabel(r"$-\mathrm{Im}(M\omega)$")
    for ax in axes:
        ax.grid(alpha=0.25)
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_deformation_trend(
    baselines: dict[float, BaselineResult],
    rows: list[ModeResult],
    output: Path,
) -> None:
    final_fund = final_fundamental_rows(rows)
    final_fund.sort(key=lambda row: row.a)
    a_values = np.array([row.a for row in final_fund])
    td = np.array([baselines[a].omega_fit for a in a_values])
    pencil = np.array([baselines[a].omega_pencil for a in a_values])
    spectral = np.array([row.omega for row in final_fund])
    residual = np.array([row.omega_residual for row in final_fund])

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4))
    series = [
        ("baseline fit", td, "#1f5f8b"),
        ("matrix pencil", pencil, "#669bbc"),
        ("direct spectral", spectral, "#2a9d8f"),
        ("spectral residual", residual, "#d1495b"),
    ]
    for label, values, color in series:
        axes[0].plot(a_values, values.real, marker="o", label=label, color=color)
        axes[1].plot(a_values, -values.imag, marker="o", label=label, color=color)
    axes[0].set_xlabel(r"$a/M$")
    axes[0].set_ylabel(r"$\mathrm{Re}(M\omega)$")
    axes[1].set_xlabel(r"$a/M$")
    axes[1].set_ylabel(r"$-\mathrm{Im}(M\omega)$")
    for ax in axes:
        ax.grid(alpha=0.25)
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def run_pipeline(base_dir: Path) -> None:
    results_dir = base_dir / "outputs" / "results"
    figures_dir = base_dir / "outputs" / "figures"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)

    tests = run_self_tests()
    failed = [test for test in tests if not test.passed]
    if failed:
        raise RuntimeError("Self-tests failed: " + ", ".join(test.name for test in failed))

    baselines = run_baseline(A_VALUES, figures_dir)
    baseline_targets = {a: result.omega_fit for a, result in baselines.items()}
    spectral_rows = run_spectral_study(A_VALUES, SPECTRAL_SIZES, baseline_targets)
    leaver_rows = run_leaver_validation(A_VALUES)
    assert_leaver_validation(leaver_rows)
    catalogue_rows = run_catalogue(A_VALUES)
    assert_catalogue_validation(catalogue_rows)

    write_spectral_results(results_dir / "spectral_results.csv", baselines, spectral_rows)
    write_spectral_results(results_dir / "exploratory_spectral_results.csv", baselines, spectral_rows, publication_only=False)
    write_convergence_table(results_dir / "convergence_results.csv", spectral_rows)
    write_convergence_table(results_dir / "convergence_table.csv", spectral_rows)
    write_resource_estimates(results_dir / "resource_estimates.csv", spectral_rows)
    write_leaver_comparison(results_dir / "leaver_comparison.csv", leaver_rows)
    write_catalogue(results_dir / "qnm_catalogue.csv", catalogue_rows)
    write_catalogue_report(results_dir / "qnm_catalogue_report.md", catalogue_rows)
    physics_outputs = write_physics_analysis(results_dir, figures_dir, catalogue_rows)
    write_residual_operator_report(results_dir / "residual_operator_report.md", baselines, spectral_rows)
    write_run_metadata(results_dir / "run_metadata.json")
    plot_convergence(spectral_rows, figures_dir / "spectral_convergence.png")
    plot_deformation_trend(baselines, spectral_rows, figures_dir / "spectral_deformation_trend.png")
    trajectory_plots = plot_mode_trajectories(catalogue_rows, figures_dir)

    print("Spectral residual QNM pipeline complete.")
    print(f"Baseline module: src/qnm/baseline.py")
    print(f"Spectral module: src/qnm/spectral.py")
    print(f"Leaver validation module: src/qnm/leaver.py")
    print(f"Results directory: {results_dir}")
    print(f"Figures directory: {figures_dir}")
    worst_leaver = max(leaver_rows, key=lambda row: row.relative_difference)
    print(f"Leaver comparison CSV: {results_dir / 'leaver_comparison.csv'}")
    print(f"QNM catalogue CSV: {results_dir / 'qnm_catalogue.csv'}")
    print(f"QNM catalogue report: {results_dir / 'qnm_catalogue_report.md'}")
    print(f"Catalogue physics report: {physics_outputs['report']}")
    print(f"Publication spectral CSV: {results_dir / 'spectral_results.csv'}")
    print(f"Exploratory spectral CSV: {results_dir / 'exploratory_spectral_results.csv'}")
    print(f"Mode trajectory plots: {len(trajectory_plots)} files")
    print(
        "Leaver validation worst relative difference: "
        f"{worst_leaver.relative_difference:.3e} "
        f"(a/M={worst_leaver.a:g}, {worst_leaver.mode}, spectral N={worst_leaver.spectral_n})"
    )
    worst_catalogue = max(catalogue_rows, key=lambda row: row.spectral_leaver_relative_difference)
    print(
        "Catalogue worst spectral/Leaver relative difference: "
        f"{worst_catalogue.spectral_leaver_relative_difference:.3e} "
        f"({worst_catalogue.perturbation_type}, ell={worst_catalogue.ell}, "
        f"n={worst_catalogue.overtone}, a/M={worst_catalogue.a:g})"
    )
    for row in publication_mode_rows(spectral_rows):
        baseline = baselines[row.a]
        baseline_text = ""
        if row.mode == "fundamental":
            baseline_text = (
                f", baseline={baseline.omega_fit.real:.9f}{baseline.omega_fit.imag:+.9f}i, "
                f"pencil={baseline.omega_pencil.real:.9f}{baseline.omega_pencil.imag:+.9f}i"
            )
        print(
            f"a/M={row.a:g}, N={row.n}, {row.mode} [{row.branch_status}]: "
            f"spectral={row.omega.real:.9f}{row.omega.imag:+.9f}i, "
            f"residual={row.omega_residual.real:.9f}{row.omega_residual.imag:+.9f}i"
            f"{baseline_text}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tests-only", action="store_true", help="run verification tests without producing outputs")
    args = parser.parse_args()

    if args.tests_only:
        tests = run_self_tests()
        for test in tests:
            status = "PASS" if test.passed else "FAIL"
            print(f"{status}: {test.name}: value={test.value:.3e}, threshold={test.threshold:.3e}")
        if not all(test.passed for test in tests):
            raise SystemExit(1)
        leaver_rows = run_leaver_validation(A_VALUES)
        for row in leaver_rows:
            status = "PASS" if row.validation_passed else "FAIL"
            print(
                f"{status}: Leaver validation a/M={row.a:g} {row.mode}: "
                f"relative_difference={row.relative_difference:.3e}, "
                f"threshold={row.validation_threshold:.3e}"
            )
        try:
            assert_leaver_validation(leaver_rows)
        except AssertionError as exc:
            print(str(exc))
            raise SystemExit(1) from exc
        catalogue_rows = run_catalogue(A_VALUES)
        for row in catalogue_rows:
            status = "PASS" if row.spectral_leaver_validation_passed else "FAIL"
            print(
                f"{status}: Catalogue validation {row.perturbation_type} "
                f"ell={row.ell} n={row.overtone} a/M={row.a:g}: "
                f"relative_difference={row.spectral_leaver_relative_difference:.3e}"
            )
        try:
            assert_catalogue_validation(catalogue_rows)
        except AssertionError as exc:
            print(str(exc))
            raise SystemExit(1) from exc
        return

    run_pipeline(ROOT_DIR)


if __name__ == "__main__":
    main()
