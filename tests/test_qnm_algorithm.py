#!/usr/bin/env python3
"""Pytest-compatible validation for the KS QNM spectral-residual pipeline.

Default pytest runs are intentionally fast. The Leaver and catalogue checks are
marked ``slow`` and can be run with:

    pytest -m slow

The module can still be executed as a script. Script mode runs the fast checks
by default and accepts ``--full`` for the slow validation pass.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

try:
    import pytest
except ImportError:  # pragma: no cover - only used outside the test environment.
    pytest = None

from qnm.analysis import compute_branch_summaries, compute_spectroscopic_ratios, compute_trend_rows
from qnm.catalogue import CatalogueRow, assert_catalogue_validation, run_catalogue
from qnm.common import SCHWARZSCHILD_SCALAR_L2
from qnm.leaver import assert_leaver_validation, run_leaver_validation
from qnm.normalization import (
    alpha_from_horizon_beta,
    fixed_horizon_to_fixed_mass,
    fixed_mass_to_fixed_horizon,
    horizon_beta_from_alpha,
    horizon_over_mass_from_alpha,
    mass_over_horizon_from_alpha,
)
from qnm.pseudospectrum import compute_pseudospectrum_grid, summarize_grid
from qnm.prl_instability import active_branches
from qnm.spectral import (
    build_spectral_problem,
    generalized_eigenvalues,
    run_self_tests,
    select_physical_mode,
)
from qnm.universality import MODELS, build_model_spectral_problem, horizon_bardeen, horizon_hayward


def _slow_marker(func):
    if pytest is None:
        return func
    return pytest.mark.slow(func)


def _assert_self_tests_pass() -> None:
    tests = run_self_tests()
    failed = [test for test in tests if not test.passed]
    assert not failed, ", ".join(
        f"{test.name}: value={test.value:.3e}, threshold={test.threshold:.3e}" for test in failed
    )


def test_core_self_tests_pass() -> None:
    _assert_self_tests_pass()


def test_scaled_spectral_solver_recovers_schwarzschild_reference() -> None:
    problem = build_spectral_problem(0.0, 32)
    values = generalized_eigenvalues(problem)
    omega = select_physical_mode(values, SCHWARZSCHILD_SCALAR_L2)
    relative_error = abs(omega - SCHWARZSCHILD_SCALAR_L2) / abs(SCHWARZSCHILD_SCALAR_L2)
    assert relative_error < 1.0e-8


def test_catalogue_physics_analysis_tracks_endpoint_shifts() -> None:
    rows = [
        CatalogueRow(
            perturbation_type="scalar",
            ell=2,
            overtone=0,
            mode="fundamental",
            a=0.0,
            spectral_n=32,
            omega_spectral=1.0 - 0.1j,
            omega_leaver=1.0 - 0.1j,
            spectral_leaver_relative_difference=1.0e-12,
            continued_fraction_abs=1.0e-15,
            schwarzschild_literature=1.0 - 0.1j,
            literature_relative_error=0.0,
            literature_tolerance=1.0e-5,
            spectral_leaver_validation_passed=True,
            literature_validation_passed=True,
        ),
        CatalogueRow(
            perturbation_type="scalar",
            ell=2,
            overtone=0,
            mode="fundamental",
            a=1.0,
            spectral_n=32,
            omega_spectral=0.9 - 0.09j,
            omega_leaver=0.9 - 0.09j,
            spectral_leaver_relative_difference=2.0e-12,
            continued_fraction_abs=2.0e-15,
            schwarzschild_literature=None,
            literature_relative_error=None,
            literature_tolerance=None,
            spectral_leaver_validation_passed=True,
            literature_validation_passed=None,
        ),
        CatalogueRow(
            perturbation_type="scalar",
            ell=2,
            overtone=1,
            mode="first_overtone",
            a=0.0,
            spectral_n=32,
            omega_spectral=0.8 - 0.2j,
            omega_leaver=0.8 - 0.2j,
            spectral_leaver_relative_difference=3.0e-12,
            continued_fraction_abs=3.0e-15,
            schwarzschild_literature=0.8 - 0.2j,
            literature_relative_error=0.0,
            literature_tolerance=1.0e-5,
            spectral_leaver_validation_passed=True,
            literature_validation_passed=True,
        ),
        CatalogueRow(
            perturbation_type="scalar",
            ell=2,
            overtone=1,
            mode="first_overtone",
            a=1.0,
            spectral_n=32,
            omega_spectral=0.7 - 0.18j,
            omega_leaver=0.7 - 0.18j,
            spectral_leaver_relative_difference=4.0e-12,
            continued_fraction_abs=4.0e-15,
            schwarzschild_literature=None,
            literature_relative_error=None,
            literature_tolerance=None,
            spectral_leaver_validation_passed=True,
            literature_validation_passed=None,
        ),
    ]

    trend_rows = compute_trend_rows(rows)
    summaries = compute_branch_summaries(rows)
    ratio_rows = compute_spectroscopic_ratios(rows)

    endpoint_fundamental = next(row for row in trend_rows if row.overtone == 0 and row.a == 1.0)
    assert endpoint_fundamental.fractional_real_shift == pytest.approx(-0.1)
    assert endpoint_fundamental.fractional_damping_shift == pytest.approx(-0.1)
    assert summaries[0].real_trend == "decreasing"
    assert summaries[0].damping_trend == "decreasing"
    assert summaries[0].max_spectral_leaver_relative_difference == pytest.approx(2.0e-12)
    omega_ratio = next(row for row in ratio_rows if row.ratio_name == "omega0_over_omega1" and row.a == 1.0)
    assert omega_ratio.fractional_shift_abs > 0.0
    real_to_damping = next(row for row in ratio_rows if row.ratio_name == "real_to_damping_n0" and row.a == 1.0)
    assert real_to_damping.fractional_real_shift == pytest.approx(0.0)


def test_literature_normalization_roundtrip() -> None:
    alpha = 1.0
    omega_m = 0.448362409 - 0.094248179j

    beta, omega_rh = fixed_mass_to_fixed_horizon(alpha, omega_m)
    alpha_roundtrip, omega_m_roundtrip = fixed_horizon_to_fixed_mass(beta, omega_rh)

    assert beta == pytest.approx(1.0 / 5.0**0.5)
    assert alpha_from_horizon_beta(beta) == pytest.approx(alpha)
    assert horizon_beta_from_alpha(alpha) == pytest.approx(beta)
    assert horizon_over_mass_from_alpha(alpha) == pytest.approx(5.0**0.5)
    assert mass_over_horizon_from_alpha(alpha) == pytest.approx(1.0 / 5.0**0.5)
    assert alpha_roundtrip == pytest.approx(alpha)
    assert omega_m_roundtrip == pytest.approx(omega_m)


def test_pseudospectrum_grid_smoke() -> None:
    grid = compute_pseudospectrum_grid(
        a=0.0,
        leaver_target=SCHWARZSCHILD_SCALAR_L2,
        spectral_n=16,
        grid_size=7,
        half_width=0.01,
    )
    summary = summarize_grid(grid, thresholds=(-7.0,), quantiles=(0.10, 0.50))

    assert grid.log10_relative_sigma.shape == (7, 7)
    assert np.isfinite(grid.log10_relative_sigma).all()
    assert summary.min_log10_relative_sigma <= summary.median_log10_relative_sigma
    assert summary.quantiles[0.10] <= summary.quantiles[0.50]
    assert summary.threshold_areas[-7.0] >= 0.0


def test_prl_instability_scan_branch_policy() -> None:
    branches = active_branches()
    assert (0, 0) in branches
    assert (0, 1) in branches
    assert (0, 2) not in branches
    assert (4, 2) in branches


def test_hayward_comparator_spectral_problem_smoke() -> None:
    assert horizon_hayward(0.5) > 0.0
    problem = build_model_spectral_problem(MODELS["hayward"], parameter=0.2, n=12, ell=2)
    assert problem.a0.shape == (12, 12)
    assert problem.left.shape == (24, 24)
    assert np.isfinite(problem.r_nodes).all()


def test_bardeen_comparator_spectral_problem_smoke() -> None:
    assert horizon_bardeen(0.5) > 0.0
    problem = build_model_spectral_problem(MODELS["bardeen"], parameter=0.2, n=12, ell=2)
    assert problem.a0.shape == (12, 12)
    assert problem.left.shape == (24, 24)
    assert np.isfinite(problem.r_nodes).all()


@_slow_marker
def test_leaver_validation_reference_grid() -> None:
    rows = run_leaver_validation(a_values=[0.0])
    assert_leaver_validation(rows)


@_slow_marker
def test_catalogue_scalar_reference_grid() -> None:
    rows = run_catalogue(a_values=[0.0], perturbation_types=["scalar"], ell_values=[2])
    assert_catalogue_validation(rows)


def run_fast_checks() -> None:
    _assert_self_tests_pass()
    test_scaled_spectral_solver_recovers_schwarzschild_reference()
    print("PASS: fast spectral-residual checks")


def run_full_checks() -> None:
    run_fast_checks()

    leaver_rows = run_leaver_validation()
    assert_leaver_validation(leaver_rows)
    worst_leaver = max(leaver_rows, key=lambda row: row.relative_difference)
    print(
        "PASS: Leaver validation: "
        f"worst_relative_difference={worst_leaver.relative_difference:.3e} "
        f"(a/M={worst_leaver.a:g}, {worst_leaver.mode})"
    )

    catalogue_rows = run_catalogue()
    assert_catalogue_validation(catalogue_rows)
    worst_catalogue = max(catalogue_rows, key=lambda row: row.spectral_leaver_relative_difference)
    print(
        "PASS: Catalogue scalar/gravitational validation: "
        f"worst_relative_difference={worst_catalogue.spectral_leaver_relative_difference:.3e} "
        f"({worst_catalogue.perturbation_type}, ell={worst_catalogue.ell}, "
        f"n={worst_catalogue.overtone}, a/M={worst_catalogue.a:g})"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true", help="run slow Leaver and catalogue validation")
    args = parser.parse_args()

    if args.full:
        run_full_checks()
    else:
        run_fast_checks()
        print("Use `python tests/test_qnm_algorithm.py --full` for Leaver/catalogue validation.")


if __name__ == "__main__":
    main()
