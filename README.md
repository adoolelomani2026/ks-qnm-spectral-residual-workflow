# Chebyshev-Leaver Spectral Residual Workflow for KS Quasinormal Modes

This project computes quasinormal modes (QNMs) of Kazakov-Solodukhin (KS)
quantum-deformed Schwarzschild black holes with a direct Chebyshev spectral
workflow, Leaver-style continued-fraction validation, and fixed-mass
spectroscopy diagnostics. It also includes a finite-dimensional
pseudospectrum diagnostic for the validated scalar fundamental branch.

The revision's strategic scope is:

> The scalar low-lying modes are the least assumption-dependent physics result;
> the axial sector is gauge invariant under an explicit frozen-source closure.
> The numerical contribution is auditability and validation, not a new
> high-overtone solver.

Repository URL:
<https://github.com/adoolelomani2026/ks-qnm-spectral-residual-workflow>

## What Is Implemented

- A scalar time-domain baseline workflow:
  `time-domain evolution -> ringdown fit -> matrix-pencil diagnostic`.
- A direct Chebyshev pseudospectral solver for the compactified QNM equation.
- A quadratic polynomial eigenvalue problem:

  ```text
  P_N(omega) u = (A0 + omega A1 + omega^2 A2) u = 0
  ```

- A Hermitian residual operator:

  ```text
  R_N(omega) = P_N(omega)^dagger P_N(omega)
  ```

- A Leaver-style Frobenius continued-fraction validation layer.
- Scalar and gauge-invariant inverse-Cowling axial catalogues for `ell = 2, 3, 4` and
  overtone indices `n = 0, 1, 2`.
- Dimensionless spectroscopic-ratio diagnostics such as `omega0/omega1`,
  `omega0/omega2`, and `Re(omega)/[-Im(omega)]`.
- Scalar `ell=2,n=0` pseudospectrum diagnostics based on
  `eta_N = sigma_min(P_N)/sigma_max(P_N)`.
- Literature positioning against Konoplya (2020),
  Bolokhov-Bronnikov-Konoplya (2025): earlier KS work established QNM
  deformation and overtone sensitivity, while this project adds fixed-`M`
  Chebyshev-Leaver validation, dimensionless ratios, and quality-factor
  shifts. A normalization-matched scalar `ell=0` side comparison with
  Konoplya's fixed-horizon table is generated separately.
- Direct comparison of 18 overlapping conditional axial modes with the public
  2026 Batic-Dutykh-Sukaiti high-precision Chebyshev calculation.

The catalogue Leaver configuration is explicit and reproducible: Taylor
coefficient order `96`, Gaussian-elimination/continued-fraction depth `240`,
double-precision arithmetic, SciPy's MINPACK hybrid root finder (`hybr`), root
tolerance `1e-11`, at most `1000` function evaluations, and continuation in
`a/M` from branch-matched Schwarzschild seeds. The default fundamental-mode
continued-fraction residual threshold is `1e-7`; the deliberately cautious
catalogue threshold is `1e-4` because it includes exploratory `n=2` rows.
These settings are stored in `config/leaver_revision.json`. The generalized-
eigenvalue filtering, equilibration, clustering, continuation, refinement, and
acceptance settings are stored in `config/spectral_selection.json`.

## Repository Map

- `src/qnm/` - importable QNM solver package.
- `src/qnm/common.py` - shared constants, KS metric functions, potentials, and mode selection.
- `src/qnm/baseline.py` - finite-difference waveform evolution, ringdown fitting, matrix-pencil diagnostic.
- `src/qnm/spectral.py` - Chebyshev collocation, generalized eigenvalue solve, and residual diagnostics.
- `src/qnm/leaver.py` - Frobenius recurrence and continued-fraction validation.
- `src/qnm/catalogue.py` - scalar/gravitational catalogue generation and the manuscript trajectory plot.
- `src/qnm/analysis.py` - Schwarzschild-relative catalogue physics diagnostics.
- `src/qnm/pseudospectrum.py` - scalar finite-`N` pseudospectrum grids,
  quantile diagnostics, contour-area estimates, and resolution checks.
- `src/qnm/normalization.py` - fixed-horizon/fixed-mass conversion helpers for
  quantitative literature comparisons.
- `scripts/` - command-line entry points.
- `data/literature/` - transcribed-literature CSV files for normalization-matched comparison tables.
- `tests/` - pytest-compatible validation checks.
- `outputs/results/` - generated CSV tables and selected Markdown summaries.
- `outputs/figures/` - manuscript-supporting convergence, catalogue, trajectory, and pseudospectrum figures.
- `papers/manuscript/` - current manuscript TeX/PDF; figures are loaded from `outputs/figures/`.
- `docs/` - literature-normalization protocol.

## Installation

The current local environment used for validation was:

```text
Python 3.14.3
numpy 2.4.3
scipy 1.17.1
matplotlib 3.10.8
```

Create and activate a virtual environment, then install the pinned runtime
dependencies and the local package in editable mode:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
pip install -e .
```

For pytest-based development checks, install:

```powershell
pip install -r requirements-dev.txt
pip install -e .
```

After installation, a quick import check should succeed:

```powershell
python -c "import qnm; print(qnm.__file__)"
```

## Validation

Run the full validation suite:

```powershell
python tests/test_qnm_algorithm.py --full
```

This checks:

- Hermiticity of `R_N`.
- Numerical consistency with the algebraic Hermiticity and positive
  semidefiniteness of `R_N`.
- Schwarzschild scalar reference recovery.
- Leaver/spectral agreement for scalar validation cases.
- Catalogue-level scalar and inverse-Cowling axial validation.

The default script mode and default pytest mode run only fast checks:

```powershell
python tests/test_qnm_algorithm.py
python -m pytest
```

The full Leaver/catalogue validation is computationally nontrivial. On the
current machine, the full test run took about 2.5 minutes.

## Regenerating Outputs

Run the full pipeline:

```powershell
python scripts/run_hybrid_qnm_algorithm.py
```

Run only verification from the pipeline driver:

```powershell
python scripts/run_hybrid_qnm_algorithm.py --tests-only
```

Regenerate just the catalogue:

```powershell
python scripts/run_catalogue.py
```

Analyze catalogue-level physics trends from the generated catalogue:

```powershell
python scripts/analyze_catalogue_physics.py
```

Reproduce the normalization-matched Konoplya (2020) scalar `ell=0` side
comparison:

```powershell
python scripts/compare_konoplya2020_scalar_l0.py
```

Run the scalar `ell=2,n=0` endpoint `N=128` spot check:

```powershell
python scripts/check_n128_spot.py
```

This writes:

- `outputs/results/n128_spot_check.csv`
- `outputs/results/n128_spot_check_report.md`

Generate the continued-fraction convergence table:

```powershell
python scripts/analyze_leaver_convergence.py
```

Generate Taylor-order and Chebyshev-to-continued-fraction convergence audits:

```powershell
python scripts/analyze_solver_convergence.py
```

Generate the axial potential, positivity, growing-root, and external-public-data audits:

```powershell
python scripts/audit_axial_model.py
```

Generate the scalar-potential peak analysis:

```powershell
python scripts/analyze_scalar_potential_peaks.py
```

Generate the line-numbered marked manuscript by comparing the revision with
submitted commit `9090c8c`:

```powershell
python scripts/build_marked_manuscript.py
```

Compile the generated TeX with `markedrevision` defined to color only changed
lines and enable line numbers.

Regenerate the scalar `ell=2,n=0` pseudospectrum grids, summaries, report, and
figures:

```powershell
python scripts/analyze_pseudospectrum.py
python scripts/audit_pseudospectrum_robustness.py
```

This writes:

- `outputs/results/scalar_l2_pseudospectrum_grid.csv`
- `outputs/results/scalar_l2_pseudospectrum_summary.csv`
- `outputs/results/scalar_l2_pseudospectrum_resolution_check.csv`
- `outputs/results/scalar_pseudospectrum_report.md`
- `outputs/figures/scalar_l2_pseudospectrum_contours.png`
- `outputs/figures/scalar_l2_pseudospectrum_sensitivity.png`
- `outputs/figures/scalar_l2_pseudospectrum_resolution_check.png`

For a clean submission-facing regeneration, run the main pipeline, the scalar
pseudospectrum diagnostic, the normalization-matched literature comparison, and
the automated tests:

```powershell
python scripts/run_hybrid_qnm_algorithm.py
python scripts/check_n128_spot.py
python scripts/analyze_leaver_convergence.py
python scripts/analyze_solver_convergence.py
python scripts/analyze_scalar_potential_peaks.py
python scripts/audit_axial_model.py
python scripts/analyze_pseudospectrum.py
python scripts/audit_pseudospectrum_robustness.py
python scripts/compare_konoplya2020_scalar_l0.py
pytest
pytest -m slow
```

## Scientific Scope

The scalar sector is the least assumption-dependent physics target.  For the
axial sector, the repository adopts an additional inverse-Cowling closure of
the general gauge-invariant sourced odd-parity equations: the background
effective stress is retained while both odd source amplitudes are set to zero.

```text
V_ax,ell = f_a [ell(ell+1)/r^2 + 2(f_a-1)/r^2 - f_a'/r].
```

This reduces to the standard Schwarzschild Regge--Wheeler potential at `a=0`
and differs from simple lapse substitution when `a>0`.  Gauge invariance
removes coordinate artifacts from the metric variable; it does not make this
closure the unique nonspherical perturbation theory of the spherical KS
effective action.

## Interpretation Boundaries

- Scalar fundamental modes are the safest physics claim in the current project.
- Axial entries are conditional odd-parity modes of the gauge-invariant
  inverse-Cowling effective-source model, not predictions of an unspecified
  quantum nonspherical completion.
- Catalogue shifts are not observational forecasts; observability would require
  waveform modeling, detector-noise weighting, and parameter-degeneracy studies.
- First overtones are cross-validated low-lying modes on the Leaver-checked `N=32`
  grid. Second overtones are exploratory branch diagnostics; higher overtones
  are outside the revision's scope.
- Pseudospectrum contours are finite-dimensional Chebyshev diagnostics, not
  proofs about the infinite-dimensional KS wave operator. Absolute contour
  levels depend on `N` and operator normalization.

Claim hierarchy:

| Sector or branch | Status | Use |
|---|---|---|
| Scalar `n=0` fundamentals | Robust quantitative | Main physics claim |
| Scalar `n=1` first overtones | Validated low-lying at `N=32` | Secondary catalogue and ratio diagnostics |
| `n=2` overtones | Exploratory/diagnostic | Branch diagnostics only |
| Axial rows | Gauge invariant, closure dependent | Conditional modes with frozen axial effective-source current |

## Current Numerical Highlights

- The external Schwarzschild scalar `ell=2,n=0` benchmark is
  `M omega = 0.483643872211 - 0.096758775978 i` (Cavalcante and Carneiro da
  Cunha, 2021). The regenerated `N=96` spectral and Leaver differences are
  about `5.60e-11` and `8.24e-13`, respectively.
- Catalogue spectral/Leaver disagreement is worst for second overtones and
  is currently `1.202e-05`, below the project threshold of `1e-4`.
- Catalogue-level physics diagnostics show successive decreases across the sampled grid in both
  oscillation frequency and damping magnitude across the validated deformation
  grid. For the scalar `ell=2` fundamental branch at `a/M=1`, the shifts are
  `-7.29%` in `Re(M omega)`, `-2.59%` in `-Im(M omega)`, and `-4.83%`
  in the dimensionless quality factor.
- Dimensionless scalar `ell=2` spectroscopy diagnostics shift by `1.92%`
  for `omega0/omega1`, `3.56%` for `omega0/omega2`, and `-4.83%`
  for the fundamental `Re(omega)/[-Im(omega)]` ratio at `a/M=1`.
- Scalar `ell=2,n=0` pseudospectrum diagnostics at `N=64` show increasing
  local finite-`N` sensitivity with deformation: the 10% quantile
  susceptibility `-Q10(log10 eta_N)` increases by `0.161` from `a/M=0` to
  `a/M=1`, and the `log10 eta_N <= -10` local area grows by a factor `5.06`
  within the chosen window.
- Reported first-overtone rows are frozen at the Leaver-validated
  `N=32` grid; tracked high-`N` overtone rows are kept in
  `outputs/results/exploratory_spectral_results.csv`.

The exact environment for regenerated outputs is recorded in
`outputs/results/run_metadata.json`.
