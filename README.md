# Chebyshev-Leaver Spectral Residual Workflow for KS Quasinormal Modes

This project computes quasinormal modes (QNMs) of Kazakov-Solodukhin (KS)
quantum-deformed Schwarzschild black holes with a direct Chebyshev spectral
workflow, Leaver-style continued-fraction validation, and fixed-mass
spectroscopy diagnostics. It also includes a finite-dimensional
pseudospectrum diagnostic for the validated scalar fundamental branch.

The strongest current interpretation is:

> A Chebyshev-Leaver spectral residual workflow for KS black-hole
> quasinormal-mode spectroscopy, with explicit branch-status discipline and
> dimensionless catalogue diagnostics.

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
- Scalar and axial gravitational catalogues for `ell = 2, 3, 4` and
  overtone indices `n = 0, 1, 2`.
- Dimensionless spectroscopic-ratio diagnostics such as `omega0/omega1`,
  `omega0/omega2`, and `Re(omega)/[-Im(omega)]`.
- Scalar `ell=2,n=0` pseudospectrum diagnostics based on
  `eta_N = sigma_min(P_N)/sigma_max(P_N)`.
- Literature positioning against Konoplya (2020) and
  Bolokhov-Bronnikov-Konoplya (2025): earlier KS work established QNM
  deformation and overtone sensitivity, while this project adds fixed-`M`
  Chebyshev-Leaver validation, dimensionless ratios, and quality-factor
  shifts. A normalization-matched scalar `ell=0` side comparison with
  Konoplya's fixed-horizon table is generated separately.

## Repository Map

- `src/qnm/` - importable QNM solver package.
- `src/qnm/common.py` - shared constants, KS metric functions, potentials, and mode selection.
- `src/qnm/baseline.py` - finite-difference waveform evolution, ringdown fitting, matrix-pencil diagnostic.
- `src/qnm/spectral.py` - Chebyshev collocation, generalized eigenvalue solve, and residual diagnostics.
- `src/qnm/leaver.py` - Frobenius recurrence and continued-fraction validation.
- `src/qnm/catalogue.py` - scalar/gravitational catalogue generation and trajectory plots.
- `src/qnm/analysis.py` - Schwarzschild-relative catalogue physics diagnostics.
- `src/qnm/pseudospectrum.py` - scalar finite-`N` pseudospectrum grids,
  quantile diagnostics, contour-area estimates, and resolution checks.
- `src/qnm/normalization.py` - fixed-horizon/fixed-mass conversion helpers for
  quantitative literature comparisons.
- `scripts/` - command-line entry points.
- `data/literature/` - transcribed-literature CSV files for normalization-matched comparison tables.
- `tests/` - pytest-compatible validation checks.
- `outputs/results/` - generated CSV tables and Markdown reports.
- `outputs/figures/` - generated waveform, convergence, deformation, sensitivity, and trajectory plots.
- `papers/manuscript/` - current manuscript TeX/PDF; figures are loaded from `outputs/figures/`.
- `papers/followup/` - outline for the next physics-first spectroscopy paper.
- `docs/` - literature-normalization protocol and project notes.

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
- Positive semidefiniteness of `R_N` up to numerical roundoff.
- Schwarzschild scalar reference recovery.
- Leaver/spectral agreement for scalar validation cases.
- Catalogue-level scalar and axial gravitational validation.

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

Regenerate the scalar `ell=2,n=0` pseudospectrum grids, summaries, report, and
figures:

```powershell
python scripts/analyze_pseudospectrum.py
```

This writes:

- `outputs/results/scalar_l2_pseudospectrum_grid.csv`
- `outputs/results/scalar_l2_pseudospectrum_summary.csv`
- `outputs/results/scalar_l2_pseudospectrum_resolution_check.csv`
- `outputs/results/scalar_pseudospectrum_report.md`
- `outputs/figures/scalar_l2_pseudospectrum_contours.png`
- `outputs/figures/scalar_l2_pseudospectrum_sensitivity.png`
- `outputs/figures/scalar_l2_pseudospectrum_resolution_check.png`

Run the broader PRL-level scalar instability stress test:

```powershell
python scripts/run_prl_instability_scan.py
```

This scan tests whether the stronger claim "KS deformation universally
amplifies scalar QNM pseudospectral sensitivity" survives lower multipoles,
overtones, denser deformation sampling, endpoint Leaver checks, finite-`N`
checks through `N=96` for fundamentals, window-size checks, and mode-pair
diagnostics. The current output is deliberately conservative: it recommends
PRD/CQG rather than PRL.

Main outputs:

- `outputs/results/prl_instability_assessment.md`
- `outputs/results/prl_instability_branch_verdicts.csv`
- `outputs/results/prl_instability_scan_summary.csv`
- `outputs/results/prl_instability_endpoint_leaver_checks.csv`
- `outputs/results/prl_instability_mode_pair_diagnostics.csv`
- `outputs/results/prl_instability_barrier_metrics.csv`
- `outputs/results/prl_instability_window_sensitivity.csv`
- `outputs/figures/prl_instability_central_heatmap.png`
- `outputs/figures/prl_instability_robustness.png`
- `outputs/figures/prl_instability_softening_scatter.png`
- `papers/prl/ks_instability_prl_letter_draft.md`
- `papers/prl/supplemental_material_outline.md`

Classify the KS scalar stress-test branches by pseudospectral sensitivity:

```powershell
python scripts/analyze_branch_sensitivity.py
```

This joins the existing stress-test branch verdicts, endpoint Leaver checks,
barrier metrics, condition indicators, and finite-`N` spread diagnostics. It
does not assume universal amplification; it separates robustly increasing,
weakly increasing, non-monotonic, and numerically inconclusive branches.

Main outputs:

- `outputs/results/branch_sensitivity_classes.csv`
- `outputs/results/branch_sensitivity_correlations.csv`
- `outputs/results/branch_sensitivity_predictor_search.csv`
- `outputs/results/branch_sensitivity_report.md`
- `outputs/figures/branch_sensitivity_predictors.png`

Run the model-zoo universality test:

```powershell
python scripts/run_universality_scan.py
```

This compares the KS deformation with Hayward and Bardeen regular black holes
using the same finite-dimensional Chebyshev residual diagnostics. The result is
a direct falsification of the broadest PRL-style claim: KS shows positive
endpoint susceptibility gains on the usable scalar branches, while the Hayward
and Bardeen scans show negative endpoint gains across the tested branches. The
current assessment is therefore that the effect is a deformation-response
taxonomy rather than a generic quantum-correction principle.

Main outputs:

- `outputs/results/universality_assessment.md`
- `outputs/results/universality_branch_verdicts.csv`
- `outputs/results/universality_model_scan.csv`
- `outputs/results/universality_barrier_metrics.csv`
- `outputs/results/universality_correlations.csv`
- `outputs/results/universality_mode_pair_diagnostics.csv`
- `outputs/figures/universality_endpoint_gain_heatmap.png`
- `outputs/figures/universality_softening_vs_sensitivity.png`
- `outputs/figures/universality_barrier_correlation.png`
- `outputs/figures/universality_width_correlation.png`

## Scientific Scope

The scalar sector is the cleanest physics target. The axial gravitational sector
uses a KS-lapse-deformed Regge-Wheeler potential,

```text
V_RW,l(r; a) = f_a(r) [ l(l+1)/r^2 - 6M/r^3 ],
```

which reduces to the standard Schwarzschild Regge-Wheeler model at `a = 0`.
This is useful for validation and catalogue exploration, but it is not yet a
complete gauge-invariant treatment of gravitational perturbations in the KS
spacetime.

## Interpretation Boundaries

- Scalar fundamental modes are the safest physics claim in the current project.
- Axial gravitational entries are phenomenological KS-lapse-deformed
  Regge-Wheeler diagnostics, not final gauge-invariant KS gravitational
  predictions.
- Catalogue shifts are not detectability forecasts; observability would require
  waveform modeling, detector-noise weighting, and parameter-degeneracy studies.
- Overtones remain the least robust numerical sector. First overtones are
  publication-facing only on the Leaver-validated `N=32` grid, while `n=2`
  branches should be treated as branch diagnostics until further stress-tested.
- Pseudospectrum contours are finite-dimensional Chebyshev diagnostics, not
  proofs about the infinite-dimensional KS wave operator. Absolute contour
  levels depend on `N` and operator normalization.

Claim hierarchy:

| Sector or branch | Status | Use |
|---|---|---|
| Scalar `n=0` fundamentals | Strongest validated sector | Main physics claim |
| Scalar `n=1` first overtones | Leaver-validated at `N=32` | Secondary catalogue and ratio diagnostics |
| `n=2` overtones | Informative but most delicate | Branch diagnostics only |
| Axial gravitational rows | Phenomenological KS-lapse-deformed Regge-Wheeler model | Validation and trend comparison only |

## Current Numerical Highlights

- Direct spectral Schwarzschild scalar `ell=2` fundamental relative error:
  `1.462e-10` at `N=96` in the current regenerated outputs.
- Leaver-style Schwarzschild scalar fundamental relative error:
  about `1e-10`.
- Catalogue spectral/Leaver disagreement is worst for second overtones and
  is currently `1.097e-05`, below the project threshold of `1e-4`.
- Catalogue-level physics diagnostics show monotonic decreases in both
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
- A broader PRL-level scalar scan finds positive endpoint susceptibility gains
  for the active oscillatory branches, but only `8/14` branches satisfy the
  stricter usable/supporting criteria. The current verdict is therefore
  PRD/CQG-strength, not PRL-strength.
- Branch-dependence analysis classifies those active scalar branches into
  `8` robustly increasing, `3` weakly increasing, `1` non-monotonic, and `2`
  numerically inconclusive cases. The strongest numerical correlate of gain
  magnitude is condition-indicator growth, while the empirical overtone-load
  separator `n/(ell+1/2) <= 2/3` identifies the reliable increasing cases in
  the current finite sample. This is not a derived scaling law.
- A KS/Hayward/Bardeen model-zoo universality scan finds only `8/30`
  model-branch verdicts satisfy the stricter support criteria. KS contributes
  all eight supporting branches, while Hayward and Bardeen each have `0/10`
  supporting branches and negative endpoint susceptibility gains on every
  tested branch. The broad "generic quantum amplification" claim is therefore
  not supported.
- Publication-facing first-overtone rows are frozen at the Leaver-validated
  `N=32` grid; tracked high-`N` overtone rows are kept in
  `outputs/results/exploratory_spectral_results.csv`.

The exact environment for regenerated outputs is recorded in
`outputs/results/run_metadata.json`.

## Next Research Steps

The next stage should focus on physics and robustness rather than adding more
infrastructure:

- Develop `papers/followup/spectroscopic_signatures_outline.md` into a scalar-first physics paper using the generated dimensionless ratio diagnostics, the normalization protocol in `docs/LITERATURE_NORMALIZATION_PROTOCOL.md`, the Konoplya side comparison in `scripts/compare_konoplya2020_scalar_l0.py`, and the scaffold in `scripts/prepare_literature_comparison.py`.
- Axial-polar or gauge-invariant gravitational perturbation analysis for the KS background.
- Dedicated overtone branch tracking, especially for `n=2`.
- Conditioning/scaling improvements for high-resolution generalized eigenvalue problems.
- Extend pseudospectrum diagnostics to scalar overtones only after stronger branch tracking is in place.
- Compare the finite-dimensional KS pseudospectrum against other quantum-corrected Schwarzschild models under matched normalization and discretization choices.
