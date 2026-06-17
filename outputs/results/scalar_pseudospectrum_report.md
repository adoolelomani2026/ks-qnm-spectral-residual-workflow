# Scalar Fundamental Pseudospectrum Upgrade Report

## Repository Audit

- Core numerical code lives under `src/qnm/`: compactified geometry helpers,
  Chebyshev spectral operators, Leaver-style validation, catalogue generation,
  normalization utilities, and catalogue-level physics diagnostics.
- Validated generated products live under `outputs/results/` and `outputs/figures/`.
- Manuscript sources live under `papers/manuscript/`.
- The existing pipeline safely supports scalar-sector extensions because the scalar
  Chebyshev and Leaver branches are already cross-validated. The phenomenological
  axial sector was not extended in this upgrade.

## Upgrade Selection

| candidate direction | impact | feasibility | decision |
|---|---|---|---|
| Pseudospectrum and spectral instability | high: connects the residual workflow to QNM instability literature | high: uses existing P_N(omega) and sigma_min machinery | selected |
| Exceptional points or branch interactions | potentially high, but needs denser branch/eigenvector tracking and stronger mathematical evidence | medium | deferred |
| Literature-matched KS benchmarking | useful, but incremental after the Konoplya side comparison | high | secondary future work |
| Gauge-invariant gravitational sector | very high, but requires a new perturbation derivation beyond the current codebase | low for this iteration | deferred |
| Spectroscopy/detectability | useful, but risks unsupported detector claims without a full waveform/noise model | medium | deferred |

## What Was Attempted

A finite-dimensional pseudospectrum diagnostic was added for the scalar
ell=2 fundamental KS QNM branch. For each deformation value, the code
builds the Chebyshev residual operator P_N(omega), centers a local grid
on the residual-minimized spectral frequency, and evaluates
eta(omega)=sigma_min(P_N)/sigma_max(P_N).

## Reproduction Commands

```bash
python scripts/analyze_pseudospectrum.py
python -m pytest
python scripts/run_hybrid_qnm_algorithm.py --tests-only
```

The pseudospectrum command writes the grid, summary table, resolution check,
figures, and this report. The test commands check that the existing validated
pipeline remains intact.

## New Implementation Artifacts

- `src/qnm/pseudospectrum.py`
- `scripts/analyze_pseudospectrum.py`
- `outputs/results/scalar_l2_pseudospectrum_grid.csv`
- `outputs/results/scalar_l2_pseudospectrum_summary.csv`
- `outputs/results/scalar_l2_pseudospectrum_resolution_check.csv`
- `outputs/figures/scalar_l2_pseudospectrum_contours.png`
- `outputs/figures/scalar_l2_pseudospectrum_sensitivity.png`
- `outputs/figures/scalar_l2_pseudospectrum_resolution_check.png`

## What Worked

- Main grid: N=64, grid=81x81, half-width=0.025 in both Re(M omega) and Im(M omega).
- Maximum center-to-Leaver relative difference: 6.816e-12.
- The 10% quantile susceptibility, -Q10(log10 eta), increases by 0.161 from a/M=0 to a/M=1.
- The area fraction satisfying log10(eta)<=-10 grows by a factor 5.06 from Schwarzschild to a/M=1.
- The sign of the Q10 susceptibility trend is stable across N=32, 48, and 64.
- The contour-area diagnostic is secondary to the quantile diagnostic because fixed
  epsilon contour areas depend more strongly on N and on the chosen plotting window.

## Resolution Check

- N=32: Q10 susceptibility gain from a/M=0 to 1 is 0.097.
- N=48: Q10 susceptibility gain from a/M=0 to 1 is 0.130.
- N=64: Q10 susceptibility gain from a/M=0 to 1 is 0.161.

## What Failed Or Was Limited

- The log10(eta)<=-10 contour at a/M=1 touches the local-window boundary: True. The reported area factor is therefore a finite-window diagnostic,
  not a global contour area.
- Absolute contour levels shift with Chebyshev size N, so the finite-N robustness check
  uses the sign and monotonicity of the Q10 susceptibility gain rather than exact equality
  of epsilon-contour areas.
- No exceptional-point search was attempted in this upgrade; mode coalescence would require
  a separate eigenvector and mode-pair condition analysis.

## Publishable Claim

Within the finite-dimensional Chebyshev residual normalization used here,
the scalar ell=2 fundamental KS branch shows increasing local pseudospectral
sensitivity as a/M grows. This complements the frequency-softening result:
the branch moves to lower oscillation frequency while the surrounding
relative-smallest-singular-value basin expands.

## What Remains Speculative

- This is not a proof about the infinite-dimensional KS wave operator.
- Absolute epsilon-contour values depend on N and on the chosen operator normalization.
- The result is local to the scalar ell=2 fundamental branch and should not be
  generalized to overtones or the phenomenological axial sector without separate checks.
- The analysis is not a detector forecast.

## Realistic Journal Target

The new result makes the paper more suitable for a numerics-focused GR journal
or a strong mathematical-physics venue. EPJ Plus and Universe remain realistic;
CQG becomes more plausible if the pseudospectrum section is presented as a
finite-dimensional diagnostic and not as an infinite-dimensional stability theorem.
