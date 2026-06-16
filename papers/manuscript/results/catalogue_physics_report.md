# Catalogue Physics Analysis

This report analyzes the Leaver-validated catalogue by comparing each KS
branch against its Schwarzschild endpoint at a/M=0.

## Main Observations

- Largest endpoint fractional frequency shift: `7.23%` for `scalar ell=4 n=0` at a/M=1.
- Scalar ell=2 fundamental endpoint shift: real part `-7.29%`, damping magnitude `-2.59%`.
- Scalar ell=2 fundamental quality-factor shift: `-4.83%`.
- Scalar ell=2 spectroscopic ratio shift `omega0/omega1`: `1.92%`.
- Scalar ell=2 spectroscopic ratio shift `omega0/omega2`: `3.56%`.
- Scalar ell=2 `Re(omega)/[-Im(omega)]` shifts: n=0 `-4.83%`, n=1 `-5.71%`, n=2 `-7.31%`.
- Scalar ell=2 second overtone endpoint shift: real part `-9.16%`, damping magnitude `-1.99%`.
- Axial gravitational ell=2 fundamental endpoint shift: real part `-5.61%`, damping magnitude `-2.63%`.
- Real-part trends across the grid: `decreasing`.
- Damping-magnitude trends across the grid: `decreasing`.
- Worst spectral/Leaver mismatch remains `1.097e-05` for `gravitational ell=2 n=2`.

## Dimensionless Spectroscopic Ratios

The ratio table records mode ratios such as `omega0/omega1`,
`omega0/omega2`, and `Re(omega)/[-Im(omega)]`. These diagnostics are
less directly absorbable into a simple mass rescaling than individual
frequencies.

| branch | ratio | endpoint value | Schwarzschild value | shift |
|---|---|---:|---:|---:|
| scalar ell=2 | omega0/omega1 | 0.823241+0.335659i | 0.836060+0.324207i | 1.92% |
| scalar ell=2 | omega0/omega2 | 0.553881+0.464905i | 0.579814+0.460140i | 3.56% |
| scalar ell=2 | Re/[-Im], n=0 | 4.757253 | 4.998450 | -4.83% |

## Sensitivity Ranking

| rank | branch | endpoint shift | real shift | damping shift | validation max |
|---:|---|---:|---:|---:|---:|
| 1 | scalar ell=4 n=0 | 7.23% | -7.27% | -2.64% | 2.957e-13 |
| 2 | scalar ell=3 n=0 | 7.21% | -7.27% | -2.63% | 3.638e-13 |
| 3 | scalar ell=2 n=0 | 7.17% | -7.29% | -2.59% | 7.914e-13 |
| 4 | scalar ell=4 n=1 | 7.12% | -7.47% | -2.57% | 7.533e-11 |
| 5 | scalar ell=3 n=1 | 7.04% | -7.62% | -2.51% | 1.956e-10 |
| 6 | scalar ell=4 n=2 | 6.91% | -7.89% | -2.44% | 1.225e-08 |
| 7 | gravitational ell=4 n=0 | 6.83% | -6.87% | -2.61% | 3.156e-13 |
| 8 | scalar ell=2 n=1 | 6.83% | -7.95% | -2.38% | 1.814e-09 |

## Interpretation Guardrails

- The scalar sector is the cleanest physics target.
- The axial gravitational rows use a KS-lapse-deformed Regge-Wheeler model;
  they are useful phenomenological diagnostics, not a full gauge-invariant KS
  gravitational perturbation derivation.
- Second overtones carry the largest validation and branch-selection uncertainty.
- The endpoint shifts are catalogue diagnostics, not claims about detectability.
