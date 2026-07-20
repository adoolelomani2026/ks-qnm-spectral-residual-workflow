# QNM Catalogue Report

This low-lying catalogue uses the scalar sector as the least assumption-dependent physics result
and a gauge-invariant odd-parity sector under an explicit frozen-source closure.

## Scope

- Perturbation types: scalar, gravitational.
- Multipoles: ell = 2, 3, 4.
- Modes: fundamental, first overtone, second overtone.
- Confidence tiers: n=0 robust quantitative; n=1 validated low-lying; n=2 exploratory diagnostic.
- Deformations: a/M = 0, 0.2, 0.5, 1.
- Spectral comparison size: N = 32.

For a/M > 0 the gravitational potential is the gauge-invariant frozen-source form
`V=f_a[ell(ell+1)/r^2+2(f_a-1)/r^2-f_a'/r]`.

## Validation

- Worst spectral/Leaver relative difference: `1.202e-05` (gravitational, ell=2, n=2, a/M=0.5).
- Worst Schwarzschild literature relative error: `1.833e-05` (gravitational, ell=4, n=2).

The automated catalogue validation fails if any spectral/Leaver relative
difference exceeds `1.0e-04`. Literature checks use
rounded table tolerances because several source tables report six significant figures.
The Leaver solver is collocation-independent: it uses no Chebyshev grid, matrix-pencil data,
and residual minimization, but it intentionally shares the same perturbation
equation, compact coordinate, endpoint factorization, and potential model.
The continued-fraction residual is reported row-by-row; high-deformation second
overtones can have larger CF residuals because the finite Frobenius recurrence
reduction is less well conditioned there.
