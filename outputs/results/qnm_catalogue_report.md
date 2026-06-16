# QNM Catalogue Report

This is the Leaver-validated catalogue. It extends the scalar Chebyshev spectral
workflow to the axial gravitational Regge-Wheeler-type sector while preserving
the scalar sector.

## Scope

- Perturbation types: scalar, gravitational.
- Multipoles: ell = 2, 3, 4.
- Modes: fundamental, first overtone, second overtone.
- Deformations: a/M = 0, 0.2, 0.5, 1.
- Spectral comparison size: N = 32.

For a/M > 0 the gravitational potential is the KS-lapse-deformed axial
Regge-Wheeler model, `V=f_a(r)[ell(ell+1)/r^2 - 6M/r^3]`.

## Validation

- Worst spectral/Leaver relative difference: `1.097e-05` (gravitational, ell=2, n=2, a/M=0.5).
- Worst Schwarzschild literature relative error: `1.833e-05` (gravitational, ell=4, n=2).

The automated catalogue validation fails if any spectral/Leaver relative
difference exceeds `1.0e-04`. Literature checks use
rounded table tolerances because several source tables report six significant figures.
The Leaver solver is independent of Chebyshev collocation, matrix-pencil data,
and residual minimization, but it intentionally shares the same perturbation
equation, compact coordinate, endpoint factorization, and potential model.
The continued-fraction residual is reported row-by-row; high-deformation second
overtones can have larger CF residuals because the finite Frobenius recurrence
reduction is less well conditioned there.
