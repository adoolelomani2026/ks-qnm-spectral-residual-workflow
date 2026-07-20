# Continued-Fraction Truncation Check

- Taylor order: `96`.
- Root solver: SciPy MINPACK hybrid method with tolerance `1.0e-11`.
- Maximum residual-function evaluations: `1000`.
- Accepted root residual threshold: `1.0e-07`.
- Arithmetic: IEEE-754 double precision.
- Initial guesses: branch-tracked `N=32` Chebyshev eigenvalues.

| case | CF depth | Re(M omega) | Im(M omega) | difference from depth 320 | |CF| |
| --- | ---: | ---: | ---: | ---: | ---: |
| schwarzschild_scalar_l2_n0 | 60 | 0.483643872211 | -0.096758775978 | 2.679e-13 | 4.578e-16 |
| schwarzschild_scalar_l2_n0 | 90 | 0.483643872211 | -0.096758775978 | 4.889e-16 | 7.850e-17 |
| schwarzschild_scalar_l2_n0 | 120 | 0.483643872211 | -0.096758775978 | 0.000e+00 | 4.389e-16 |
| schwarzschild_scalar_l2_n0 | 180 | 0.483643872211 | -0.096758775978 | 0.000e+00 | 4.743e-16 |
| schwarzschild_scalar_l2_n0 | 240 | 0.483643872211 | -0.096758775978 | 0.000e+00 | 4.743e-16 |
| schwarzschild_scalar_l2_n0 | 320 | 0.483643872211 | -0.096758775978 | 0.000e+00 | 4.743e-16 |
| ks_a1_scalar_l2_n0 | 60 | 0.448362409002 | -0.094248179159 | 1.761e-13 | 2.831e-16 |
| ks_a1_scalar_l2_n0 | 90 | 0.448362409002 | -0.094248179159 | 1.347e-15 | 3.925e-16 |
| ks_a1_scalar_l2_n0 | 120 | 0.448362409002 | -0.094248179159 | 3.925e-16 | 4.558e-15 |
| ks_a1_scalar_l2_n0 | 180 | 0.448362409002 | -0.094248179159 | 2.776e-17 | 7.804e-15 |
| ks_a1_scalar_l2_n0 | 240 | 0.448362409002 | -0.094248179159 | 5.551e-17 | 4.066e-15 |
| ks_a1_scalar_l2_n0 | 320 | 0.448362409002 | -0.094248179159 | 0.000e+00 | 3.925e-15 |
| ks_a1_scalar_l2_n1 | 60 | 0.426971489725 | -0.288572913307 | 2.090e-10 | 4.475e-16 |
| ks_a1_scalar_l2_n1 | 90 | 0.426971489574 | -0.288572913453 | 1.115e-12 | 4.847e-16 |
| ks_a1_scalar_l2_n1 | 120 | 0.426971489576 | -0.288572913449 | 2.838e-12 | 5.985e-12 |
| ks_a1_scalar_l2_n1 | 180 | 0.426971489574 | -0.288572913453 | 1.409e-12 | 3.536e-12 |
| ks_a1_scalar_l2_n1 | 240 | 0.426971489574 | -0.288572913453 | 2.052e-12 | 2.589e-12 |
| ks_a1_scalar_l2_n1 | 320 | 0.426971489574 | -0.288572913452 | 0.000e+00 | 9.545e-13 |
