# N=128 Scalar Fundamental Spot Check

This report tests the scalar `ell=2,n=0` fundamental branch at the
Schwarzschild and endpoint-deformed cases, `a/M=0` and `a/M=1`.
It is a robustness check for the high-condition-number spectral pencils,
not a claim of additional significant digits.

## Endpoint Differences

| a/M | |omega_128 - omega_96| | relative difference | |omega_128 - omega_Leaver| |
| --- | ---: | ---: | ---: |
| 0 | 6.949e-11 | 1.409e-10 | 6.867e-11 |
| 1 | 4.639e-10 | 1.013e-09 | 4.699e-10 |

## Interpretation

- The `N=128` endpoint frequencies remain within `5e-10` of the
  `N=96` publication-facing values for the tested scalar fundamental branch.
- The movement is small compared with the displayed table precision and does
  not change any catalogue trend, percent shift, or physics conclusion.
- The `N=96 -> N=128` movement is not used as a claim of extra digits because
  the high-`N` sequence is on a double-precision plateau for these
  ill-conditioned compactified pencils.
