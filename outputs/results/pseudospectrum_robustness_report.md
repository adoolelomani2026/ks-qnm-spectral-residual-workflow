# Pseudospectrum robustness audit

The endpoint change in the scalar fundamental local susceptibility is
$\Delta[-q_{10}(\log_{10}\eta_N)]$ between $a/M=0$ and $a/M=1$.
Each one-factor check holds the other settings at the stated reference value.

| test | N | grid | half-width | endpoint change | sign |
|---|---:|---:|---:|---:|---|
| grid | 64 | 81 | 0.025 | 0.1613 | positive |
| grid | 64 | 121 | 0.025 | 0.1603 | positive |
| grid | 64 | 161 | 0.025 | 0.1607 | positive |
| window | 64 | 121 | 0.020 | 0.1396 | positive |
| window | 64 | 121 | 0.030 | 0.1647 | positive |
| spectral_N | 32 | 121 | 0.025 | 0.1002 | positive |
| spectral_N | 48 | 121 | 0.025 | 0.1287 | positive |

The sign is positive in 7/7 configurations; the observed changes span 0.1002--0.1647.
The sign and order of magnitude, rather than a single window-dependent decimal,
are the defensible robustness statement.
