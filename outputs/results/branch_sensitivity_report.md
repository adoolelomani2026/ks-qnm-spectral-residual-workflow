# Branch Dependence of KS Pseudospectral Sensitivity

## Branch Classes

- Robustly increasing sensitivity: 8.
- Weakly increasing sensitivity: 3.
- Non-monotonic sensitivity: 1.
- Numerically inconclusive: 2.

## Strongest Predictors

- Strongest continuous predictor of gain magnitude: `log10_condition_growth` with Pearson r=0.961.
- Strongest continuous predictor of reliable increasing class: `nearest_eigenvalue_distance_endpoint` with Pearson r=0.888.
- Best one-threshold separator in this finite sample: `overtone_load <= 0.666667` for target `reliable_increasing`, accuracy=1.000.
- Number of perfect one-threshold separators for reliable increase in this 14-branch sample: 6.

## Interpretation

The numerical condition-growth indicator tracks the magnitude of the
finite-N susceptibility gain most strongly, but it is not an independent
physical scaling law because it is itself a non-Hermitian diagnostic of
the same discretized operator.  The simplest empirical branch separator
is the dimensionless overtone load n/(ell+1/2): all reliable increasing
branches satisfy n/(ell+1/2) <= 2/3, while the nonmonotonic or
inconclusive branches have larger load.  The endpoint damping-rate shift
gives an equivalent finite-sample split, but neither threshold should be
read as a derived law.

No predictive law is established.  The conservative conclusion is that
KS pseudospectral sensitivity is branch dependent and remains an open
problem beyond the robust fundamental and first-overtone sectors.

## Class Table

| ell | n | class | gain | damping shift | condition growth | comment |
|---:|---:|---|---:|---:|---:|---|
| 0 | 0 | weakly increasing sensitivity | +0.310 | -0.020 | +0.260 | deep contour touches local window |
| 0 | 1 | non-monotonic sensitivity | +0.756 | +0.000 | +0.694 | deep contour touches local window; N=64 susceptibility not monotonic in a/M |
| 1 | 0 | robustly increasing sensitivity | +0.202 | -0.025 | +0.141 | deep contour touches local window |
| 1 | 1 | robustly increasing sensitivity | +0.455 | -0.020 | +0.402 | deep contour touches local window |
| 1 | 2 | numerically inconclusive | +0.449 | +0.296 | +0.593 | deep contour touches local window; N=64 susceptibility not monotonic in a/M; endpoint center spread across N is large |
| 2 | 0 | robustly increasing sensitivity | +0.161 | -0.026 | +0.094 | deep contour touches local window |
| 2 | 1 | robustly increasing sensitivity | +0.347 | -0.024 | +0.287 | deep contour touches local window |
| 2 | 2 | numerically inconclusive | +0.550 | -0.019 | +0.485 | deep contour touches local window; N=64 susceptibility not monotonic in a/M |
| 3 | 0 | robustly increasing sensitivity | +0.141 | -0.026 | +0.054 | trend passes finite-N scalar checks |
| 3 | 1 | robustly increasing sensitivity | +0.285 | -0.025 | +0.221 | deep contour touches local window |
| 3 | 2 | weakly increasing sensitivity | +0.428 | -0.023 | +0.368 | deep contour touches local window |
| 4 | 0 | robustly increasing sensitivity | +0.123 | -0.026 | +0.028 | trend passes finite-N scalar checks |
| 4 | 1 | robustly increasing sensitivity | +0.247 | -0.026 | +0.177 | deep contour touches local window |
| 4 | 2 | weakly increasing sensitivity | +0.371 | -0.024 | +0.302 | deep contour touches local window |
