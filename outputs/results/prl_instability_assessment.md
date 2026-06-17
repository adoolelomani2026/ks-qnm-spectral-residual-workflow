# PRL Instability Assessment

## Short Verdict: not PRL-level yet

The scan tests the proposed Letter claim that KS quantum deformation softens
scalar ringdown while universally amplifying finite-dimensional pseudospectral
sensitivity.  The result is interesting, but the evidence is not yet strong
enough for a Physical Review Letters claim.

## Strongest Supported One-Sentence Claim

In the current Chebyshev residual normalization, several scalar KS branches
show simultaneous frequency softening and increased local finite-N
pseudospectral susceptibility, but the trend is branch- and reliability-dependent
rather than a demonstrated universal instability law.

## Branch-Level Outcome

- Active oscillatory scalar branches scanned: 14.
- Usable branches by finite-N criteria: 8.
- Usable branches supporting positive monotonic gain: 8.
- Branches with non-positive N=64 endpoint gain: 0.
- Branches with nonmonotonic N=64 susceptibility in a/M: 3.
- Endpoint Leaver checks that failed to converge: 2.
- Largest endpoint Leaver relative difference among successful checks: 8.244e-02.
- Exceptional-point-like pair candidates found: 0.

## Why This Is Not A PRL Yet

- The claim is still based on finite Chebyshev matrices, not a continuum
  pseudospectrum theorem for the KS wave operator.
- The low-ell and second-overtone sectors require caution; ell=0,n=2 is not
  included as an oscillatory branch in the current selector.
- Endpoint Leaver checks expose overtone risk: some checks fail to converge,
  and the largest successful exploratory discrepancy is too large for a
  Letter-level universal claim.
- No independent quantum-corrected or regular black-hole comparator has been
  implemented, so the effect cannot yet be called universal across quantum
  black-hole models.
- No analytical mechanism has been derived; barrier metrics are diagnostic
  rather than explanatory.
- No exceptional point is supported: frequency distances and mode-shape
  overlaps do not demonstrate eigenvalue/eigenvector coalescence.

## Publishable Direction

This is strong PRD/CQG material if presented as a scalar finite-N
pseudospectrum audit of KS ringdown softening.  It should not be submitted
as a PRL unless a stronger mechanism, continuum robustness argument, and at
least one comparator metric are added.

## Branch Table

| ell | n | gain | monotonic | all-N positive | reliability | comment |
|---:|---:|---:|:---:|:---:|---|---|
| 0 | 0 | +0.310 | True | True | caution | deep contour touches local window |
| 0 | 1 | +0.756 | False | True | caution | deep contour touches local window; N=64 susceptibility not monotonic in a/M |
| 1 | 0 | +0.202 | True | True | usable | deep contour touches local window |
| 1 | 1 | +0.455 | True | True | usable | deep contour touches local window |
| 1 | 2 | +0.449 | False | True | exploratory | deep contour touches local window; N=64 susceptibility not monotonic in a/M; endpoint center spread across N is large |
| 2 | 0 | +0.161 | True | True | usable | deep contour touches local window |
| 2 | 1 | +0.347 | True | True | usable | deep contour touches local window |
| 2 | 2 | +0.550 | False | True | caution | deep contour touches local window; N=64 susceptibility not monotonic in a/M |
| 3 | 0 | +0.141 | True | True | usable | trend passes finite-N scalar checks |
| 3 | 1 | +0.285 | True | True | usable | deep contour touches local window |
| 3 | 2 | +0.428 | True | True | caution | deep contour touches local window |
| 4 | 0 | +0.123 | True | True | usable | trend passes finite-N scalar checks |
| 4 | 1 | +0.247 | True | True | usable | deep contour touches local window |
| 4 | 2 | +0.371 | True | True | caution | deep contour touches local window |
