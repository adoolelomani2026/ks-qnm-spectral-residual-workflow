# Axial inverse-Cowling audit

- Local source commit: `21e14ca6c1eee13bb265911c2771643b94b8ea4f` (the run metadata separately records whether the tree was dirty).
- Public Batic--Dutykh--Sukaiti repository commit inspected: `f53435ebdb8d1124d13bee75fa54972ac1613a03`.
- Public repository first commit and date: `dcdc3f6c07584b17c83662eae3b8ac4a424dae30 2026-06-07T06:49:31+04:00`.
- The public files use the same mass-normalised deformation `a/M` and frequency `M omega`.
- Largest external relative difference over the 18 overlapping modes: `1.618e-05`.
- Smallest sampled exterior potential: `6.444e-11`.
- The growing-mode search was performed on the unfiltered finite generalized spectrum before the ordinary damped-mode window was applied.
- Search domain: `Im(omega) > 1e-8` and `|omega| < 5`, at `N=32,48,64`.
- A candidate had to persist within `1e-4` at all three sizes; any survivor was then locally residual-refined within radius `0.02` at `N=64` and required polynomial backward error below `1e-8`.
- Three-resolution candidates before residual acceptance: `0`.
- Accepted growing-mode survivors: `0`.

Potential positivity is a sufficient mode-stability diagnostic for the source-free one-dimensional problem; the finite raw-spectrum search is an additional numerical check, not a proof for a different source closure.
