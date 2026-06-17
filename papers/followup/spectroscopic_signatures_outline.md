# Spectroscopic Signatures of Quantum Deformations in Kazakov-Solodukhin Black Holes

## Working Thesis

Kazakov-Solodukhin deformation produces a coherent softening of the ringdown spectrum: validated scalar QNM branches shift to lower oscillation frequency and lower damping magnitude as `a/M` increases. The scalar `ell=2` fundamental mode is the safest leading signature.

## Core Result

For the scalar `ell=2` fundamental branch at `a/M=1`:

- `Delta Re(M omega) / Re(M omega_0) = -7.29%`
- `Delta[-Im(M omega)] / [-Im(M omega_0)] = -2.59%`
- `Delta Q / Q_0 = -4.83%`
- `omega0/omega1` complex-ratio shift = `1.92%`
- `omega0/omega2` complex-ratio shift = `3.56%`

Interpretation: the mode rings at a lower tone and damps more slowly in absolute time, but the frequency shift is larger than the damping-rate shift, so the dimensionless quality factor decreases.

## Proposed Abstract Shape

1. Black-hole spectroscopy constrains deviations from Schwarzschild/Kerr through QNM frequencies and damping times.
2. Use the existing Chebyshev-Leaver catalogue to isolate how KS deformation changes scalar QNM spectra.
3. Show monotonic softening across sampled deformation values.
4. Emphasize the scalar sector as the clean physical result.
5. Present axial gravitational results only as phenomenological diagnostics.

## Paper Structure

1. Introduction: why quantum-deformed ringdown signatures matter.
2. KS background and scalar perturbation model.
3. Validated catalogue data set and uncertainty conventions.
4. Fundamental-mode deformation signatures.
5. Multipole and overtone dependence.
6. Mode-ratio diagnostics and mass-scaling degeneracy.
7. Comparison with existing KS QNM literature and parameter conventions.
8. Phenomenological axial-gravitational comparison.
9. Observational interpretation without detectability claims.
10. Limitations and next theory steps.

## Figures To Build

- Scalar-only heatmap of endpoint fractional shifts.
- `ell=2` scalar branch plot: real shift, damping shift, quality-factor shift.
- Mode-ratio plots, starting from the generated scalar `ell=2`
  `omega0/omega1`, `omega0/omega2`, and `Re(omega)/[-Im(omega)]`
  diagnostics.
- Overtone-sensitivity plot with validation-error overlays.
- Optional axial/scalar comparison figure clearly labeled phenomenological.

## Tables To Include

- Publication-safe scalar endpoint shifts.
- Sensitivity ranking by branch.
- Validation and uncertainty table for all branches used in claims.
- Literature comparison table: previous KS WKB/time-domain and
  Frobenius/continued-fraction results versus the fixed-`M`
  Chebyshev-Leaver catalogue.
- Caveat table separating scalar, axial phenomenological, and exploratory overtone claims.

## Claims That Are Safe Now

- Scalar fundamental KS QNMs soften monotonically on the sampled grid.
- The scalar `ell=2` fundamental shift is percent-level at `a/M=1`.
- Real-frequency shifts dominate damping-rate shifts.
- Dimensionless scalar mode ratios shift coherently and help reduce simple
  mass-rescaling ambiguity.
- Second overtones are more deformation-sensitive but less numerically secure.

## Claims To Avoid For Now

- Do not claim a full KS gravitational perturbation result.
- Do not claim detectability with current detectors.
- Do not claim high-`N` overtone stability without independent branch validation.
- Keep the residual workflow framed as a numerical diagnostic, not as a separate platform claim.

## Minimum New Work Before Drafting

- Add scalar-only analysis plots that do not rely on the phenomenological axial sector.
- Expand mode-ratio diagnostics beyond the current scalar `ell=2` generated table.
- Use `docs/LITERATURE_NORMALIZATION_PROTOCOL.md` and `src/qnm/normalization.py`
  to translate existing KS QNM literature conventions onto the same fixed-mass
  comparison grid, or state clearly where only qualitative comparison is
  possible.
- Transcribe the relevant published scalar rows into
  `data/literature/ks_qnm_literature_template.csv`, enable them, and run
  `scripts/prepare_literature_comparison.py` to generate the first
  normalization-matched comparison table.
- Add a short uncertainty convention for second overtones.
- Optionally run a denser `a/M` grid for scalar fundamentals only.
