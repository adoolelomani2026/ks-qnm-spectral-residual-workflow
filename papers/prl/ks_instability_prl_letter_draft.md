# PRL-Style Letter Draft: Not Yet Submission Ready

## Title

Does Quantum Deformation Amplify Black-Hole Ringdown Spectral Instability?

## 600-word abstract/introduction draft

Quasinormal modes are often treated as the spectral fingerprints of compact
objects, but black-hole ringdown spectra are eigenvalues of non-Hermitian
operators and can be highly sensitive to perturbations.  This point has become
central in recent pseudospectrum studies of black-hole quasinormal modes: a
small movement of the operator can produce a large movement of the spectrum,
so a frequency catalogue alone is not a complete stability diagnostic.  A
natural question is whether quantum-corrected black-hole geometries merely
shift the ringdown frequencies, or whether they also modify the non-Hermitian
sensitivity of the QNM operator itself.

Here we test that question for the Kazakov--Solodukhin deformation of the
Schwarzschild black hole.  The proposed Letter-level claim would be simple:
quantum deformation softens the ringdown spectrum while amplifying
pseudospectral instability.  We stress-test this claim using a Chebyshev
spectral residual discretization of the scalar perturbation equation.  For
each branch we compute local grids of the relative smallest singular value
eta_N(omega)=sigma_min(P_N)/sigma_max(P_N), where P_N is the finite Chebyshev
operator obtained after compactification and quasinormal-mode asymptotic
factorization.  We scan scalar multipoles ell=0,...,4, overtones n=0,1,2 where
the branch is numerically trackable, deformations a/M=0,...,1, and Chebyshev
sizes N=32,48,64, with N=96 checks for fundamentals.

The result is not yet a Physical Review Letters result.  Several branches do
show the expected paired behavior: the real frequency decreases and the local
finite-N pseudospectral susceptibility, measured by -Q10(log10 eta_N),
increases.  The scalar fundamental branches are the cleanest part of the
evidence.  However, the effect is not yet established as a universal law of
quantum-deformed ringdown.  Low multipoles and second overtones require
additional branch discipline, contour areas depend on finite-N normalization
and local-window choices, and no continuum operator theorem or independent
quantum-black-hole comparator is yet available.  The data therefore support a
strong PRD/CQG-style statement about finite-dimensional KS scalar
pseudospectra, but not a PRL-level claim of universal quantum amplification of
black-hole spectral instability.

## Four-page Letter skeleton

### 1. Opening result

Black-hole QNM spectra are non-Hermitian objects.  A quantum deformation may
change not only the QNM frequencies but also the operator's susceptibility to
perturbation.  We test whether the KS deformation produces a universal
increase in pseudospectral sensitivity.

### 2. Method in one paragraph

For each scalar branch we form P_N(omega) from the compactified Chebyshev
collocation equation and evaluate eta_N=sigma_min(P_N)/sigma_max(P_N) on local
complex-frequency grids centered on tracked QNM roots.  Branches are compared
across a/M and N using quantiles, fixed-threshold area fractions, absolute
singular-value checks, eigenvalue-condition indicators, and mode-pair
diagnostics.

### 3. Main figure

Use `outputs/figures/prl_instability_central_heatmap.png`.  This figure should
be the honesty test: if the heatmap does not show a clean, robust, universal
positive effect across reliable branches, the Letter claim fails.

### 4. Robustness figure

Use `outputs/figures/prl_instability_robustness.png`.  It shows whether the
fundamental-branch trend survives changes in Chebyshev size through N=96.

### 5. Physical interpretation

The KS deformation lowers the scalar barrier scale and shifts the ringdown
frequencies.  Barrier diagnostics are consistent with softening, but they do
not yet provide a derivation of pseudospectral amplification.  No exceptional
point is established by the current mode-pair diagnostics.

### 6. Conclusion

The current defensible conclusion is narrower than the PRL claim: KS scalar
ringdown softening is accompanied on several branches by increased local
finite-N pseudospectral susceptibility.  This is important and publishable,
but not yet universal enough, mechanistic enough, or model-comparative enough
for PRL.
