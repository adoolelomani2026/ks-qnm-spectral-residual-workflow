# PRL Discovery-Hunt Report

## Verdict

The project is stronger after the model-zoo stress test, but the result is not
PRL-level. It is best framed as a PRD/CQG-quality negative universality result:
the KS deformation broadens the finite-dimensional scalar residual
pseudospectrum on usable branches, but the same diagnostic narrows for two
regular-black-hole comparators.

## Strongest Defensible Claim

KS scalar branches show simultaneous frequency softening and increased
finite-N pseudospectral susceptibility, whereas Hayward and Bardeen regular
black holes mostly harden the scalar spectrum and reduce the same residual
susceptibility; pseudospectral response is therefore deformation-specific, not
a generic consequence of black-hole regularization or quantum-inspired
correction.

## Evidence Against PRL Universality

- Models tested: KS, Hayward, Bardeen.
- Branch verdicts tested: 30 scalar model-branch combinations.
- Universality-supporting verdicts: 8/30.
- KS support: 8/10.
- Hayward support: 0/10.
- Bardeen support: 0/10.
- Hayward and Bardeen show negative endpoint susceptibility gains on every
  tested branch.

This falsifies the broad claim that quantum-corrected or regular black holes
generically amplify non-Hermitian ringdown sensitivity.

## Mechanism Evidence

The most predictive quantities in the three-model scan are damping shift,
horizon expansion, eikonal Lyapunov change, tortoise-width change, and
photon-sphere frequency softening. The strongest correlations are useful
diagnostics but do not yet constitute a universal scaling law.

The qualitative mechanism is:

- KS lowers the scalar barrier and broadens the local residual susceptibility.
- Hayward and Bardeen raise the scalar barrier for most ell >= 1 branches and
  narrow the same diagnostic.
- Barrier width is anticorrelated with susceptibility gain in the current
  model zoo.

This is a deformation-response taxonomy, not a discovery-level principle.

## Exceptional-Point Search

No eigenvalue/eigenvector coalescence criterion is met. The mode-pair scan
finds no near-coalescent candidates, and Petermann-like condition indicators do
not reveal a non-Hermitian critical point.

## Journal Recommendation

Recommended target: PRD or CQG.

EPJ C is plausible if the manuscript remains compact and focused. EPJ Plus
remains safe. PRL is not recommended unless a future scan finds a simple,
comparator-independent scaling law or genuine non-Hermitian critical behavior.

## Files Generated

- `outputs/results/universality_model_scan.csv`
- `outputs/results/universality_branch_verdicts.csv`
- `outputs/results/universality_barrier_metrics.csv`
- `outputs/results/universality_correlations.csv`
- `outputs/results/universality_assessment.md`
- `outputs/figures/universality_endpoint_gain_heatmap.png`
- `outputs/figures/universality_softening_vs_sensitivity.png`
- `outputs/figures/universality_barrier_correlation.png`
- `outputs/figures/universality_width_correlation.png`
