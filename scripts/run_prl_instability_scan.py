#!/usr/bin/env python3
"""Run the PRL-level scalar KS pseudospectral-instability stress test."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from qnm.prl_instability import (
    compute_barrier_metrics,
    compute_pair_diagnostics,
    compute_window_sensitivity,
    plot_central_heatmap,
    plot_robustness,
    plot_softening_scatter,
    run_endpoint_leaver_checks,
    run_scan,
    summarize_branches,
    write_assessment_report,
    write_barrier_csv,
    write_pair_csv,
    write_scan_csv,
    write_verdict_csv,
    _write_dict_csv,
)


def write_letter_draft(output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        """# PRL-Style Letter Draft: Not Yet Submission Ready

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
""",
        encoding="utf-8",
    )


def write_supplement_outline(output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        """# Supplemental Material Outline

1. KS scalar perturbation equation and compactification.
2. Chebyshev collocation construction of P_N(omega).
3. Residual singular-value diagnostics and eta_N normalization.
4. Branch tracking protocol for ell=0,...,4 and n=0,1,2.
5. Leaver endpoint checks and failure modes.
6. Finite-N convergence and N=96 fundamental checks.
7. Absolute versus relative residual pseudospectra.
8. Window-size sensitivity of local pseudospectrum quantiles.
9. Effective-potential barrier metrics.
10. Mode-pair distances, overlaps, and exceptional-point exclusion criteria.
11. Full CSV tables and reproduction commands.

Reproduction command:

```bash
python scripts/run_prl_instability_scan.py
```
""",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=ROOT_DIR / "outputs" / "results")
    parser.add_argument("--figures-dir", type=Path, default=ROOT_DIR / "outputs" / "figures")
    parser.add_argument("--paper-dir", type=Path, default=ROOT_DIR / "papers" / "prl")
    parser.add_argument("--main-grid-size", type=int, default=31)
    parser.add_argument("--resolution-grid-size", type=int, default=21)
    args = parser.parse_args()

    args.results_dir.mkdir(parents=True, exist_ok=True)
    args.figures_dir.mkdir(parents=True, exist_ok=True)
    args.paper_dir.mkdir(parents=True, exist_ok=True)

    print("Running broad scalar pseudospectrum scan...")
    rows = run_scan(main_grid_size=args.main_grid_size, resolution_grid_size=args.resolution_grid_size)
    verdicts = summarize_branches(rows)
    barriers = compute_barrier_metrics()
    pairs = compute_pair_diagnostics(rows)
    windows = compute_window_sensitivity(rows)
    print("Running endpoint Leaver checks...")
    leaver_checks = run_endpoint_leaver_checks(rows)

    scan_csv = args.results_dir / "prl_instability_scan_summary.csv"
    verdict_csv = args.results_dir / "prl_instability_branch_verdicts.csv"
    barrier_csv = args.results_dir / "prl_instability_barrier_metrics.csv"
    pair_csv = args.results_dir / "prl_instability_mode_pair_diagnostics.csv"
    window_csv = args.results_dir / "prl_instability_window_sensitivity.csv"
    leaver_csv = args.results_dir / "prl_instability_endpoint_leaver_checks.csv"
    report = args.results_dir / "prl_instability_assessment.md"

    write_scan_csv(scan_csv, rows)
    write_verdict_csv(verdict_csv, verdicts)
    write_barrier_csv(barrier_csv, barriers)
    write_pair_csv(pair_csv, pairs)
    _write_dict_csv(window_csv, windows)
    _write_dict_csv(leaver_csv, leaver_checks)
    write_assessment_report(report, verdicts, pairs, leaver_checks)

    central = args.figures_dir / "prl_instability_central_heatmap.png"
    robustness = args.figures_dir / "prl_instability_robustness.png"
    scatter = args.figures_dir / "prl_instability_softening_scatter.png"
    plot_central_heatmap(central, verdicts)
    plot_robustness(robustness, rows)
    plot_softening_scatter(scatter, verdicts)

    letter = args.paper_dir / "ks_instability_prl_letter_draft.md"
    supplement = args.paper_dir / "supplemental_material_outline.md"
    write_letter_draft(letter)
    write_supplement_outline(supplement)

    support = [verdict for verdict in verdicts if verdict.prl_support]
    print(f"Wrote scan CSV: {scan_csv}")
    print(f"Wrote verdict CSV: {verdict_csv}")
    print(f"Wrote report: {report}")
    print(f"Wrote central figure: {central}")
    print(f"Wrote robustness figure: {robustness}")
    print(f"Wrote PRL-style draft: {letter}")
    print(f"PRL-supporting usable branches: {len(support)} / {len(verdicts)}")
    print("Verdict: not PRL-level yet; see assessment report.")


if __name__ == "__main__":
    main()
