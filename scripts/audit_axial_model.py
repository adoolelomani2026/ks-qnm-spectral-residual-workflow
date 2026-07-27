#!/usr/bin/env python3
"""Audit the inverse-Cowling axial potential and its public-data overlap."""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar

from qnm.common import A_VALUES, df_ks, f_ks, horizon_radius, regge_wheeler_potential
from qnm.spectral import (
    BACKWARD_ERROR_ACCEPTANCE_THRESHOLD,
    build_spectral_problem,
    generalized_eigenvalues,
    minimize_residual,
    polynomial_backward_error,
)


PUBLIC = ROOT / "papers" / "revision" / "references" / "KS-quantum-public"
PUBLIC_URL = "https://github.com/dutykh/KS-quantum.git"
PUBLIC_COMMIT = "f53435ebdb8d1124d13bee75fa54972ac1613a03"
RESULTS = ROOT / "outputs" / "results"
FIGURES = ROOT / "outputs" / "figures"


def ensure_public_checkout() -> None:
    """Materialize and pin the external comparison data when absent."""

    if not (PUBLIC / ".git").is_dir():
        PUBLIC.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(["git", "clone", PUBLIC_URL, str(PUBLIC)], check=True)
    subprocess.run(["git", "checkout", "--detach", PUBLIC_COMMIT], cwd=PUBLIC, check=True)
    actual = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=PUBLIC, capture_output=True, text=True, check=True
    ).stdout.strip()
    if actual != PUBLIC_COMMIT:
        raise RuntimeError(f"External checkout mismatch: expected {PUBLIC_COMMIT}, found {actual}")


def lapse_substitution_potential(r: np.ndarray, ell: int, a: float) -> np.ndarray:
    f = f_ks(r, a)
    return f * (ell * (ell + 1.0) / r**2 - 6.0 / r**3)


def plot_potentials() -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    colors = {0.0: "#333333", 0.5: "#2878b5", 1.0: "#d95319"}
    for a in (0.0, 0.5, 1.0):
        rh = horizon_radius(a)
        x = np.geomspace(1.0 + 1.0e-7, 15.0, 5000)
        r = rh * x
        new = regge_wheeler_potential(r, 2, a)
        lapse = lapse_substitution_potential(r, 2, a)
        delta = new - lapse
        peak_index = int(np.argmax(new))
        peak = float(new[peak_index])
        mask = new >= 1.0e-3 * peak
        max_peak_scaled = float(np.max(np.abs(delta)) / peak)
        max_pointwise = float(np.max(np.abs(delta[mask] / new[mask])))
        rows.append(
            {
                "a_over_M": a,
                "r_h_over_M": rh,
                "r_peak_over_M": float(r[peak_index]),
                "V_peak_M2": peak,
                "max_abs_deltaV_over_Vpeak": max_peak_scaled,
                "max_pointwise_fraction_where_Vnew_ge_1e-3_peak": max_pointwise,
            }
        )
        label = rf"$a/M={a:g}$"
        axes[0].plot(x, new, color=colors[a], label=label)
        axes[0].plot(x, lapse, color=colors[a], linestyle="--", alpha=0.65)
        axes[0].plot(x[peak_index], peak, "o", color=colors[a], ms=4)
        axes[1].plot(x, delta / peak if peak else delta, color=colors[a], label=label)

    axes[0].set(xlabel=r"$r/r_h$", ylabel=r"$M^2 V$", xlim=(1, 8))
    axes[0].set_title("Solid: inverse-Cowling; dashed: lapse substitution")
    axes[1].set(
        xlabel=r"$r/r_h$",
        ylabel=r"$(V_{\rm ax}-V_{\rm lapse})/V_{\rm ax}^{\rm peak}$",
        xlim=(1, 8),
    )
    axes[1].set_title("Comparison-potential difference")
    for ax in axes:
        ax.axvline(1.0, color="k", lw=0.8, alpha=0.35)
        ax.grid(alpha=0.2)
    axes[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURES / "axial_potential_comparison.png", dpi=220)
    plt.close(fig)
    return rows


def positivity_rows() -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for ell in (2, 3, 4):
        for a in A_VALUES:
            rh = horizon_radius(a)

            def potential_logx(logx: float) -> float:
                r = rh * (1.0 + np.exp(logx))
                return float(regge_wheeler_potential(np.array([r]), ell, a)[0])

            # Exclude the two limiting zeros and search the exterior over
            # r-r_h in [1e-10 r_h, 1e4 r_h].
            result = minimize_scalar(potential_logx, bounds=(-23.02585, 9.21034), method="bounded")
            grid_x = np.geomspace(1.0e-10, 1.0e4, 20000)
            grid_r = rh * (1.0 + grid_x)
            grid_v = regge_wheeler_potential(grid_r, ell, a)
            rows.append(
                {
                    "ell": ell,
                    "a_over_M": a,
                    "minimum_sampled_V_M2": float(min(result.fun, np.min(grid_v))),
                    "minimum_bracket_sampled_M2": float(
                        np.min(
                            ell * (ell + 1.0) / grid_r**2
                            + 2.0 * (f_ks(grid_r, a) - 1.0) / grid_r**2
                            - df_ks(grid_r, a) / grid_r
                        )
                    ),
                }
            )
    return rows


def growing_candidate_audit() -> list[dict[str, object]]:
    spectra: dict[tuple[int, float, int], np.ndarray] = {}
    problems = {}
    for ell in (2, 3, 4):
        for a in A_VALUES:
            for n in (32, 48, 64):
                problem = build_spectral_problem(a, n, ell=ell, perturbation_type="gravitational")
                problems[(ell, a, n)] = problem
                values = generalized_eigenvalues(problem)
                spectra[(ell, a, n)] = np.array(
                    [z for z in values if np.isfinite(z) and z.imag > 1.0e-8 and abs(z) < 5.0],
                    dtype=complex,
                )

    rows: list[dict[str, object]] = []
    for ell in (2, 3, 4):
        for a in A_VALUES:
            base = spectra[(ell, a, 32)]
            persistent = []
            accepted = []
            for z in base:
                if all(
                    len(spectra[(ell, a, n)])
                    and np.min(np.abs(spectra[(ell, a, n)] - z)) < 1.0e-4
                    for n in (48, 64)
                ):
                    persistent.append(z)
                    n64_values = spectra[(ell, a, 64)]
                    n64_candidate = n64_values[int(np.argmin(np.abs(n64_values - z)))]
                    refined, _ = minimize_residual(
                        problems[(ell, a, 64)], n64_candidate, radius=2.0e-2
                    )
                    if (
                        refined.imag > 1.0e-8
                        and abs(refined) < 5.0
                        and polynomial_backward_error(problems[(ell, a, 64)], refined)
                        < BACKWARD_ERROR_ACCEPTANCE_THRESHOLD
                    ):
                        accepted.append(refined)
            rows.append(
                {
                    "ell": ell,
                    "a_over_M": a,
                    "raw_positive_imag_counts_N32_N48_N64": ";".join(
                        str(len(spectra[(ell, a, n)])) for n in (32, 48, 64)
                    ),
                    "positive_imaginary_floor": 1.0e-8,
                    "modulus_ceiling": 5.0,
                    "three_resolution_matching_tolerance": 1.0e-4,
                    "residual_refinement_radius": 2.0e-2,
                    "backward_error_acceptance_threshold": BACKWARD_ERROR_ACCEPTANCE_THRESHOLD,
                    "persistent_before_residual": len(persistent),
                    "accepted_growing_mode_survivors": len(accepted),
                    "accepted_survivor_frequencies": ";".join(
                        f"{z.real:.12g}{z.imag:+.12g}i" for z in accepted
                    ),
                }
            )
    return rows


def public_modes(path: Path) -> list[complex]:
    modes: list[complex] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) < 3 or not fields[0].isdigit():
            continue
        real, imag = float(fields[1]), float(fields[2])
        if real > 1.0e-10 and imag < 0.0:
            modes.append(complex(real, imag))
    return sorted(modes, key=lambda z: -z.imag)


def public_comparison() -> list[dict[str, object]]:
    catalogue = list(csv.DictReader((RESULTS / "qnm_catalogue.csv").open(newline="")))
    rows: list[dict[str, object]] = []
    for a, tag in ((0.5, "a0_5"), (1.0, "a1_0")):
        for ell in (2, 3, 4):
            external = public_modes(PUBLIC / "reports" / "s2-axial" / f"s2_L{ell}_{tag}.txt")
            for overtone in range(3):
                ours_row = next(
                    row for row in catalogue
                    if row["perturbation_type"] == "gravitational"
                    and int(row["ell"]) == ell
                    and int(row["overtone"]) == overtone
                    and float(row["a_over_M"]) == a
                )
                ours = complex(float(ours_row["leaver_real"]), float(ours_row["leaver_imag"]))
                theirs = external[overtone]
                rows.append(
                    {
                        "a_over_M": a,
                        "ell": ell,
                        "overtone": overtone,
                        "our_leaver_real": ours.real,
                        "our_leaver_imag": ours.imag,
                        "public_chebyshev_real": theirs.real,
                        "public_chebyshev_imag": theirs.imag,
                        "relative_difference": abs(ours - theirs) / abs(theirs),
                    }
                )
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    ensure_public_checkout()
    potential = plot_potentials()
    positivity = positivity_rows()
    growing = growing_candidate_audit()
    external = public_comparison()
    write_csv(RESULTS / "axial_potential_comparison.csv", potential)
    write_csv(RESULTS / "axial_potential_positivity.csv", positivity)
    write_csv(RESULTS / "axial_growing_mode_audit.csv", growing)
    write_csv(RESULTS / "axial_external_comparison.csv", external)

    public_hash = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=PUBLIC, capture_output=True, text=True, check=True
    ).stdout.strip()
    root_hash = subprocess.run(
        ["git", "rev-list", "--max-parents=0", "HEAD"],
        cwd=PUBLIC,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    public_first = subprocess.run(
        ["git", "show", "-s", "--format=%H %aI", root_hash],
        cwd=PUBLIC,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    local_hash = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout.strip()
    report = [
        "# Axial inverse-Cowling audit",
        "",
        f"- Local source commit: `{local_hash}` (the run metadata separately records whether the tree was dirty).",
        f"- Public Batic--Dutykh--Sukaiti repository commit inspected: `{public_hash}`.",
        f"- Public repository first commit and date: `{public_first}`.",
        "- The public files use the same mass-normalised deformation `a/M` and frequency `M omega`.",
        f"- Largest external relative difference over the 18 overlapping modes: `{max(float(r['relative_difference']) for r in external):.3e}`.",
        f"- Smallest sampled exterior potential: `{min(float(r['minimum_sampled_V_M2']) for r in positivity):.3e}`.",
        "- The growing-mode search was performed on the unfiltered finite generalized spectrum before the ordinary damped-mode window was applied.",
        "- Search domain: `Im(omega) > 1e-8` and `|omega| < 5`, at `N=32,48,64`.",
        "- A candidate had to persist within `1e-4` at all three sizes; any survivor was then locally residual-refined within radius `0.02` at `N=64` and required polynomial backward error below `1e-8`.",
        f"- Three-resolution candidates before residual acceptance: `{sum(int(r['persistent_before_residual']) for r in growing)}`.",
        f"- Accepted growing-mode survivors: `{sum(int(r['accepted_growing_mode_survivors']) for r in growing)}`.",
        "",
        "Potential positivity is a sufficient mode-stability diagnostic for the source-free one-dimensional problem; the finite raw-spectrum search is an additional numerical check, not a proof for a different source closure.",
    ]
    (RESULTS / "axial_model_audit_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
