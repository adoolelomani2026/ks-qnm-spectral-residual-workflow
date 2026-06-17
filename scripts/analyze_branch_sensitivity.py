#!/usr/bin/env python3
"""Classify KS scalar branches by finite-N pseudospectral sensitivity.

This script deliberately treats universal pseudospectral amplification as the
null hypothesis to be rejected.  It joins the existing PRL stress-test outputs
and asks which scalar branches show reliable sensitivity growth, which do not,
and whether any simple dimensionless predictor separates the classes.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Callable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT_DIR / "outputs" / "results"
FIGURES_DIR = ROOT_DIR / "outputs" / "figures"


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _as_float(value: str) -> float:
    try:
        return float(value)
    except ValueError:
        return float("nan")


def _as_bool(value: str) -> bool:
    return value.strip().lower() == "true"


def _pearson(x_values: list[float], y_values: list[float]) -> float:
    x = np.asarray(x_values, dtype=float)
    y = np.asarray(y_values, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(mask) < 3:
        return float("nan")
    return float(np.corrcoef(x[mask], y[mask])[0, 1])


def classify_branch(row: dict[str, object]) -> str:
    """Assign the conservative branch-sensitivity class."""

    if row["reliability"] == "exploratory" or row["leaver_failed"] or row["max_center_spread_endpoint"] > 5.0e-2:
        return "numerically inconclusive"
    if not row["monotonic_q10_n64"]:
        return "non-monotonic sensitivity"
    if row["reliability"] == "usable" and row["positive_gain_all_tested_n"] and row["q10_gain_n64"] > 0.0:
        return "robustly increasing sensitivity"
    if row["positive_gain_all_tested_n"] and row["q10_gain_n64"] > 0.0:
        return "weakly increasing sensitivity"
    return "numerically inconclusive"


def build_branch_rows() -> list[dict[str, object]]:
    verdicts = _read_csv(RESULTS_DIR / "prl_instability_branch_verdicts.csv")
    scan_rows = _read_csv(RESULTS_DIR / "prl_instability_scan_summary.csv")
    barrier_rows = _read_csv(RESULTS_DIR / "prl_instability_barrier_metrics.csv")
    leaver_rows = _read_csv(RESULTS_DIR / "prl_instability_endpoint_leaver_checks.csv")

    output: list[dict[str, object]] = []
    for verdict in verdicts:
        ell = int(verdict["ell"])
        overtone = int(verdict["overtone"])
        local_scan = [
            row
            for row in scan_rows
            if int(row["ell"]) == ell and int(row["overtone"]) == overtone and int(row["spectral_N"]) == 64
        ]
        baseline = next(row for row in local_scan if float(row["a_over_M"]) == 0.0)
        endpoint = next(row for row in local_scan if float(row["a_over_M"]) == 1.0)

        base_real = _as_float(baseline["center_real"])
        base_damping = -_as_float(baseline["center_imag"])
        endpoint_real = _as_float(endpoint["center_real"])
        endpoint_damping = -_as_float(endpoint["center_imag"])
        base_quality = base_real / (2.0 * base_damping)
        endpoint_quality = endpoint_real / (2.0 * endpoint_damping)

        local_barrier = [row for row in barrier_rows if int(row["ell"]) == ell]
        barrier_base = next(row for row in local_barrier if float(row["a"]) == 0.0)
        barrier_endpoint = next(row for row in local_barrier if float(row["a"]) == 1.0)

        base_height = _as_float(barrier_base["peak_height"])
        endpoint_height = _as_float(barrier_endpoint["peak_height"])
        base_width = _as_float(barrier_base["rstar_width_halfmax"])
        endpoint_width = _as_float(barrier_endpoint["rstar_width_halfmax"])
        base_curvature = abs(_as_float(barrier_base["curvature_rstar_at_peak"]))
        endpoint_curvature = abs(_as_float(barrier_endpoint["curvature_rstar_at_peak"]))

        local_leaver = [
            row for row in leaver_rows if int(row["ell"]) == ell and int(row["overtone"]) == overtone
        ]
        leaver_failed = any("failed" in row["status"] for row in local_leaver)
        leaver_differences = [
            _as_float(row["relative_difference"])
            for row in local_leaver
            if np.isfinite(_as_float(row["relative_difference"]))
        ]
        endpoint_status = next(
            (row["status"] for row in local_leaver if float(row["a_over_M"]) == 1.0),
            "missing",
        )

        base_condition = _as_float(baseline["eigen_condition_indicator"])
        endpoint_condition = _as_float(endpoint["eigen_condition_indicator"])
        base_distance = _as_float(baseline["nearest_eigenvalue_distance"])
        endpoint_distance = _as_float(endpoint["nearest_eigenvalue_distance"])

        row: dict[str, object] = {
            "ell": ell,
            "overtone": overtone,
            "q10_gain_n64": _as_float(verdict["q10_gain_n64"]),
            "q50_gain_n64": _as_float(verdict["q50_gain_n64"]),
            "monotonic_q10_n64": _as_bool(verdict["monotonic_q10_n64"]),
            "positive_gain_all_tested_n": _as_bool(verdict["positive_gain_all_tested_n"]),
            "max_center_spread_endpoint": _as_float(verdict["max_center_spread_endpoint"]),
            "reliability": verdict["reliability"],
            "real_frequency_softening": -_as_float(verdict["real_frequency_shift_n64"]),
            "damping_shift": _as_float(verdict["damping_shift_n64"]),
            "quality_factor_shift": endpoint_quality / base_quality - 1.0,
            "quality_factor_baseline": base_quality,
            "potential_height_softening": -(endpoint_height / base_height - 1.0),
            "potential_width_growth": endpoint_width / base_width - 1.0,
            "potential_curvature_softening": -(endpoint_curvature / base_curvature - 1.0),
            "log10_condition_baseline": math.log10(base_condition),
            "log10_condition_endpoint": math.log10(endpoint_condition),
            "log10_condition_growth": math.log10(endpoint_condition / base_condition),
            "nearest_eigenvalue_distance_endpoint": endpoint_distance,
            "nearest_eigenvalue_distance_growth": endpoint_distance / base_distance - 1.0,
            "overtone_load": overtone / (ell + 0.5),
            "branch_ratio": (overtone + 1.0) / (ell + 1.0),
            "leaver_endpoint_status": endpoint_status,
            "leaver_failed": leaver_failed,
            "max_endpoint_leaver_relative_difference": max(leaver_differences) if leaver_differences else float("nan"),
            "comment": verdict["comment"],
        }
        row["sensitivity_class"] = classify_branch(row)
        output.append(row)
    return output


def compute_correlations(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    predictors = [
        "ell",
        "overtone",
        "overtone_load",
        "branch_ratio",
        "real_frequency_softening",
        "damping_shift",
        "quality_factor_shift",
        "quality_factor_baseline",
        "potential_height_softening",
        "potential_width_growth",
        "potential_curvature_softening",
        "log10_condition_endpoint",
        "log10_condition_growth",
        "nearest_eigenvalue_distance_endpoint",
        "nearest_eigenvalue_distance_growth",
    ]
    q10 = [float(row["q10_gain_n64"]) for row in rows]
    robust = [1.0 if row["sensitivity_class"] == "robustly increasing sensitivity" else 0.0 for row in rows]
    reliable_increase = [
        1.0
        if row["sensitivity_class"] in {"robustly increasing sensitivity", "weakly increasing sensitivity"}
        else 0.0
        for row in rows
    ]

    output = []
    for predictor in predictors:
        values = [float(row[predictor]) for row in rows]
        output.append(
            {
                "predictor": predictor,
                "pearson_r_with_q10_gain": _pearson(values, q10),
                "pearson_r_with_robust_class": _pearson(values, robust),
                "pearson_r_with_reliable_increase_class": _pearson(values, reliable_increase),
                "sample_count": len(rows),
            }
        )
    return output


def threshold_search(
    rows: list[dict[str, object]],
    target_name: str,
    target: Callable[[dict[str, object]], bool],
) -> list[dict[str, object]]:
    predictors = [
        "ell",
        "overtone",
        "overtone_load",
        "branch_ratio",
        "real_frequency_softening",
        "damping_shift",
        "quality_factor_shift",
        "quality_factor_baseline",
        "potential_height_softening",
        "potential_width_growth",
        "potential_curvature_softening",
        "log10_condition_endpoint",
        "log10_condition_growth",
        "nearest_eigenvalue_distance_endpoint",
    ]
    truth = np.asarray([target(row) for row in rows], dtype=bool)
    output: list[dict[str, object]] = []
    for predictor in predictors:
        values = np.asarray([float(row[predictor]) for row in rows], dtype=float)
        finite_values = sorted(set(values[np.isfinite(values)]))
        best: dict[str, object] | None = None
        for threshold in finite_values:
            for direction in ("<=", ">="):
                prediction = values <= threshold if direction == "<=" else values >= threshold
                accuracy = float(np.mean(prediction == truth))
                tp = int(np.count_nonzero(prediction & truth))
                fp = int(np.count_nonzero(prediction & ~truth))
                tn = int(np.count_nonzero(~prediction & ~truth))
                fn = int(np.count_nonzero(~prediction & truth))
                candidate = {
                    "target": target_name,
                    "predictor": predictor,
                    "direction": direction,
                    "threshold": threshold,
                    "accuracy": accuracy,
                    "true_positive": tp,
                    "false_positive": fp,
                    "true_negative": tn,
                    "false_negative": fn,
                }
                if best is None or accuracy > float(best["accuracy"]):
                    best = candidate
        if best is not None:
            output.append(best)
    return sorted(output, key=lambda row: float(row["accuracy"]), reverse=True)


def plot_predictors(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    colors = {
        "robustly increasing sensitivity": "#2b6cb0",
        "weakly increasing sensitivity": "#38a169",
        "non-monotonic sensitivity": "#d69e2e",
        "numerically inconclusive": "#c53030",
    }
    markers = {0: "o", 1: "s", 2: "^"}
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.35), constrained_layout=True)
    for row in rows:
        color = colors[str(row["sensitivity_class"])]
        marker = markers[int(row["overtone"])]
        axes[0].scatter(
            float(row["log10_condition_growth"]),
            float(row["q10_gain_n64"]),
            color=color,
            marker=marker,
            s=55,
            alpha=0.85,
        )
        axes[1].scatter(
            100.0 * float(row["damping_shift"]),
            float(row["q10_gain_n64"]),
            color=color,
            marker=marker,
            s=55,
            alpha=0.85,
        )
    axes[0].set_xlabel(r"$\log_{10}[\kappa(a/M=1)/\kappa(0)]$")
    axes[0].set_ylabel(r"$\Delta[-Q_{10}(\log_{10}\eta_N)]$")
    axes[1].set_xlabel(r"damping-rate shift $\Delta[-\mathrm{Im}\,\omega]$ (%)")
    axes[1].set_ylabel(r"$\Delta[-Q_{10}(\log_{10}\eta_N)]$")
    for axis in axes:
        axis.axhline(0.0, color="0.2", linewidth=0.8, alpha=0.35)
        axis.grid(alpha=0.25)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def write_report(
    path: Path,
    rows: list[dict[str, object]],
    correlations: list[dict[str, object]],
    threshold_rows: list[dict[str, object]],
) -> None:
    counts = Counter(str(row["sensitivity_class"]) for row in rows)
    best_gain = max(correlations, key=lambda row: abs(float(row["pearson_r_with_q10_gain"])))
    best_class = max(correlations, key=lambda row: abs(float(row["pearson_r_with_reliable_increase_class"])))
    best_threshold = threshold_rows[0]
    perfect_thresholds = [
        row
        for row in threshold_rows
        if row["target"] == "reliable_increasing" and abs(float(row["accuracy"]) - 1.0) < 1.0e-12
    ]

    lines = [
        "# Branch Dependence of KS Pseudospectral Sensitivity",
        "",
        "## Branch Classes",
        "",
        f"- Robustly increasing sensitivity: {counts['robustly increasing sensitivity']}.",
        f"- Weakly increasing sensitivity: {counts['weakly increasing sensitivity']}.",
        f"- Non-monotonic sensitivity: {counts['non-monotonic sensitivity']}.",
        f"- Numerically inconclusive: {counts['numerically inconclusive']}.",
        "",
        "## Strongest Predictors",
        "",
        f"- Strongest continuous predictor of gain magnitude: `{best_gain['predictor']}` "
        f"with Pearson r={float(best_gain['pearson_r_with_q10_gain']):.3f}.",
        f"- Strongest continuous predictor of reliable increasing class: `{best_class['predictor']}` "
        f"with Pearson r={float(best_class['pearson_r_with_reliable_increase_class']):.3f}.",
        f"- Best one-threshold separator in this finite sample: `{best_threshold['predictor']} "
        f"{best_threshold['direction']} {float(best_threshold['threshold']):.6g}` for target "
        f"`{best_threshold['target']}`, accuracy={float(best_threshold['accuracy']):.3f}.",
        f"- Number of perfect one-threshold separators for reliable increase in this 14-branch sample: {len(perfect_thresholds)}.",
        "",
        "## Interpretation",
        "",
        "The numerical condition-growth indicator tracks the magnitude of the",
        "finite-N susceptibility gain most strongly, but it is not an independent",
        "physical scaling law because it is itself a non-Hermitian diagnostic of",
        "the same discretized operator.  The simplest empirical branch separator",
        "is the dimensionless overtone load n/(ell+1/2): all reliable increasing",
        "branches satisfy n/(ell+1/2) <= 2/3, while the nonmonotonic or",
        "inconclusive branches have larger load.  The endpoint damping-rate shift",
        "gives an equivalent finite-sample split, but neither threshold should be",
        "read as a derived law.",
        "",
        "No predictive law is established.  The conservative conclusion is that",
        "KS pseudospectral sensitivity is branch dependent and remains an open",
        "problem beyond the robust fundamental and first-overtone sectors.",
        "",
        "## Class Table",
        "",
        "| ell | n | class | gain | damping shift | condition growth | comment |",
        "|---:|---:|---|---:|---:|---:|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['ell']} | {row['overtone']} | {row['sensitivity_class']} | "
            f"{float(row['q10_gain_n64']):+.3f} | {float(row['damping_shift']):+.3f} | "
            f"{float(row['log10_condition_growth']):+.3f} | {row['comment']} |"
        )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    global RESULTS_DIR, FIGURES_DIR

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--figures-dir", type=Path, default=FIGURES_DIR)
    args = parser.parse_args()

    RESULTS_DIR = args.results_dir
    FIGURES_DIR = args.figures_dir

    rows = build_branch_rows()
    correlations = compute_correlations(rows)
    thresholds = []
    thresholds.extend(
        threshold_search(
            rows,
            "robustly_increasing",
            lambda row: row["sensitivity_class"] == "robustly increasing sensitivity",
        )
    )
    thresholds.extend(
        threshold_search(
            rows,
            "reliable_increasing",
            lambda row: row["sensitivity_class"]
            in {"robustly increasing sensitivity", "weakly increasing sensitivity"},
        )
    )
    thresholds = sorted(thresholds, key=lambda row: float(row["accuracy"]), reverse=True)

    _write_csv(args.results_dir / "branch_sensitivity_classes.csv", rows)
    _write_csv(args.results_dir / "branch_sensitivity_correlations.csv", correlations)
    _write_csv(args.results_dir / "branch_sensitivity_predictor_search.csv", thresholds)
    plot_predictors(args.figures_dir / "branch_sensitivity_predictors.png", rows)
    write_report(args.results_dir / "branch_sensitivity_report.md", rows, correlations, thresholds)

    counts = Counter(str(row["sensitivity_class"]) for row in rows)
    print("Branch sensitivity classification complete.")
    for name, count in counts.items():
        print(f"{name}: {count}")
    print(f"Wrote {args.results_dir / 'branch_sensitivity_report.md'}")


if __name__ == "__main__":
    sys.exit(main())
