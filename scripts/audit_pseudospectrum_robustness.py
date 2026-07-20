#!/usr/bin/env python3
"""Audit the scalar q10 pseudospectral trend against grid, window, and N."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from qnm.pseudospectrum import (  # noqa: E402
    compute_pseudospectrum_grid,
    scalar_fundamental_targets,
    summarize_grid,
)


CONFIGURATIONS = (
    ("grid", 64, 81, 0.025),
    ("grid", 64, 121, 0.025),
    ("grid", 64, 161, 0.025),
    ("window", 64, 121, 0.020),
    ("window", 64, 121, 0.030),
    ("spectral_N", 32, 121, 0.025),
    ("spectral_N", 48, 121, 0.025),
)


def main() -> None:
    catalogue = ROOT / "outputs" / "results" / "qnm_catalogue.csv"
    output = ROOT / "outputs" / "results" / "pseudospectrum_robustness.csv"
    report = ROOT / "outputs" / "results" / "pseudospectrum_robustness_report.md"
    targets = scalar_fundamental_targets(catalogue)
    endpoint_targets = {a: targets[a] for a in (0.0, 1.0)}
    rows: list[dict[str, object]] = []

    for test, spectral_n, grid_size, half_width in CONFIGURATIONS:
        summaries = []
        for a, target in endpoint_targets.items():
            grid = compute_pseudospectrum_grid(
                a=a,
                leaver_target=target,
                spectral_n=spectral_n,
                grid_size=grid_size,
                half_width=half_width,
            )
            summaries.append(summarize_grid(grid))
        summaries.sort(key=lambda item: item.a)
        gain = -summaries[1].quantiles[0.10] + summaries[0].quantiles[0.10]
        rows.append(
            {
                "test": test,
                "spectral_N": spectral_n,
                "grid_size": grid_size,
                "half_width": half_width,
                "minus_q10_a0": -summaries[0].quantiles[0.10],
                "minus_q10_a1": -summaries[1].quantiles[0.10],
                "endpoint_gain": gain,
                "positive_gain": gain > 0.0,
            }
        )
        print(test, spectral_n, grid_size, half_width, f"gain={gain:.6f}", flush=True)

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    gains = [float(row["endpoint_gain"]) for row in rows]
    lines = [
        "# Pseudospectrum robustness audit",
        "",
        "The endpoint change in the scalar fundamental local susceptibility is",
        r"$\Delta[-q_{10}(\log_{10}\eta_N)]$ between $a/M=0$ and $a/M=1$.",
        "Each one-factor check holds the other settings at the stated reference value.",
        "",
        "| test | N | grid | half-width | endpoint change | sign |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['test']} | {row['spectral_N']} | {row['grid_size']} | "
            f"{float(row['half_width']):.3f} | {float(row['endpoint_gain']):.4f} | "
            f"{'positive' if row['positive_gain'] else 'non-positive'} |"
        )
    lines += [
        "",
        f"The sign is positive in {sum(gain > 0 for gain in gains)}/{len(gains)} configurations; "
        f"the observed changes span {min(gains):.4f}--{max(gains):.4f}.",
        "The sign and order of magnitude, rather than a single window-dependent decimal,",
        "are the defensible robustness statement.",
    ]
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {output}")
    print(f"Wrote {report}")


if __name__ == "__main__":
    main()
