#!/usr/bin/env python3
"""Compare the scalar ell=0 fundamental side diagnostic with Konoplya 2020.

Konoplya's published table fixes the horizon radius, r_h=1, and reports
frequencies in the horizon-scaled convention r_h omega. The main catalogue in
this repository is fixed-mass and starts at ell=2, so this script is deliberately
a side diagnostic rather than part of the publication-facing catalogue.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from qnm.leaver import solve_leaver_mode  # noqa: E402
from qnm.normalization import alpha_from_horizon_beta, rhomega_from_momega  # noqa: E402
from qnm.spectral import build_spectral_problem, generalized_eigenvalues  # noqa: E402
from qnm.common import select_physical_mode  # noqa: E402


DEFAULT_INPUT = ROOT_DIR / "data" / "literature" / "konoplya2020_scalar_l0_table.csv"
DEFAULT_OUTPUT = ROOT_DIR / "outputs" / "results" / "konoplya2020_scalar_l0_comparison.csv"

OUTPUT_COLUMNS = [
    "source",
    "reference",
    "a_over_rh",
    "a_over_M",
    "spectral_N",
    "this_work_rhomega_real",
    "this_work_rhomega_imag",
    "konoplya_wkb_rhomega_real",
    "konoplya_wkb_rhomega_imag",
    "konoplya_time_domain_rhomega_real",
    "konoplya_time_domain_rhomega_imag",
    "relative_difference_vs_wkb",
    "relative_difference_vs_time_domain",
    "continued_fraction_abs",
    "notes",
]


def _truthy(value: str | None) -> bool:
    return (value or "").strip().lower() in {"1", "true", "yes", "y"}


def _float(value: str | None, field: str) -> float:
    stripped = (value or "").strip()
    if not stripped:
        raise ValueError(f"Missing required field {field}.")
    return float(stripped)


def _optional_complex(real: str | None, imag: str | None) -> complex | None:
    if not (real or "").strip() and not (imag or "").strip():
        return None
    return complex(_float(real, "real"), _float(imag, "imag"))


def build_rows(input_path: Path, spectral_n: int) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    target = 0.110455 - 0.104896j

    with input_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if not _truthy(row.get("enabled")):
                continue

            beta = _float(row.get("a_over_rh"), "a_over_rh")
            alpha = alpha_from_horizon_beta(beta)
            wkb_rh = complex(
                _float(row.get("wkb_real_rh"), "wkb_real_rh"),
                _float(row.get("wkb_imag_rh"), "wkb_imag_rh"),
            )
            td_rh = _optional_complex(row.get("time_domain_real_rh"), row.get("time_domain_imag_rh"))

            problem = build_spectral_problem(alpha, spectral_n, ell=0, perturbation_type="scalar")
            eigenvalues = generalized_eigenvalues(problem)
            guess = select_physical_mode(eigenvalues, target)
            leaver_omega, cf_abs = solve_leaver_mode(
                alpha,
                guess,
                ell=0,
                perturbation_type="scalar",
                residual_tolerance=1.0e-4,
            )
            target = leaver_omega
            this_work_rh = rhomega_from_momega(alpha, leaver_omega)

            rel_wkb = abs(this_work_rh - wkb_rh) / abs(wkb_rh)
            rel_td = "" if td_rh is None else abs(this_work_rh - td_rh) / abs(td_rh)
            rows.append(
                {
                    "source": row["source"],
                    "reference": row["reference"],
                    "a_over_rh": f"{beta:.12g}",
                    "a_over_M": f"{alpha:.12g}",
                    "spectral_N": str(spectral_n),
                    "this_work_rhomega_real": f"{this_work_rh.real:.12g}",
                    "this_work_rhomega_imag": f"{this_work_rh.imag:.12g}",
                    "konoplya_wkb_rhomega_real": f"{wkb_rh.real:.12g}",
                    "konoplya_wkb_rhomega_imag": f"{wkb_rh.imag:.12g}",
                    "konoplya_time_domain_rhomega_real": "" if td_rh is None else f"{td_rh.real:.12g}",
                    "konoplya_time_domain_rhomega_imag": "" if td_rh is None else f"{td_rh.imag:.12g}",
                    "relative_difference_vs_wkb": f"{rel_wkb:.12g}",
                    "relative_difference_vs_time_domain": "" if rel_td == "" else f"{rel_td:.12g}",
                    "continued_fraction_abs": f"{cf_abs:.12g}",
                    "notes": row.get("notes", ""),
                }
            )
    return rows


def write_rows(output_path: Path, rows: list[dict[str, str]]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--spectral-n", type=int, default=64)
    args = parser.parse_args()

    rows = build_rows(args.input, args.spectral_n)
    write_rows(args.output, rows)
    print(f"Wrote {len(rows)} Konoplya 2020 comparison row(s) to {args.output}.")


if __name__ == "__main__":
    main()
