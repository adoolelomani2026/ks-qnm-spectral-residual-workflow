#!/usr/bin/env python3
"""Normalize transcribed KS literature QNM rows and compare to the catalogue.

The input CSV is intentionally a human-transcribed table. The script does not
scrape papers; it makes the convention conversion explicit and reproducible.
Rows are ignored unless ``enabled`` is true.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from qnm.normalization import (  # noqa: E402
    alpha_from_horizon_beta,
    horizon_beta_from_alpha,
    momega_from_rhomega,
    rhomega_from_momega,
)


DEFAULT_INPUT = ROOT_DIR / "data" / "literature" / "ks_qnm_literature_template.csv"
DEFAULT_CATALOGUE = ROOT_DIR / "outputs" / "results" / "qnm_catalogue.csv"
DEFAULT_OUTPUT = ROOT_DIR / "outputs" / "results" / "literature_normalized_comparison.csv"

OUTPUT_COLUMNS = [
    "source",
    "reference",
    "perturbation_type",
    "ell",
    "overtone",
    "mode",
    "literature_a_over_M",
    "literature_a_over_rh",
    "literature_Momega_real",
    "literature_Momega_imag",
    "catalogue_a_over_M",
    "catalogue_Momega_real",
    "catalogue_Momega_imag",
    "parameter_abs_difference",
    "frequency_abs_difference",
    "relative_difference",
    "comparison_status",
    "notes",
]


@dataclass(frozen=True)
class LiteratureRow:
    source: str
    reference: str
    perturbation_type: str
    ell: int
    overtone: int
    mode: str
    a_over_m: float
    a_over_rh: float
    omega_m: complex
    notes: str


@dataclass(frozen=True)
class CatalogueRow:
    perturbation_type: str
    ell: int
    overtone: int
    mode: str
    a_over_m: float
    omega_m: complex


def _truthy(value: str | None) -> bool:
    return (value or "").strip().lower() in {"1", "true", "yes", "y"}


def _optional_float(value: str | None) -> float | None:
    stripped = (value or "").strip()
    if not stripped:
        return None
    return float(stripped)


def _require_float(value: str | None, field: str, row_number: int) -> float:
    parsed = _optional_float(value)
    if parsed is None:
        raise ValueError(f"Row {row_number}: missing required field {field!r}.")
    return parsed


def _read_literature_rows(path: Path) -> list[LiteratureRow]:
    rows: list[LiteratureRow] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row_number, row in enumerate(reader, start=2):
            if not _truthy(row.get("enabled")):
                continue

            alpha = _optional_float(row.get("a_over_M"))
            beta = _optional_float(row.get("a_over_rh"))
            if alpha is None and beta is None:
                raise ValueError(f"Row {row_number}: provide a_over_M or a_over_rh.")
            if alpha is not None and beta is None:
                beta = horizon_beta_from_alpha(alpha)
            if beta is not None and alpha is None:
                alpha = alpha_from_horizon_beta(beta)

            scale = (row.get("frequency_scale") or "").strip().lower()
            omega = complex(
                _require_float(row.get("omega_real"), "omega_real", row_number),
                _require_float(row.get("omega_imag"), "omega_imag", row_number),
            )
            if scale in {"m", "mass", "momega", "m_omega"}:
                omega_m = omega
            elif scale in {"rh", "r_h", "horizon", "rhomega", "rh_omega"}:
                assert beta is not None
                omega_m = momega_from_rhomega(beta, omega)
            else:
                raise ValueError(
                    f"Row {row_number}: frequency_scale must be one of M or rh; got {scale!r}."
                )

            assert alpha is not None and beta is not None
            rows.append(
                LiteratureRow(
                    source=(row.get("source") or "").strip(),
                    reference=(row.get("reference") or "").strip(),
                    perturbation_type=(row.get("perturbation_type") or "").strip(),
                    ell=int((row.get("ell") or "").strip()),
                    overtone=int((row.get("overtone") or "").strip()),
                    mode=(row.get("mode") or "").strip(),
                    a_over_m=alpha,
                    a_over_rh=beta,
                    omega_m=omega_m,
                    notes=(row.get("notes") or "").strip(),
                )
            )
    return rows


def _read_catalogue(path: Path) -> list[CatalogueRow]:
    rows: list[CatalogueRow] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(
                CatalogueRow(
                    perturbation_type=row["perturbation_type"],
                    ell=int(row["ell"]),
                    overtone=int(row["overtone"]),
                    mode=row["mode"],
                    a_over_m=float(row["a_over_M"]),
                    omega_m=complex(float(row["spectral_real"]), float(row["spectral_imag"])),
                )
            )
    return rows


def _nearest_catalogue_row(lit: LiteratureRow, catalogue: list[CatalogueRow]) -> CatalogueRow | None:
    branch = [
        row
        for row in catalogue
        if row.perturbation_type == lit.perturbation_type
        and row.ell == lit.ell
        and row.overtone == lit.overtone
    ]
    if not branch:
        return None
    return min(branch, key=lambda row: abs(row.a_over_m - lit.a_over_m))


def build_comparison_rows(
    literature_rows: list[LiteratureRow],
    catalogue_rows: list[CatalogueRow],
    *,
    exact_tolerance: float,
) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    for lit in literature_rows:
        matched = _nearest_catalogue_row(lit, catalogue_rows)
        base = {
            "source": lit.source,
            "reference": lit.reference,
            "perturbation_type": lit.perturbation_type,
            "ell": str(lit.ell),
            "overtone": str(lit.overtone),
            "mode": lit.mode,
            "literature_a_over_M": f"{lit.a_over_m:.12g}",
            "literature_a_over_rh": f"{lit.a_over_rh:.12g}",
            "literature_Momega_real": f"{lit.omega_m.real:.12g}",
            "literature_Momega_imag": f"{lit.omega_m.imag:.12g}",
            "notes": lit.notes,
        }
        if matched is None:
            output.append(
                base
                | {
                    "catalogue_a_over_M": "",
                    "catalogue_Momega_real": "",
                    "catalogue_Momega_imag": "",
                    "parameter_abs_difference": "",
                    "frequency_abs_difference": "",
                    "relative_difference": "",
                    "comparison_status": "literature_only_no_matching_branch",
                }
            )
            continue

        parameter_difference = abs(matched.a_over_m - lit.a_over_m)
        frequency_difference = abs(matched.omega_m - lit.omega_m)
        relative_difference = frequency_difference / max(abs(lit.omega_m), 1.0e-300)
        status = "exact_catalogue_match" if parameter_difference <= exact_tolerance else "nearest_catalogue_point"
        output.append(
            base
            | {
                "catalogue_a_over_M": f"{matched.a_over_m:.12g}",
                "catalogue_Momega_real": f"{matched.omega_m.real:.12g}",
                "catalogue_Momega_imag": f"{matched.omega_m.imag:.12g}",
                "parameter_abs_difference": f"{parameter_difference:.12g}",
                "frequency_abs_difference": f"{frequency_difference:.12g}",
                "relative_difference": f"{relative_difference:.12g}",
                "comparison_status": status,
            }
        )
    return output


def write_comparison(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="transcribed literature CSV")
    parser.add_argument("--catalogue", type=Path, default=DEFAULT_CATALOGUE, help="project catalogue CSV")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="normalized comparison CSV")
    parser.add_argument(
        "--exact-tolerance",
        type=float,
        default=1.0e-9,
        help="a/M tolerance for an exact catalogue parameter match",
    )
    args = parser.parse_args()

    literature_rows = _read_literature_rows(args.input)
    catalogue_rows = _read_catalogue(args.catalogue)
    comparison_rows = build_comparison_rows(
        literature_rows,
        catalogue_rows,
        exact_tolerance=args.exact_tolerance,
    )
    write_comparison(args.output, comparison_rows)

    exact = sum(row["comparison_status"] == "exact_catalogue_match" for row in comparison_rows)
    nearest = sum(row["comparison_status"] == "nearest_catalogue_point" for row in comparison_rows)
    missing = sum(row["comparison_status"] == "literature_only_no_matching_branch" for row in comparison_rows)
    print(
        "Wrote "
        f"{len(comparison_rows)} enabled literature row(s) to {args.output}. "
        f"exact={exact}, nearest={nearest}, missing_branch={missing}"
    )
    if not comparison_rows:
        print("No rows were enabled; set enabled=true after transcribing literature table entries.")


if __name__ == "__main__":
    main()
