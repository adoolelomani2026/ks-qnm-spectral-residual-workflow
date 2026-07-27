#!/usr/bin/env python3
"""Fail if the CQG submission package contains stale or duplicate artifacts."""

from __future__ import annotations

import hashlib
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PAPERS = ROOT / "papers"
SUBMISSION = PAPERS / "CQG_116364_SUBMIT_THESE_5_FILES"
OUTER_ZIP = PAPERS / "CQG_116364_SUBMIT_THESE_5_FILES.zip"
SOURCE_TREE = PAPERS / "resubmission" / "source"
MANUSCRIPT = PAPERS / "manuscript"

SUBMISSION_FILES = {
    "CQG_116364_clean_revised_manuscript.pdf",
    "CQG_116364_cover_letter.pdf",
    "CQG_116364_marked_revised_manuscript.pdf",
    "CQG_116364_response_to_referees.docx",
    "CQG_116364_source_files.zip",
}

FIGURE_MAP = {
    "figure1_spectral_convergence.png": "spectral_convergence_to_cf.png",
    "figure2_mode_trajectories_gravitational_l2.png":
        "mode_trajectories_gravitational_l2.png",
    "figure2a_axial_potential_comparison.png": "axial_potential_comparison.png",
    "figure3_catalogue_l2_fractional_shifts.png":
        "catalogue_l2_fractional_shifts.png",
    "figure4_catalogue_l2_spectroscopic_ratios.png":
        "catalogue_l2_spectroscopic_ratios.png",
    "figure5_scalar_l2_pseudospectrum_contours.png":
        "scalar_l2_pseudospectrum_contours.png",
    "figure6_scalar_l2_pseudospectrum_sensitivity.png":
        "scalar_l2_pseudospectrum_sensitivity.png",
    "figure7_scalar_l2_pseudospectrum_resolution_check.png":
        "scalar_l2_pseudospectrum_resolution_check.png",
}

SOURCE_FILES = {
    "hybrid_qnm_research_paper.tex",
    *(f"figures/{name}" for name in FIGURE_MAP),
}


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"FAIL: {message}")


def relative_files(directory: Path) -> set[str]:
    return {
        path.relative_to(directory).as_posix()
        for path in directory.rglob("*")
        if path.is_file()
    }


def verify_zip_exact(archive: Path, expected: set[str]) -> dict[str, bytes]:
    with zipfile.ZipFile(archive) as bundle:
        names = {name.rstrip("/") for name in bundle.namelist() if not name.endswith("/")}
        require(names == expected, f"{archive.name} entries are {sorted(names)}")
        return {name: bundle.read(name) for name in expected}


def main() -> None:
    require(
        relative_files(SUBMISSION) == SUBMISSION_FILES,
        "submission folder must contain exactly the five upload files",
    )
    require(
        relative_files(SOURCE_TREE) == SOURCE_FILES,
        "source tree contains build debris, duplicates, or missing inputs",
    )

    source_zip = SUBMISSION / "CQG_116364_source_files.zip"
    packed_source = verify_zip_exact(source_zip, SOURCE_FILES)
    for name in SOURCE_FILES:
        disk_data = (SOURCE_TREE / name).read_bytes()
        require(
            digest(packed_source[name]) == digest(disk_data),
            f"stale source ZIP entry: {name}",
        )

    manuscript_tex = (MANUSCRIPT / "hybrid_qnm_research_paper.tex").read_bytes()
    source_tex = (SOURCE_TREE / "hybrid_qnm_research_paper.tex").read_bytes()
    require(
        digest(manuscript_tex) == digest(source_tex),
        "packaged TeX differs from the authoritative manuscript TeX",
    )
    disclosure_text = source_tex.lower()
    for forbidden in (b"openai", b"codex", b"gpt-5", b"model identifier"):
        require(
            forbidden not in disclosure_text,
            f"forbidden model/vendor wording remains in manuscript: {forbidden!r}",
        )
    require(
        b"artificial-intelligence tools were used" in disclosure_text,
        "generic AI-use disclosure is missing from the manuscript",
    )

    for packaged_name, generated_name in FIGURE_MAP.items():
        generated = (ROOT / "outputs" / "figures" / generated_name).read_bytes()
        manuscript = (MANUSCRIPT / "figures" / packaged_name).read_bytes()
        source = (SOURCE_TREE / "figures" / packaged_name).read_bytes()
        require(
            digest(generated) == digest(manuscript) == digest(source),
            f"figure mismatch: {packaged_name}",
        )

    pdf_pairs = {
        "CQG_116364_clean_revised_manuscript.pdf":
            "hybrid_qnm_research_paper.pdf",
        "CQG_116364_marked_revised_manuscript.pdf":
            "hybrid_qnm_research_paper_marked.pdf",
    }
    for packaged_name, manuscript_name in pdf_pairs.items():
        require(
            digest((SUBMISSION / packaged_name).read_bytes())
            == digest((MANUSCRIPT / manuscript_name).read_bytes()),
            f"stale packaged PDF: {packaged_name}",
        )

    outer = verify_zip_exact(OUTER_ZIP, SUBMISSION_FILES)
    for name in SUBMISSION_FILES:
        require(
            digest(outer[name]) == digest((SUBMISSION / name).read_bytes()),
            f"stale outer ZIP entry: {name}",
        )

    duplicate_source_zip = PAPERS / "resubmission" / "CQG_116364_source_files.zip"
    require(
        digest(duplicate_source_zip.read_bytes()) == digest(source_zip.read_bytes()),
        "resubmission source ZIP differs from the upload source ZIP",
    )

    response_docx = SUBMISSION / "CQG_116364_response_to_referees.docx"
    with zipfile.ZipFile(response_docx) as response:
        response_xml = b"".join(
            response.read(name)
            for name in response.namelist()
            if name.endswith(".xml")
        ).lower()
    for forbidden in (b"openai", b"codex", b"gpt-5", b"model identifier"):
        require(
            forbidden not in response_xml,
            f"forbidden model/vendor wording remains in response: {forbidden!r}",
        )

    print("PASS: submission folder contains exactly five files")
    print("PASS: source tree and source ZIP contain exactly one TeX file and eight figures")
    print("PASS: packaged TeX, figures, PDFs, and nested ZIP bytes are current")
    print("PASS: active submission files contain only the generic AI disclosure")


if __name__ == "__main__":
    main()
