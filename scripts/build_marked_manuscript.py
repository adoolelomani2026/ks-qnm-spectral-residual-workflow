#!/usr/bin/env python3
"""Build a line-numbered TeX review copy with changed lines colored blue."""

from __future__ import annotations

import argparse
import difflib
import subprocess
from pathlib import Path


DEFAULT_SOURCE = Path("papers/manuscript/hybrid_qnm_research_paper.tex")
DEFAULT_OUTPUT = Path("papers/manuscript/hybrid_qnm_research_paper_marked.tex")


def submitted_text(commit: str, source: Path) -> str:
    result = subprocess.run(
        ["git", "show", f"{commit}:{source.as_posix()}"],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    return result.stdout


def mark_changed_lines(original: str, revised: str) -> str:
    old_lines = original.splitlines(keepends=True)
    new_lines = revised.splitlines(keepends=True)
    matcher = difflib.SequenceMatcher(None, old_lines, new_lines, autojunk=False)
    changed: set[int] = set()
    for opcode, _old_start, _old_end, new_start, new_end in matcher.get_opcodes():
        if opcode != "equal":
            changed.update(range(new_start, new_end))

    output: list[str] = []
    document_started = False
    tabular_depth = 0
    verbatim_depth = 0
    math_depth = 0
    blue_active = False
    for index, line in enumerate(new_lines):
        enters_tabular = "\\begin{tabular" in line
        leaves_tabular = "\\end{tabular" in line
        enters_verbatim = "\\begin{verbatim}" in line
        leaves_verbatim = "\\end{verbatim}" in line
        enters_math = any(f"\\begin{{{name}}}" in line for name in ("equation", "equation*", "align", "align*"))
        leaves_math = any(f"\\end{{{name}}}" in line for name in ("equation", "equation*", "align", "align*"))
        unsafe = (
            tabular_depth > 0
            or verbatim_depth > 0
            or enters_tabular
            or leaves_tabular
            or enters_verbatim
            or leaves_verbatim
            or math_depth > 0
            or enters_math
            or leaves_math
        )
        should_be_blue = index in changed and document_started and not unsafe
        if should_be_blue and not blue_active:
            output.append("\\color{blue}\n")
            blue_active = True
        elif not should_be_blue and blue_active:
            output.append("\\color{black}\n")
            blue_active = False
        output.append(line)
        if "\\begin{document}" in line:
            document_started = True
        if enters_tabular:
            tabular_depth += 1
        if leaves_tabular:
            tabular_depth = max(0, tabular_depth - 1)
        if enters_verbatim:
            verbatim_depth += 1
        if leaves_verbatim:
            verbatim_depth = max(0, verbatim_depth - 1)
        if enters_math:
            math_depth += 1
        if leaves_math:
            math_depth = max(0, math_depth - 1)

    if blue_active:
        output.append("\\color{black}\n")

    return "".join(output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--commit", default="9090c8c")
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    original = submitted_text(args.commit, args.source)
    revised = args.source.read_text(encoding="utf-8")
    args.output.write_text(mark_changed_lines(original, revised), encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
