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

    # Color a protected LaTeX environment as one unit whenever any line inside
    # it changed.  This keeps color switches outside tabular, verbatim, and
    # display-math environments while ensuring revisions within those blocks
    # are still visibly marked.
    protected_names = ("tabular", "tabularx", "verbatim", "equation", "equation*", "align", "align*")
    protected_ranges: list[tuple[int, int]] = []
    stack: list[tuple[str, int]] = []
    for index, line in enumerate(new_lines):
        for name in protected_names:
            if f"\\begin{{{name}}}" in line:
                stack.append((name, index))
                break
        for name in protected_names:
            if f"\\end{{{name}}}" in line:
                for position in range(len(stack) - 1, -1, -1):
                    open_name, start = stack[position]
                    if open_name == name:
                        del stack[position]
                        protected_ranges.append((start, index + 1))
                        break
                break
    for start, end in protected_ranges:
        if any(index in changed for index in range(start, end)):
            changed.update(range(start, end))

    output: list[str] = []
    document_started = False
    blue_active = False
    for index, line in enumerate(new_lines):
        should_be_blue = index in changed and document_started
        if should_be_blue and not blue_active:
            output.append("\\color{blue}\n")
            blue_active = True
        elif not should_be_blue and blue_active:
            output.append("\\color{black}\n")
            blue_active = False
        output.append(line)
        if "\\begin{document}" in line:
            document_started = True

    if blue_active:
        output.append("\\color{black}\n")

    return "\\def\\markedrevision{1}\n" + "".join(output)


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
