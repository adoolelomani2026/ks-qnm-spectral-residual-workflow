# CQG-116364 resubmission package

Prepared on 21 July 2026.

This build replaces the earlier lapse-substitution axial proxy with a
gauge-invariant sourced odd-parity construction under an explicitly added
inverse-Cowling effective-source closure. It includes potential-positivity and
growing-root audits, Taylor-order and `N=128` convergence, raw-polynomial
backward errors, and an 18-mode comparison with the public 2026
Batic--Dutykh--Sukaiti calculation. The scalar sector remains the least
assumption-dependent result. The final preflight also adds a documented
candidate-selection algorithm and a seven-configuration pseudospectral
grid/window/resolution analysis. Both 26-page manuscript PDFs were visually
inspected page by page.

- `CQG-116364_response_to_referees.pdf`: point-by-point response with page and line locations.
- `CQG-116364_marked_revised_manuscript.pdf`: blue, line-numbered review copy.
- `CQG-116364_clean_revised_manuscript.pdf`: clean revised manuscript.
- `CQG-116364_revision_checklist.pdf`: completed working checklist with authenticated/manual items left open.
- `source/hybrid_qnm_research_paper.tex`: clean source; marked output is activated only when the `markedrevision` macro is defined.
- `source/figures/`: manuscript figure files.
- `CQG-116364_source_files.zip`: upload-ready TeX and figure bundle.

The source folder was compiled independently. Auxiliary build files in that
folder are not submission inputs and should not be uploaded.
Every documented generator and test command was also run from a fresh sparse
clone with no pre-existing `outputs/` directory; the regenerated metadata
records the clean source commit used for that pass.

Manual actions still required before upload:

1. Verify the 3 September 2026 deadline in the authenticated Author Centre.
2. Confirm author details and whether double-anonymous review is enabled.
3. Apply the journal's required file designations during upload.
4. Push the existing local `cqg-r1` and `cqg-major-revision-r1` tags only after author approval.
