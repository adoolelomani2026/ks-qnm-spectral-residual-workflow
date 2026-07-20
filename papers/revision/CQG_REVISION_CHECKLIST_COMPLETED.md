# CQG revised-submission checklist

Manuscript: CQG-116364
Prepared: 21 July 2026

## Author response

- [x] Point-by-point response addresses every substantive referee comment.
- [x] Each response states the change and gives marked-manuscript page/line locations.
- [ ] Confirm whether the submission uses double-anonymous review; anonymize if required.
- [ ] Upload on Author Centre Step 1, “View and Respond to Decision Letter.”

File: `papers/revision/CQG_RESPONSE_TO_REFEREES.pdf`

## Highlighted PDF

- [x] PDF format.
- [x] Revised text is blue and the manuscript is line numbered.
- [x] Figures and tables are included.
- [x] All 28 pages were visually inspected after the final compile; no clipped
  tables, broken paths, blank pages, or unreadable legends were found.
- [ ] Confirm anonymization requirements.
- [ ] Designate “Complete Document for Review (PDF Only)” during upload.

File: `papers/manuscript/hybrid_qnm_research_paper_marked.pdf`

## Clean source file

- [x] TeX source is clean when built without the marked-revision macro.
- [x] Full author name, affiliation, ORCID, corresponding-author marker, and email are present.
- [x] Funding, competing-interest, data-availability, code-availability, and author-contribution statements are present.
- [x] Tables, captions, and equations are editable TeX.
- [x] Tables do not rely on colored text.
- [x] No ethics statement is required for this computational study.
- [ ] Compare the author name, affiliation, ORCID, and email with the Author Centre fields.

File: `papers/manuscript/hybrid_qnm_research_paper.tex`

## Additional source files

- [x] High-resolution figure files are present under `papers/manuscript/figures/`.
- [x] Reproduction scripts, configuration, and generated tables are present.
- [ ] Upload relevant TeX dependencies and figures on Step 3 with designation “Source Files.”

## Clean PDF

- [x] Built from the clean revised TeX source.
- [x] Unmarked and without line numbers.
- [x] Compiled twice without LaTeX warnings or unresolved references.
- [x] Re-ran every documented generator and all tests from a fresh sparse clone with no pre-existing `outputs/` directory.
- [ ] Upload with designation “Source Files.”

File: `papers/manuscript/hybrid_qnm_research_paper.pdf`

## Supplementary material

- [x] No separate publication-as-is supplementary file is currently required.
- [ ] Confirm this choice in the Author Centre.

## Final authenticated checks

- [ ] Verify the 3 September 2026 deadline in the authenticated Author Centre.
- [ ] Confirm all submission metadata and file designations.
- [x] Create local `cqg-r2-final` tag at the exact resubmission-package commit.
- [x] Record the source commit used for clean regeneration in `run_metadata.json`; the annotated release tag resolves the final artifact commit.
- [ ] Push the release tags only after author approval.
