# EPJ Plus Submission Checklist

Target manuscript: *A Chebyshev-Leaver Spectral Residual Workflow for Kazakov-Solodukhin Quasinormal Modes*

## Implemented manuscript requirements

- Title page includes the concise title, author name, full affiliation, corresponding-author email, and ORCID.
- Abstract has a visible heading and is 227 words, within the EPJ Plus 150-250 word range.
- Abstract defines QNM and KS at first use.
- Abstract includes the requested core results:
  - Schwarzschild scalar `ell=2` validation relative error: `1.462e-10`.
  - KS scalar `ell=2,n=0` endpoint shifts at `a/M=1`: `-7.29%` in `Re(M omega)`, `-2.59%` in `-Im(M omega)`, and `-4.83%` in quality factor.
  - Worst Leaver-spectral relative difference: `1.097e-5`.
  - Axial gravitational caveat: phenomenological KS-lapse-deformed Regge-Wheeler sector.
- Keywords were added below the abstract.
- Section headings use decimal numbering and no more than three levels.
- The weak "Frequency comparison" subsection now contains explanatory text and citations to its table and figure.
- Figures and tables use Arabic numbering through LaTeX and journal-style `Fig.`/`Table` captions.
- Every figure and table label is cited in the manuscript.
- Figure-generation code was updated to remove embedded plot titles where avoidable.
- Numbered equations are labeled and referenced; supporting displayed definitions are unnumbered.
- References are numbered consecutively in a single-reference-per-number bibliography.
- Citation audit passed: 20 citation keys, 20 bibliography entries, no missing or uncited bibliography items.
- DOI links are included for all bibliography entries.
- Reproducibility section now names the spectral solver, Leaver-style validation layer, catalogue-generation scripts, generated tables, automated tests, and repository.
- Declarations section includes:
  - Funding.
  - Competing interests.
  - Data availability.
  - Code availability.
  - Author contribution.
  - Ethics approval.
  - Consent to participate.
  - Consent for publication.
  - AI-tool use disclosure for Springer/EPJ transparency.
- Claim discipline was preserved:
  - No invention claim for residual singular-value methods.
  - No quantum-computing claim.
  - No detectability forecast.
  - Axial gravitational caveat retained.
  - Overtone caveats retained.

## Repository and files

- Repository link used in manuscript:
  <https://github.com/adoolelomani2026/ks-qnm-spectral-residual-workflow>
- Revised LaTeX manuscript:
  `papers/manuscript/hybrid_qnm_research_paper.tex`
- Compiled PDF:
  `papers/manuscript/hybrid_qnm_research_paper.pdf`

## Validation performed

- Full pipeline regenerated generated tables and figures with `python scripts/run_hybrid_qnm_algorithm.py`.
- Fast automated tests passed with `python -m pytest`.
- Manuscript compiled successfully with `pdflatex` run twice.
- LaTeX log scan found no unresolved references, undefined citations, fatal errors, or overfull boxes.
- Abstract word-count check returned 227 words.
- Figure/table label audit found all figures and tables cited.
- Equation audit found all numbered equations labeled and referenced.

## Remaining submission-portal items

- No repository placeholder remains; the GitHub link is active.
- The official Springer Nature `sn-jnl.cls` class was not available in the local TeX installation, so the manuscript keeps the existing LaTeX article layout while implementing EPJ/Springer structural requirements.
- EPJ Plus may request a separate graphical abstract upload through the submission portal. That file is not part of the manuscript PDF.
