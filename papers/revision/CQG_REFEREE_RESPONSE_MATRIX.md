# CQG-116364 Major-Revision Response Matrix

Revision branch: `cqg-major-revision`

Submitted-version baseline: Git commit `9090c8c` (`Update manuscript document`).
The submitted-version files at that commit are not to be modified; all revision
work is made on the revision branch. The decision letter gives a resubmission
deadline of **3 September 2026**. Verification in the authenticated Author
Centre remains a manual action.

Strategic spine:

> The scalar low-lying QNMs remain the least assumption-dependent headline;
> the axial metric variable is gauge invariant under an explicitly added
> inverse-Cowling effective-source closure; the numerical novelty is auditability and validation, not a new
> high-overtone solver.

| Referee comment | Required action | Manuscript location | Status |
| --- | --- | --- | --- |
| Referee 1.1: Equations (1), (3), and (4) lack sufficient provenance or derivation, especially the axial-gravitational potential. | Cite and derive the metric and scalar potential; give the sourced Gerlach--Sengupta/Karlovini equation and exact equation numbers; define the invariant current and added inverse-Cowling closure; derive the matter-corrected KS potential and Schwarzschild limit. | Section 2.1, background, scalar, and odd-parity derivations | Complete; corrected 1979 citation and convention mapping included |
| Referee 1.2: The Schwarzschild reference value and Leaver benchmark procedure are insufficiently documented, and reproducing a standard value does not establish methodological advantage. | Cite the exact reference frequency and conventions; document all continued-fraction and root-solver settings; add truncation convergence; explain that novelty is auditability, cross-validation, backward-error diagnostics, and traceable branch acceptance rather than raw precision. | Introduction; numerical-method subsection; Section 3.1; new convergence table | Complete; compiled and reproduced |
| Referee 1.3: The method should not be presented as solving the general high-overtone problem; recent Chen--Heun work should be acknowledged. | Adopt explicit confidence tiers; describe the catalogue as low lying; keep second overtones diagnostic; state higher overtones are out of scope; add a respectful Chen--Heun comparison and citations. | Abstract; Introduction; catalogue results; limitations; conclusion | Complete; compiled |
| Referee 1.4: Axial results require a gauge-invariant derivation or must remain phenomenological. | Define the odd gauge transformation and invariant master variable, include the sourced equation and effective matter term, and distinguish gauge invariance from the separate inverse-Cowling closure. | Section 2.1; axial results; captions; limitations; conclusion | Complete; closure is explicitly an added model assumption |
| Referee 2.1: Clarify that simple lapse substitution is not first-principles perturbation theory. | Go beyond the requested disclaimer by removing the lapse-substitution model; derive the gauge-invariant frozen-source equation and display its correction relative to the old proxy. | Section 2.1, directly after the axial derivation | Complete; correction vanishes only in Schwarzschild limit |
| Referee 2.2: Clarify whether higher-overtone errors indicate diagnostics or equal-confidence predictions. | State directly beside the catalogue table that first overtones are validated low-lying modes, second overtones are exploratory diagnostics, and higher overtones are outside scope; distinguish literature-reference error from spectral--Leaver disagreement. | Catalogue table and accompanying text; confidence-tier table | Complete; catalogue CSV also records tiers |
| Referee 2.3: Explain the slightly larger scalar `ell=4` endpoint shift. | Verify all shift definitions; compute scalar-potential peak radius, height, and tortoise-coordinate curvature for `ell=2,3,4`; use WKB/eikonal intuition only if supported, otherwise describe the difference as a modest catalogue trend. | Catalogue-level trends subsection; new supporting table and repository figure | Complete; analysis regenerated |
| Referee 2.4: Consider plotting the quality factor against deformation. | Add `Q=Re(omega)/(2[-Im(omega)])` versus `a/M` for the scalar `ell=2` fundamental, with optional visibly qualified overtone curves; use the same definition in text and captions. | Figure 4 and spectroscopy subsection | Complete; Figure 4 regenerated |
| Referee 2.5: Clarify the physical significance of `a/M=1`. | State that it is the upper endpoint of the sampled interval, not extremality or horizon degeneracy; verify the original convention and explain the range choice; replace ambiguous uses of "endpoint." | Section 2.1; literature-normalization subsection; results and captions | Complete; compiled |
| Referee 2.6: Expand the introduction's overview of standard QNM methods. | Add a concise background paragraph covering WKB, time-domain evolution, Frobenius/Leaver, AIM, spectral methods, pseudospectra, and Chen--Heun/high-overtone approaches; verify the suggested references. | Introduction, literature context | Complete; representative sources verified |
| Referee 2.7: Mathematical symbols in the abstract are not consistently typeset. | Rewrite the abstract after scientific revisions and put every mathematical expression in math mode, including `N`, `ell`, `n`, `a/M`, `Q`, and scientific notation. | Abstract | Complete; compiled |
| Referee 2.8: Capitalize "Figure" consistently. | Replace lowercase prose references with "Figure"; audit all figure, table, and equation cross-references and scientific hyphenation. | Entire manuscript | Complete; text audit passed |
| Adversarial theory audit: foundational citation, sourced equation, closure provenance, and public overlap must be unambiguous. | Correct Gerlach--Sengupta to PRD 19, 2268 (1979); map harmonic and epsilon conventions; display source conservation and `L_a=L=0`; cite and compare the public Batic--Dutykh--Sukaiti calculation without a priority claim. | Section 2.1; axial results; references | Complete; 18-mode external comparison generated |
| Adversarial numerical audit: replace weak conditioning and self-convergence evidence. | Replace near-singular `kappa(P)` with raw-polynomial backward error; plot error to CF through `N=128`; add Taylor orders 64--128 and first-overtone deterioration; remove uneven literature-error column. | Solver validation; Figures 1 and 2a; Tables 3, 4, and residual diagnostics | Complete; scripts and CSV evidence added |
| Adversarial stability audit: quantify the potential change and test for instability. | Plot new versus old potentials, minimize the new potential for all sampled multipoles/deformations, and search raw spectra for growing roots at three resolutions. | Axial results and repository audit report | Complete; non-negative potential and zero surviving growing candidates |
| Final-preflight reproducibility audit: raw generalized eigenvalues must map to branches through an operational algorithm. | Record the frequency window, non-finite filtering, duplicate clustering, equilibration, continuation gate, match score, minimizer settings, backward-error threshold, CF threshold, and overtone assignment. | Candidate filtering, refinement, and branch assignment subsection | Complete; values synchronized with code |
| Final-preflight pseudospectrum audit: the exact local quantile number is window dependent. | Rename the statistic to `q_p`; vary grids `81^2,121^2,161^2`, half-widths `0.020,0.025,0.030`, and `N=32,48,64`; report the robust sign and magnitude range rather than one decimal in the abstract. | Pseudospectrum subsection and robustness table | Complete; positive in 7/7 configurations, range 0.1002--0.1647 |
| Final-preflight terminology and presentation audit. | Replace legacy conditioning language; use collocation-independent CF wording; separate scalar-QNM and GW language; standardize axial inverse-Cowling figure labels; remove manuscript priority chronology; disambiguate `Q`; rename Table 4. | Manuscript-wide; figures; references | Complete; text and figure sources updated |
| Editorial decision: Submit all required revision files and ensure author metadata agree. | Prepare a point-by-point response, highlighted review PDF, clean TeX source, additional source/figure files, and clean PDF; check author name, affiliation, email, funding, ethics, data, code, and contribution statements against the submission form. | Resubmission package and Author Centre | Package assembled; authenticated metadata check remains manual |

## Repository and numerical work linked to the matrix

| Work item | Referee link | Status |
| --- | --- | --- |
| Continued-fraction truncation convergence script and table | R1.2 | Complete |
| Taylor-order convergence and Chebyshev-to-CF `N=16...128` plots | R1.2 | Complete |
| Axial potential/positivity/growing-root/public-data audit | R1.1, R1.4, R2.1 | Complete |
| Polynomial backward-error diagnostic on raw matrices | R1.2 | Complete |
| Raw-candidate rejection/auditability example, or an explicit traceability-only scope statement if no defensible example is found | R1.2 | Complete via explicit traceability-only scope statement; no contrived rejection example added |
| Confidence-tier column/status in catalogue outputs | R1.3, R2.2 | Complete |
| Higher-precision Schwarzschild references and reference-precision audit | R1.2, R2.2 | Complete for headline scalar reference; uneven axial literature-error column removed |
| Scalar-potential peak analysis | R2.3 | Complete |
| Quality-factor plot | R2.4 | Complete |
| Regenerated tables, figures, and metadata | All numerical comments | Complete |
| Full automated test suite and clean independent LaTeX compilation | Editorial | Complete: 8/8 fast checks and full validation passed; 27-page PDF compiled without warnings |
| Local revision tags `cqg-r1` and `cqg-major-revision-r1` at the final-preflight commit | Editorial | Complete; remote push remains pending author approval |
| Final-preflight candidate-score notation and mathematical audit | Reproducibility | Complete; distances, shortlist, `sigma_best`, pre-refinement residual, equilibration, eigenvector mapping, normalization, phase invariance, and raw backward-error gate stated explicitly |
| Final-preflight confidence and continued-fraction terminology | Claim discipline | Complete; robust/cross-validated/exploratory/conditional labels synchronized and collocation-independent wording used |
| Final-preflight external-reference and AI-disclosure audit | Editorial integrity | Complete; no persistent identifier found for the external work as of 20 July 2026; exact public commit retained; exact Codex model identifier verified from session metadata |
| Fresh-clone bootstrap and command audit | Reproducibility | Complete; every documented generator and test passed with no pre-existing `outputs/` tree, and recursive output-directory creation was fixed |

## Submission files

| Required file | Status |
| --- | --- |
| Point-by-point author response | Complete with marked-manuscript page/line locations |
| Highlighted revised PDF | Complete: blue, line-numbered review copy |
| Clean revised TeX source | Complete |
| Additional TeX components and high-resolution figures | Complete |
| Clean revised PDF | Complete |
| Supplementary material, if designated for publication | Not designated; repository remains the reproducibility record |
