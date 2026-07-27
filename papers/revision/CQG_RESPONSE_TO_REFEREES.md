# Response to the Referees

Manuscript: **CQG-116364**, *Low-lying quasinormal modes of Kazakov--Solodukhin black holes from Chebyshev and continued-fraction methods*

We thank the editor and referees for their careful and constructive reports. The
manuscript has been extensively revised. We have derived the minimally coupled
scalar equation, replaced the former lapse-substitution axial proxy with a
gauge-invariant sourced odd-parity formulation under an explicitly stated
inverse-Cowling closure, regenerated the complete axial catalogue, documented
the spectral candidate-selection algorithm, added continued-fraction depth and
Taylor-order convergence studies, replaced condition-number reporting with
scale-aware backward errors, added external high-precision axial comparisons,
documented all 18 Schwarzschild branch seeds, audited the positive-imaginary
raw spectrum before damped-mode filtering, clarified the confidence hierarchy
of the overtone catalogue, and added the requested physical and graphical
diagnostics. The minimally coupled scalar
sector remains the least assumption-dependent KS result, and the
methodological contribution is an audited low-lying-mode workflow rather than
a new high-overtone solver.
The final presentation now foregrounds the physical chain connecting the KS
geometry, the displacement and softening of the scalar potential barrier, and
the resulting changes in oscillation frequency, damping time, and quality
factor; the numerical audits are retained as support for that interpretation.
Definitions and physical explanations are now introduced naturally at first
use: the mode labels and WKB picture in the Introduction, the geometric and
boundary-condition language with the wave equation, gauge and closure concepts
with the odd-parity derivation, and each numerical diagnostic in the Methods
section where it enters the calculation.  The Methods and pseudospectrum
sections state explicitly what question each tool answers, why it is used, and
what it cannot establish.

We have revised the title to reflect more accurately the manuscript's low-lying spectral scope and its use of two complementary numerical methods. The new title is *Low-lying quasinormal modes of Kazakov--Solodukhin black holes from Chebyshev and continued-fraction methods*.

Page and line references below refer to the blue, line-numbered marked
manuscript.

## Referee 1

### Referee 1, Comment 1

> My first concern is about the theoretical setup in Sec. 2.1. The spacetime
> form in Eq. (1) is introduced without any reference or derivation, although it
> is not, at least to me, a completely standard expression. The same issue
> occurs for the effective potentials in Eqs. (3) and (4). The authors simply
> write down these potentials, but do not explain where they come from. This is
> particularly important for the axial-gravitational case, because the
> potential is later used more like a lapse-deformed Regge-Wheeler potential
> than a perturbation equation derived from a full gauge-invariant analysis.
> The authors should either provide a derivation or clearly cite the original
> source and state the assumptions behind these equations.

**Response:** We agree that the provenance and assumptions were previously
underexplained. Section 2.1 now cites the original Kazakov--Solodukhin
construction and a later paper using the same lapse convention, explains the
areal-radius and deformation-parameter notation without introducing an unused
renormalized coupling, states the coordinate domain \(r\geq a\), derives
the horizon relation \(r_h^2=4M^2+a^2\), and displays the Schwarzschild limit. We now derive
the scalar potential from \(\Box\Phi=0\) using
\(\Phi=e^{-i\omega t}Y_{\ell m}\Psi/r\). For the axial equation, we went beyond
retaining a phenomenological disclaimer. Following Gerlach--Sengupta,
Gundlach--Martínez-García, and Karlovini, the revision maps the harmonic and
volume-form conventions; defines the odd gauge transformation, invariant
one-form \(k_a\), invariant curl \(\Pi\), and invariant source current \(L_a\);
and displays the general sourced master equation before imposing a closure.
We corrected the foundational citation to Physical Review D 19, 2268 (1979)
and give the exact equation numbers used. Setting \(L_a=L=0\) both removes the
source and satisfies its conservation equation. This is stated as an
additional inverse-Cowling closure motivated by, but not uniquely prescribed
by, the original KS mode split. The resulting potential contains the effective
density and radial pressure and differs from simple lapse substitution for
\(a>0\). We also compare the same potential and 18 overlapping frequencies
with the published Batic--Dutykh--Sukaiti calculation at an exactly cited
repository commit.

**Changes in the manuscript:** Introduction, pages 1--3, lines 8--74;
Section 2.1, pages 4--7, lines 115--208; axial external validation and
stability analysis, pages 14--16, lines 420--461.

### Referee 1, Comment 2

> The Schwarzschild benchmark in Sec. 3.1 also needs clarification. The authors
> quote a “known” scalar Schwarzschild fundamental frequency, but it is not
> clear where this number comes from or how it was obtained. Similarly, when
> they say that “the Leaver-style solver gives” a certain value, the truncation
> order, convergence criterion, and numerical stability are not specified. For
> the scalar Schwarzschild fundamental mode, many standard methods can reach
> very high precision once the order is high enough. Therefore, reproducing
> this number does not by itself show the advantage of the present workflow.
> The authors need to explain more clearly what is actually improved by their
> method.

**Response:** We now identify the external value
\(M\omega=0.483643872211-0.096758775978i\) as the spin-zero,
\(\ell=2,n=0\) Schwarzschild result reported by Cavalcante and Carneiro da
Cunha (Table I, page 6), with \(M=1\), \(e^{-i\omega t}\), and
\(\operatorname{Im}\omega<0\). The manuscript now also defines
Δrel(ω₁, ω₂) = |ω₁ − ω₂|/|ω₂|,
where \(\omega_2\) is the reference value.
The revised methods give the actual implementation settings: Taylor order 96,
continued-fraction/Gaussian-elimination depth 240, IEEE-754 double precision,
SciPy/MINPACK hybrid root finding, root tolerance \(10^{-11}\), at most 1000
function evaluations, the residual thresholds, branch-matched \(N=32\)
initial guesses, and continuation in \(a/M\). A new depth-truncation table covers the
Schwarzschild fundamental, one KS-deformed fundamental, and a KS first
overtone. A separate Taylor-order table varies orders 64, 80, 96, 112, and
128 at fixed depth 320. Figure 1 now plots error from an order-128,
depth-320 continued-fraction root through \(N=128\), for both the fundamental
and first overtone. It exposes the overtone's high-\(N\) deterioration and
supports the more precise label “cross-validated low-lying first overtones at
\(N=32\).” The same
settings are stored in `config/leaver_revision.json`.

We also rewrote the novelty claim. Agreement with Schwarzschild is presented as
a cross-discretization check, not a claim of superior precision. The
contribution is the traceable chain from equilibrated generalized-eigenvalue
candidates through explicit filtering, clustering, continuation, branch
scoring, singular-value refinement, spectral-size convergence, polynomial
backward-error diagnostics on the raw matrices, and collocation-independent
continued-fraction validation. No illustrative rejected candidate was included
because the available audit data did not provide a representative example;
instead, the acceptance and rejection criteria are now stated explicitly.

**Changes in the manuscript:** Introduction, pages 1--3, lines 8--74;
root-selection algorithm, pages 9--10, lines 286--335;
continued-fraction methods, pages 10--11, lines 336--364; physical
interpretation, Schwarzschild benchmark, and convergence evidence, pages
11--14, lines 365--419; residual diagnostics, pages 20--23, lines 539--599;
conclusion, pages 26--27, lines 678--708.

### Referee 1, Comment 3

> I am also not fully convinced by the claimed capability of the method for
> higher overtones. The paper mainly focuses on low-lying modes, and even the
> second overtone already shows larger discrepancies and requires some caution.
> This suggests that the present method may be more reliable for low-lying modes
> than for high-\(n\) modes. This limitation should be stated explicitly. In
> particular, recent Chen-Heun-type methods have already made substantial
> progress on the accuracy of high overtones and have obtained complete
> quasinormal-mode spectra with very high precision. In comparison, the present
> method does not seem to solve the high-overtone problem, and the authors
> should avoid giving the impression that it provides a general
> complete-spectrum method.

**Response:** We agree and have narrowed the scope explicitly. Fundamentals are
classified as robust fundamental modes, scalar \(n=1\) modes as cross-validated
low-lying first overtones at \(N=32\), second overtones as exploratory diagnostics, and higher
overtones as outside scope. These tiers appear in the abstract, introduction,
results, claim-hierarchy table, conclusion, and generated catalogue CSV. We
also cite Chen *et al.* (2025), explain that confluent-Heun analytic
continuation targets complete and highly damped type-D spectra, and state that
our workflow serves the different purpose of auditing low-lying KS branches.

**Changes in the manuscript:** Abstract, page 1, lines 1--7; Introduction,
pages 2--3, lines 59--74; catalogue discussion, pages 14--18, lines 420--485;
branch-status table and limitations, page 25, lines 633--655; conclusion,
pages 26--27, lines 678--708.

### Referee 1, Comment 4

> Finally, I think the axial-gravitational results should be presented more
> carefully. Since the axial potential is not derived from a full
> gauge-invariant perturbation theory for the KS spacetime, these results should
> not be put on exactly the same physical footing as the scalar results. In my
> view, the scalar sector is the more solid part of the manuscript, while the
> axial sector is better regarded as a phenomenological extension unless a more
> complete derivation is supplied.

**Response:** We agreed with the referee's conditional phrase “unless a more
complete derivation is supplied” and supplied a more complete gauge-invariant
derivation under an explicit inverse-Cowling closure. The revision now uses the
general sourced gauge-invariant odd-parity formalism for spherical
backgrounds, includes the effective-source contribution, and replaces the old
lapse-substitution potential in both numerical solvers. The master variable is
gauge invariant by construction, while the inverse-Cowling closure remains an
additional physical assumption. It is not equated with a unique KS
quantum-source perturbation theory. We separately state the dynamical
closure—zero invariant axial effective-source current—so the reader can see
exactly what follows mathematically and what remains an assumption about the
unprovided quantum nonspherical sector. This inverse-Cowling closure is our
added model assumption. The scalar result remains the least
assumption-dependent headline, while the axial frequencies are conditional
odd-parity modes of this explicit model.

All axial frequencies, tables, figures, deformation shifts, stability diagnostics, and external comparisons in the revised manuscript were regenerated using the inverse-Cowling effective-source potential. No accepted numerical results from the earlier lapse-substitution ansatz remain.

**Changes in the manuscript:** Abstract, page 1; Section 2.1, pages 5--7,
lines 149--208; axial subsection, pages 14--16, lines 420--461; limitations,
page 25, lines 633--655; conclusion, pages 26--27, lines 696--703.

## Referee 2

### Referee 2, Comment 1

> The lapse-deformed Regge-Wheeler potential introduced in Equation (4) is a
> helpful phenomenological toy model for testing the numerical limits of your
> Chebyshev-Leaver workflow. However, as correctly noted later in the
> discussion, it is not derived from first-principles gauge-invariant
> perturbations of the true quantum-corrected field equations. In a full
> effective field theory framework, modifications to the spacetime geometry
> usually alter the linearized field equations themselves, introducing extra
> dynamic terms beyond a simple swap of the lapse function. To ensure physical
> clarity for the reader, please add a brief, explicit disclaimer directly
> below Equation (4) in Section 2.1 stating that this sector is intentionally
> phenomenological and serves primarily as a numerical benchmark. This will
> properly align the physical claims with the excellent numerical strengths of
> the paper.

**Response:** We addressed the underlying concern more fully than the requested
disclaimer: the lapse-substitution toy model has been removed as the working
axial model and is retained only as the explicitly labelled comparison
potential in Equation (17). Revised Section 2.1 starts from the general sourced
gauge-invariant odd-parity equation, imposes the additional inverse-Cowling
closure \(L_a=L=0\), derives the potential in Equation (16), and immediately
states that the resulting axial frequencies are conditional model results
rather than unique predictions of a nonspherical KS quantum theory.

**Changes in the manuscript:** Section 2.1, pages 5--7, lines 149--208.

### Referee 2, Comment 2

> In Table 2, the relative error increases for the higher overtones, reaching
> \(1.833\times10^{-5}\) for the \(\ell=4,n=2\) mode. Please clarify whether
> these higher-overtone results are intended primarily as numerical diagnostics
> or should be interpreted as quantitative physical predictions with the same
> level of confidence as the fundamental modes.

**Response:** They are not assigned equal confidence. The revision labels
\(n=0\) as robust quantitative, \(n=1\) as validated low lying, and \(n=2\)
as exploratory/diagnostic. We also separated two different error notions. The
former \(1.833\times10^{-5}\) number compared the computed Schwarzschild axial
root with an unevenly rounded transcribed reference, so it was not a clean
solver-error estimate. We removed that literature-relative-error column.
Spectral--Leaver disagreement and the uniform high-precision external
Batic--Dutykh--Sukaiti comparison are reported separately.

**Changes in the manuscript:** Axial continued-fraction table and caption,
pages 14--16, lines 420--461; branch-status statement, page 16, lines 457--461;
interpretive-status table, page 25, lines 633--655.

### Referee 2, Comment 3

> In Section 3.4, the author mentions that the largest endpoint fractional
> frequency shift observed across the full fixed-mass catalog is 7.23%,
> occurring specifically for the scalar \(\ell=4\) fundamental branch. Since
> larger multipole numbers probe different regions of the effective potential
> barrier compared to the dominant \(\ell=2\) mode, this maximal shift is an
> interesting physical result. The author should include a brief physical
> intuition explaining why the higher angular number shows a heightened
> sensitivity to the parameter \(a/M\) under a fixed mass scale.

**Response:** We first audited the percentages and now distinguish real-part
shift, damping-magnitude shift, quality-factor shift, and complex-plane
displacement. The quoted 7.23% is the \(\ell=4\) fundamental complex-plane
displacement; its real-part shift is \(-7.27\%\). We then computed
\(r_{\mathrm{peak}}\), \(V_{\mathrm{peak}}\), and the tortoise-coordinate curvature for
\(\ell=2,3,4\). At \(a/M=1\), the WKB height proxies change by
\(-7.01\%\), \(-7.13\%\), and \(-7.18\%\), respectively, with a similarly
modest multipole dependence in curvature. The revised explanation therefore
uses cautious barrier intuition and does not claim a qualitatively distinct
\(\ell=4\) mechanism. Because the axial equation itself has now been replaced
by the derived matter-corrected potential, the revised full-catalogue maximum
is the closure-dependent axial \(\ell=2\) fundamental at \(7.37\%\); the
\(7.23\%\) statement is retained only as the scalar-fundamental comparison the
referee asked us to explain.

**Changes in the manuscript:** Catalogue trends and new potential-peak table,
pages 16--18, lines 462--485.

### Referee 2, Comment 4

> In Section 3.5, the text highlights a 4.83% decrease in the quality factor
> \(Q=\operatorname{Re}(\omega)/(2|-\operatorname{Im}(\omega)|)\) as a key
> physical insight for the scalar fundamental branch. To complement this
> discussion, the author may consider adding a small inset or an extra panel to
> Figure 3 or 4 that explicitly plots \(Q\) versus \(a/M\). This will provide a
> direct visual complement for the spectroscopic trends described in the text.

**Response:** Figure 5 now contains a dedicated right panel showing
\(Q=\operatorname{Re}\omega/[2(-\operatorname{Im}\omega)]\) versus \(a/M\)
for the robust scalar \(\ell=2,n=0\) branch. The definition is identical in the
abstract, text, and caption.

**Changes in the manuscript:** Spectroscopy subsection and Figure 5,
pages 18--19, lines 486--510.

### Referee 2, Comment 5

> In Section 3.6, please clarify the physical significance of the endpoint
> \(a/M=1\). Specifically, state whether it corresponds to a physical boundary
> of the Kazakov-Solodukhin black hole (e.g., extremality or horizon
> degeneracy) or is simply a numerical cutoff adopted for the parameter scan.

**Response:** The revised text derives
\(r_h^2=4M^2+a^2\), with \(r_h>a\), and explicitly states that \(a/M=1\) is neither
extremal nor a degenerate-horizon limit. It is the upper endpoint of the
moderate deformation interval selected for this numerical study. We also
clarify the conversion to the horizon-normalized convention:
\(a/M=1\) corresponds to \((a/r_h)^2=1/5\), not \(a/r_h=1\).

**Changes in the manuscript:** Section 2.1, page 4, lines 120--130;
literature-normalization discussion, pages 19--20, lines 511--535.

### Referee 2, Comment 6

> To provide a more comprehensive background, please add a few brief sentences
> in the introduction summarizing other standard numerical methods used for
> calculating quasinormal modes such as the WKB approximation, Asymptotic
> Iteration Method (AIM), Frobenius method or direct time-domain integration
> etc. The discussion may be supported by citing representative references such
> as Rev. Mod. Phys. 83 (2011) 793, Phys. Rev. D 68, 024018 (2003), Eur. Phys.
> J. C 84, 1245 (2024), Class. Quantum Gravity 27, 155004 (2010), Eur. Phys. J.
> C 85 (2025) 1223, Phys. Rev. D 30, 295 (1984).

**Response:** We added a focused numerical-method paragraph covering WKB,
direct time-domain evolution, Frobenius/Leaver continued fractions, AIM,
spectral/pseudospectral eigenvalue methods, and recent confluent-Heun
high-overtone work. Representative references were checked and added without
turning the introduction into a general review.

**Changes in the manuscript:** Introduction, pages 1--3, lines 22--74.

### Referee 2, Comment 7

> In the abstract, several variables and parameters are rendered as normal prose
> text rather than in standard math font. Please ensure that all mathematical
> symbols throughout the abstract are consistently formatted in math font to
> maintain structural uniformity with the main text.

**Response:** The abstract was rewritten after the scientific revisions. Every
symbol and numerical expression, including \(N\), \(\ell\), \(n\), \(a/M\),
\(Q\), and scientific notation, is now in math mode. The revised abstract also
states the low-lying confidence tiers and the axial caveat.

**Changes in the manuscript:** Abstract, page 1, lines 1--7.

### Referee 2, Comment 8

> Please use “Figure 2” instead of “figure 2” (apply this consistently
> throughout the manuscript).

**Response:** All prose figure references now use “Figure.” We also audited
table and equation references, \(\ell\) notation, and the hyphenation of
“gauge-invariant,” “continued-fraction,” “high-overtone,” and
“lapse-deformed.”

**Changes in the manuscript:** Manuscript-wide editorial correction.

## Verification and reproducibility

All generated tables and figures were regenerated from the revised code. The
candidate-selection subsection now records the physical frequency window,
non-finite-root removal, \(10^{-7}\) clustering, row/column equilibration,
\(0.20\) continuation gate, competing-match score, L-BFGS-B settings, and the
backward-error and continued-fraction acceptance gates. Figure legends use
publication terminology, including “axial inverse-Cowling.” The scalar
pseudospectral endpoint diagnostic was rerun for grids \(81^2,121^2,161^2\),
half-widths \(0.020,0.025,0.030\), and \(N=32,48,64\); its endpoint change is
positive in every test and spans \(0.1002\)--\(0.1647\), so the abstract now
states the robust sign rather than a single window-dependent decimal.
The manuscript identifies OpenAI Codex, model identifier `gpt-5.6-sol`, and
lists the tasks for which it was used, together with an explicit
author-responsibility statement.

Every documented generator and test command was also run successfully from a
fresh sparse clone with no pre-existing `outputs/` directory. This audit found
and corrected one bootstrap issue in which the main pipeline assumed that the
parent output directory already existed. The fast test suite passes 8/8 checks, including a symbolic/numerical audit of the
new axial potential, and the complete validation suite passes,
with the worst spectral--Leaver catalogue difference
\(1.202\times10^{-5}\) on the deliberately exploratory axial
\(\ell=2,n=2,a/M=0.5\) row. The clean 30-page manuscript and the marked
31-page manuscript compile independently.
