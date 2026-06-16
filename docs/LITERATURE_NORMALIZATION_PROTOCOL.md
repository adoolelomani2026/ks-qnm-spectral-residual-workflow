# Literature Normalization Protocol

This note records the conversion needed before making a quantitative comparison
between this fixed-mass catalogue and KS QNM tables that use a fixed horizon
radius.

## Geometry

The implemented KS lapse is

```text
f_a(r) = sqrt(r^2 - a^2) / r - 2M / r
```

The horizon satisfies

```text
r_h = sqrt(4 M^2 + a^2).
```

This project reports fixed-mass quantities:

```text
alpha = a / M
M omega
```

For horizon-normalized tables, use

```text
beta = a / r_h
r_h omega
```

## Parameter Conversion

From this project's fixed-mass convention:

```text
r_h / M = sqrt(4 + alpha^2)
beta = alpha / sqrt(4 + alpha^2)
```

From a horizon-normalized table:

```text
M / r_h = 0.5 * sqrt(1 - beta^2)
alpha = 2 beta / sqrt(1 - beta^2)
```

For example, this project's `a/M = 1` endpoint corresponds to

```text
beta = a / r_h = 1 / sqrt(5) = 0.4472135955...
```

not to a horizon-normalized table row with `a/r_h = 1`.

## Frequency Conversion

Convert a horizon-scaled literature frequency into this project's convention:

```text
M omega = (M / r_h) * (r_h omega)
```

Convert this project's frequency into a horizon-scaled convention:

```text
r_h omega = (r_h / M) * (M omega).
```

The helper functions live in `src/qnm/normalization.py`.

## Comparison Scaffold

Transcribe published rows into:

```text
data/literature/ks_qnm_literature_template.csv
```

Then enable those rows and run:

```text
python scripts/prepare_literature_comparison.py
```

The script writes:

```text
outputs/results/literature_normalized_comparison.csv
```

It converts horizon-scaled entries into `M omega`, matches the corresponding
catalogue branch by `(perturbation_type, ell, overtone)`, and reports whether
the converted deformation parameter is an exact catalogue point or only a
nearest-neighbor comparison. It deliberately does not interpolate; interpolation
should be added only after a denser scalar grid is generated.

## Quantitative Comparison Workflow

1. Identify the exact convention used in the published table: fixed `M`,
   fixed `r_h`, or another scale choice.
2. Restrict the first comparison to scalar modes, where this project has the
   cleanest physical interpretation.
3. Convert the published deformation parameter to `alpha = a/M`.
4. Convert the published frequency to `M omega`.
5. Run the Chebyshev-Leaver solver at the converted `alpha`, or interpolate
   only if the converted point lies inside a dense enough grid.
6. Compare like with like: perturbation sector, spin, `ell`, overtone index,
   boundary convention, and numerical method.
7. Report differences separately from the spectroscopy diagnostics, because
   WKB, time-domain, Frobenius, and spectral methods have different error
   models.

## Caveats

- The axial rows in this project use a KS-lapse-deformed Regge-Wheeler model,
  so they should not be compared as final gauge-invariant gravitational
  predictions.
- Published overtone labels can be branch-sensitive; use continuation and
  residual checks before comparing high overtones.
- A qualitative literature comparison is publication-safe now. A quantitative
  comparison should wait until the converted parameter grid and branch labels
  are explicitly matched.
