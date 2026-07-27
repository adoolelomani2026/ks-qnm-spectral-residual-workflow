#!/usr/bin/env python3
"""Symbolically audit the frozen-source gauge-invariant KS axial potential."""

from __future__ import annotations

from pathlib import Path

import sympy as sp


def main() -> None:
    r, a, mass, ell = sp.symbols("r a M ell", positive=True)
    sqrt_term = sp.sqrt(r**2 - a**2)
    f = (sqrt_term - 2 * mass) / r
    mass_function = sp.simplify(r * (1 - f) / 2)
    mass_prime = sp.simplify(sp.diff(mass_function, r))

    # G^t_t=8 pi rho and G^r_r=-8 pi p_r for
    # ds^2=-f dt^2+f^{-1}dr^2+r^2dOmega^2.
    rho_times_4pi = sp.simplify(mass_prime / r**2)
    radial_pressure_times_4pi = sp.simplify(-rho_times_4pi)

    matter_form = (
        f
        * (
            ell * (ell + 1) / r**2
            - 6 * mass_function / r**3
            + rho_times_4pi
            - radial_pressure_times_4pi
        )
    )
    metric_form = (
        f
        * (
            ell * (ell + 1) / r**2
            + 2 * (f - 1) / r**2
            - sp.diff(f, r) / r
        )
    )
    proxy_form = f * (ell * (ell + 1) / r**2 - 6 * mass / r**3)
    explicit_ks_form = f * (
        ell * (ell + 1) / r**2
        + (2 * sqrt_term - 2 * r - 6 * mass - a**2 / sqrt_term) / r**3
    )
    proxy_difference = f * (2 * (sqrt_term - r) - a**2 / sqrt_term) / r**3

    checks = [(sp.Rational(5), sp.Rational(1), sp.Rational(1), 2),
              (sp.Rational(7), sp.Rational(2), sp.Rational(1, 2), 3)]
    for rv, av, mv, lv in checks:
        subs = {r: rv, a: av, mass: mv, ell: lv}
        assert abs(complex(sp.N((matter_form - metric_form).subs(subs), 30))) < 1e-28
        assert abs(complex(sp.N((metric_form - explicit_ks_form).subs(subs), 30))) < 1e-28
        assert abs(complex(sp.N((metric_form - proxy_form - proxy_difference).subs(subs), 30))) < 1e-28
    assert sp.simplify(metric_form.subs(a, 0) - proxy_form.subs(a, 0)) == 0

    report = [
        "# Gauge-invariant KS axial-potential audit",
        "",
        "The calculation assumes the gauge-invariant axial effective-source current is zero,",
        "which is the frozen-source closure of the leading KS nonspherical approximation.",
        "",
        f"- f(r) = `{sp.sstr(f)}`",
        f"- m(r)=r[1-f(r)]/2 = `{sp.sstr(mass_function)}`",
        f"- m'(r) = `{sp.sstr(mass_prime)}`",
        f"- 4*pi*rho = `{sp.sstr(rho_times_4pi)}`",
        f"- 4*pi*p_r = `{sp.sstr(radial_pressure_times_4pi)}`",
        "- Therefore p_r=-rho.",
        "",
        "The general odd-parity potential is",
        "`V=f[ell(ell+1)/r^2-6m/r^3+4*pi*(rho-p_r)]`.",
        "After the background Einstein identities it becomes",
        "`V=f[ell(ell+1)/r^2+2(f-1)/r^2-f'/r]`.",
        "",
        f"- Explicit KS potential: `{sp.sstr(explicit_ks_form)}`",
        f"- Difference from the old lapse-substitution proxy: `{sp.sstr(proxy_difference)}`",
        "- The difference vanishes identically at a=0, recovering Regge-Wheeler.",
    ]

    output = Path("outputs/results/gauge_invariant_axial_potential_report.md")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
