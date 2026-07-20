# Gauge-invariant KS axial-potential audit

The calculation assumes the gauge-invariant axial effective-source current is zero,
which is the frozen-source closure of the leading KS nonspherical approximation.

- f(r) = `(-2*M + sqrt(-a**2 + r**2))/r`
- m(r)=r[1-f(r)]/2 = `M + r/2 - sqrt(-a**2 + r**2)/2`
- m'(r) = `-r/(2*sqrt(-a**2 + r**2)) + 1/2`
- 4*pi*rho = `-1/(2*r*sqrt(-a**2 + r**2)) + 1/(2*r**2)`
- 4*pi*p_r = `1/(2*r*sqrt(-a**2 + r**2)) - 1/(2*r**2)`
- Therefore p_r=-rho.

The general odd-parity potential is
`V=f[ell(ell+1)/r^2-6m/r^3+4*pi*(rho-p_r)]`.
After the background Einstein identities it becomes
`V=f[ell(ell+1)/r^2+2(f-1)/r^2-f'/r]`.

- Explicit KS potential: `(-2*M + sqrt(-a**2 + r**2))*(ell*(ell + 1)/r**2 + (-6*M - a**2/sqrt(-a**2 + r**2) - 2*r + 2*sqrt(-a**2 + r**2))/r**3)/r`
- Difference from the old lapse-substitution proxy: `(-2*M + sqrt(-a**2 + r**2))*(-a**2/sqrt(-a**2 + r**2) - 2*r + 2*sqrt(-a**2 + r**2))/r**4`
- The difference vanishes identically at a=0, recovering Regge-Wheeler.
