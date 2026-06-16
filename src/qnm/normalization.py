"""Normalization helpers for comparing KS QNM tables.

The project catalogue is organized at fixed mass scale and reports ``a/M`` and
``M omega``. Some KS literature instead fixes the horizon radius ``r_h`` and
reports horizon-scaled quantities. These helpers keep the conversion explicit.
"""

from __future__ import annotations

import math
from typing import TypeVar


Number = TypeVar("Number", float, complex)


def _check_alpha(alpha: float) -> None:
    if alpha < 0.0:
        raise ValueError("alpha=a/M must be nonnegative.")


def _check_beta(beta: float) -> None:
    if not (0.0 <= beta < 1.0):
        raise ValueError("beta=a/r_h must satisfy 0 <= beta < 1.")


def horizon_over_mass_from_alpha(alpha: float) -> float:
    """Return ``r_h/M`` from fixed-mass deformation ``alpha=a/M``."""

    _check_alpha(alpha)
    return math.sqrt(4.0 + alpha * alpha)


def mass_over_horizon_from_alpha(alpha: float) -> float:
    """Return ``M/r_h`` from fixed-mass deformation ``alpha=a/M``."""

    return 1.0 / horizon_over_mass_from_alpha(alpha)


def horizon_beta_from_alpha(alpha: float) -> float:
    """Return ``beta=a/r_h`` from fixed-mass deformation ``alpha=a/M``."""

    _check_alpha(alpha)
    return alpha / horizon_over_mass_from_alpha(alpha)


def alpha_from_horizon_beta(beta: float) -> float:
    """Return ``alpha=a/M`` from horizon-normalized ``beta=a/r_h``."""

    _check_beta(beta)
    return 2.0 * beta / math.sqrt(1.0 - beta * beta)


def mass_over_horizon_from_beta(beta: float) -> float:
    """Return ``M/r_h`` from horizon-normalized ``beta=a/r_h``."""

    _check_beta(beta)
    return 0.5 * math.sqrt(1.0 - beta * beta)


def momega_from_rhomega(beta: float, omega_rh: Number) -> Number:
    """Convert a horizon-scaled frequency ``r_h omega`` to ``M omega``."""

    return mass_over_horizon_from_beta(beta) * omega_rh


def rhomega_from_momega(alpha: float, omega_m: Number) -> Number:
    """Convert a mass-scaled frequency ``M omega`` to ``r_h omega``."""

    return horizon_over_mass_from_alpha(alpha) * omega_m


def fixed_horizon_to_fixed_mass(beta: float, omega_rh: Number) -> tuple[float, Number]:
    """Convert ``(a/r_h, r_h omega)`` into ``(a/M, M omega)``."""

    return alpha_from_horizon_beta(beta), momega_from_rhomega(beta, omega_rh)


def fixed_mass_to_fixed_horizon(alpha: float, omega_m: Number) -> tuple[float, Number]:
    """Convert ``(a/M, M omega)`` into ``(a/r_h, r_h omega)``."""

    return horizon_beta_from_alpha(alpha), rhomega_from_momega(alpha, omega_m)
