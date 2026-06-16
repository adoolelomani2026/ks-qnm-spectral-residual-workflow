"""Shared physics constants and KS scalar-potential utilities."""

from __future__ import annotations

import math

import numpy as np


MASS = 1.0
ELL = 2
WINDOW = (60.0, 160.0)
A_VALUES = [0.0, 0.2, 0.5, 1.0]
SPECTRAL_SIZES = [16, 32, 48, 64, 96]
FINAL_SPECTRAL_N = 96
CATALOGUE_SPECTRAL_N = 32
PERTURBATION_TYPES = ["scalar", "gravitational"]
CATALOGUE_ELL_VALUES = [2, 3, 4]
CATALOGUE_OVERTONES = [0, 1, 2]

SCHWARZSCHILD_SCALAR_L2 = 0.483643872224 - 0.096758776024j
# Used only as an initial tracking target for the first overtone branch.
SCHWARZSCHILD_SCALAR_L2_OVERTONE_ESTIMATE = 0.46385058 - 0.29560394j
SCHWARZSCHILD_REFERENCES = {
    ("scalar", 2, 0): 0.483643872224 - 0.096758776024j,
    ("scalar", 2, 1): 0.463850579020 - 0.295603936988j,
    ("scalar", 2, 2): 0.430544 - 0.508558j,
    ("scalar", 3, 0): 0.675366 - 0.096500j,
    ("scalar", 3, 1): 0.660671 - 0.292285j,
    ("scalar", 3, 2): 0.633626 - 0.496008j,
    ("scalar", 4, 0): 0.867416 - 0.096392j,
    ("scalar", 4, 1): 0.855808 - 0.290876j,
    ("scalar", 4, 2): 0.833692 - 0.490325j,
    ("gravitational", 2, 0): 0.37367168 - 0.08896232j,
    ("gravitational", 2, 1): 0.346710996 - 0.273914876j,
    ("gravitational", 2, 2): 0.301053454 - 0.478276983j,
    ("gravitational", 3, 0): 0.599443288 - 0.092703048j,
    ("gravitational", 3, 1): 0.582643804 - 0.281298113j,
    ("gravitational", 3, 2): 0.551684837 - 0.479092791j,
    ("gravitational", 4, 0): 0.809178378 - 0.094163961j,
    ("gravitational", 4, 1): 0.796631 - 0.284334j,
    ("gravitational", 4, 2): 0.772695 - 0.479900j,
}


def f_ks(r: np.ndarray | float, a: float, mass: float = MASS) -> np.ndarray | float:
    r_arr = np.asarray(r)
    return (np.sqrt(r_arr * r_arr - a * a) - 2.0 * mass) / r_arr


def df_ks(r: np.ndarray | float, a: float, mass: float = MASS) -> np.ndarray | float:
    r_arr = np.asarray(r)
    s = np.sqrt(r_arr * r_arr - a * a)
    return a * a / (s * r_arr * r_arr) + 2.0 * mass / (r_arr * r_arr)


def scalar_potential(r: np.ndarray, ell: int, a: float, mass: float = MASS) -> np.ndarray:
    f = f_ks(r, a, mass)
    return f * (ell * (ell + 1.0) / (r * r) + df_ks(r, a, mass) / r)


def regge_wheeler_potential(r: np.ndarray, ell: int, a: float, mass: float = MASS) -> np.ndarray:
    """Axial Regge-Wheeler-type potential using the KS lapse.

    For a=0 this is the standard Schwarzschild axial gravitational potential.
    For a>0 it is used as the KS lapse-deformed Regge-Wheeler validation model.
    """

    f = f_ks(r, a, mass)
    return f * (ell * (ell + 1.0) / (r * r) - 6.0 * mass / (r * r * r))


def perturbation_potential(
    r: np.ndarray,
    ell: int,
    a: float,
    perturbation_type: str = "scalar",
    mass: float = MASS,
) -> np.ndarray:
    if perturbation_type == "scalar":
        return scalar_potential(r, ell, a, mass)
    if perturbation_type == "gravitational":
        return regge_wheeler_potential(r, ell, a, mass)
    raise ValueError(f"Unknown perturbation type: {perturbation_type}")


def horizon_radius(a: float, mass: float = MASS) -> float:
    return math.sqrt((2.0 * mass) ** 2 + a * a)


def select_physical_mode(values, target: complex, exclude: list[complex] | None = None) -> complex:
    exclude = exclude or []
    candidates = [
        complex(value)
        for value in values
        if np.isfinite(value)
        and value.real > 0.05
        and value.imag < -0.02
        and value.real < 2.0
        and value.imag > -3.0
        and all(abs(value - old) > 1.0e-7 for old in exclude)
    ]
    if not candidates:
        finite = [complex(value) for value in values if np.isfinite(value)]
        if not finite:
            raise RuntimeError("No finite eigenvalues found.")
        candidates = finite
    return min(candidates, key=lambda value: abs(value - target))
