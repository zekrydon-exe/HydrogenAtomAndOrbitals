"""Radial grid and finite-difference Hamiltonian for the hydrogen atom."""

from numbers import Integral, Real

import numpy as np
from scipy.sparse import diags


def _validate_integer(value, name: str, minimum: int) -> None:
    """Require a non-boolean integer greater than or equal to minimum."""
    if (
        isinstance(value, (bool, np.bool_))
        or not isinstance(value, Integral)
        or value < minimum
    ):
        raise ValueError(f"{name} must be an integer >= {minimum}.")


def _validate_positive_finite_real(value, name: str) -> None:
    """Require a finite, positive real number."""
    if (
        isinstance(value, (bool, np.bool_))
        or not isinstance(value, Real)
        or not np.isfinite(value)
        or value <= 0
    ):
        raise ValueError(f"{name} must be a finite, positive real number.")


def create_radial_grid(
    r_max: float, num_points: int
) -> tuple[np.ndarray, float]:
    """Create an evenly spaced radial grid excluding the origin.

    Parameters
    ----------
    r_max : float
        Outer radial boundary in Bohr radii.
    num_points : int
        Number of interior grid points.

    Returns
    -------
    r : numpy.ndarray
        Radial grid points, from dr to r_max - dr.
    dr : float
        Grid spacing in Bohr radii.
    """
    _validate_positive_finite_real(r_max, "r_max")
    _validate_integer(num_points, "num_points", 2)

    dr = r_max / (num_points + 1)
    r = dr * np.arange(1, num_points + 1)

    return r, dr


def effective_potential(r: np.ndarray, l: int) -> np.ndarray:
    """Calculate the hydrogen atom effective potential in atomic units.

    V_eff(r) = -1/r + l(l+1)/(2r^2)
    """
    _validate_integer(l, "l", 0)
    r = np.asarray(r, dtype=float)

    if r.ndim != 1 or r.size == 0:
        raise ValueError("r must be a non-empty one-dimensional array.")
    if not np.all(np.isfinite(r)):
        raise ValueError("All radial coordinates must be finite.")
    if np.any(r <= 0):
        raise ValueError("All radial coordinates must be positive.")

    return -1.0 / r + l * (l + 1) / (2.0 * r**2)


def build_radial_hamiltonian(r: np.ndarray, dr: float, l: int):
    """Build the sparse Hamiltonian for the reduced radial function u(r).

    The equation is

        [-1/2 d²/dr² + V_eff(r)] u(r) = E u(r)

    using a central finite-difference approximation and atomic units.
    Dirichlet boundary conditions u(0) = u(r_max) = 0 are assumed.
    The supplied grid contains only interior points.
    """
    _validate_positive_finite_real(dr, "dr")
    _validate_integer(l, "l", 0)

    r = np.asarray(r, dtype=float)

    if r.ndim != 1 or r.size < 2:
        raise ValueError("r must be a one-dimensional grid with at least 2 points.")
    if not np.all(np.isfinite(r)):
        raise ValueError("All radial coordinates must be finite.")
    if np.any(r <= 0):
        raise ValueError("All radial coordinates must be positive.")

    expected_grid = dr * np.arange(1, r.size + 1)
    if not np.allclose(r, expected_grid, rtol=1e-7, atol=1e-12):
        raise ValueError(
            "r must be the uniform interior grid dr*[1, 2, ..., N]."
        )

    diagonal = 1.0 / dr**2 + effective_potential(r, l)
    off_diagonal = np.full(r.size - 1, -0.5 / dr**2)

    return diags(
        diagonals=[off_diagonal, diagonal, off_diagonal],
        offsets=[-1, 0, 1],
        format="csr",
    )


def count_radial_nodes(
    wavefunction: np.ndarray,
    relative_threshold: float = 1e-7,
) -> int:
    """Count sign changes in a reduced radial wavefunction.

    Values smaller than the relative threshold are excluded to reduce
    false node counts caused by numerical noise in the tails.
    """
    wavefunction = np.asarray(wavefunction, dtype=float)

    if wavefunction.ndim != 1 or wavefunction.size < 2:
        raise ValueError("wavefunction must be a 1D array with at least 2 points.")
    if not np.all(np.isfinite(wavefunction)):
        raise ValueError("wavefunction must contain only finite values.")

    if (
        isinstance(relative_threshold, (bool, np.bool_))
        or not isinstance(relative_threshold, Real)
        or not np.isfinite(relative_threshold)
        or not 0 <= relative_threshold < 1
    ):
        raise ValueError("relative_threshold must be finite and in [0, 1).")

    max_amplitude = np.max(np.abs(wavefunction))

    if max_amplitude == 0:
        return 0

    significant = (
        np.abs(wavefunction) > relative_threshold * max_amplitude
    )
    signs = np.sign(wavefunction[significant])

    return int(np.count_nonzero(signs[1:] != signs[:-1]))