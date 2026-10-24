"""Radial grid and finite-difference Hamiltonian for the hydrogen atom."""

import numpy as np
from scipy.sparse import diags


def create_radial_grid(r_max: float, num_points: int) -> tuple[np.ndarray, float]:
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
    if r_max <= 0:
        raise ValueError("r_max must be positive.")
    if num_points < 2:
        raise ValueError("num_points must be at least 2.")

    dr = r_max / (num_points + 1)
    r = dr * np.arange(1, num_points + 1)

    return r, dr


def effective_potential(r: np.ndarray, l: int) -> np.ndarray:
    """Calculate the hydrogen atom effective potential in atomic units.

    V_eff(r) = -1/r + l(l+1)/(2r^2)

    Parameters
    ----------
    r : numpy.ndarray
        Positive radial coordinates in Bohr radii.
    l : int
        Orbital angular momentum quantum number (l >= 0).

    Returns
    -------
    numpy.ndarray
        Effective potential in Hartree.
    """
    r = np.asarray(r, dtype=float)

    if l < 0 or not isinstance(l, (int, np.integer)):
        raise ValueError("l must be a non-negative integer.")
    if np.any(r <= 0):
        raise ValueError("All radial coordinates must be positive.")

    return -1.0 / r + l * (l + 1) / (2.0 * r**2)


def build_radial_hamiltonian(
    r: np.ndarray, dr: float, l: int
):
    """Build the sparse Hamiltonian for the reduced radial function u(r).

    The equation is

        [-1/2 d²/dr² + V_eff(r)] u(r) = E u(r)

    using a central finite-difference approximation and atomic units.

    Dirichlet boundary conditions u(0) = u(r_max) = 0 are assumed.
    The supplied grid contains only interior points.

    Returns
    -------
    scipy.sparse.csr_matrix
        Real symmetric Hamiltonian matrix in Hartree.
    """
    r = np.asarray(r, dtype=float)

    if r.ndim != 1 or r.size < 2:
        raise ValueError("r must be a one-dimensional grid with at least 2 points.")
    if dr <= 0:
        raise ValueError("dr must be positive.")
    if l < 0 or not isinstance(l, (int, np.integer)):
        raise ValueError("l must be a non-negative integer.")
    if np.any(r <= 0):
        raise ValueError("All radial coordinates must be positive.")
    if not np.allclose(np.diff(r), dr):
        raise ValueError("The radial grid must be evenly spaced with spacing dr.")

    diagonal = 1.0 / dr**2 + effective_potential(r, l)
    off_diagonal = np.full(r.size - 1, -0.5 / dr**2)

    hamiltonian = diags(
        diagonals=[off_diagonal, diagonal, off_diagonal],
        offsets=[-1, 0, 1],
        format="csr",
    )

    return hamiltonian

def count_radial_nodes(
    wavefunction: np.ndarray,
    relative_threshold: float = 1e-7,
) -> int:
    """Count sign changes in a reduced radial wavefunction.

    Values smaller than the relative threshold are excluded to
    reduce false node counts caused by numerical noise in the tails.
    """
    wavefunction = np.asarray(wavefunction, dtype=float)

    if wavefunction.ndim != 1 or wavefunction.size < 2:
        raise ValueError("wavefunction must be a 1D array with at least 2 points.")

    if not 0 <= relative_threshold < 1:
        raise ValueError("relative_threshold must be in [0, 1).")

    max_amplitude = np.max(np.abs(wavefunction))

    if max_amplitude == 0:
        return 0

    significant = np.abs(wavefunction) > relative_threshold * max_amplitude
    signs = np.sign(wavefunction[significant])

    return int(np.count_nonzero(signs[1:] != signs[:-1]))