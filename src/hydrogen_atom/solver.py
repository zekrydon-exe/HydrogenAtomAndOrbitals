"""Numerical eigenvalue solver for the hydrogen atom."""

import numpy as np
from scipy.sparse.linalg import eigsh

from src.hydrogen_atom.radial import (
    build_radial_hamiltonian,
    create_radial_grid,
)


def solve_radial_schrodinger(
    r_max: float = 50.0,
    num_points: int = 1000,
    l: int = 0,
    num_states: int = 1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Solve the radial Schrödinger equation for the hydrogen atom.

    Parameters
    ----------
    r_max : float
        Outer boundary in Bohr radii.
    num_points : int
        Number of interior radial grid points.
    l : int
        Orbital angular momentum quantum number.
    num_states : int
        Number of lowest-energy states to compute.

    Returns
    -------
    r : numpy.ndarray
        Radial grid in Bohr radii.
    energies : numpy.ndarray
        Eigenvalues in Hartree, sorted in ascending order.
    wavefunctions : numpy.ndarray
        Reduced radial eigenfunctions u(r), normalised so that
        integral |u(r)|² dr = 1. Each row corresponds to one state.
    dr : float
        Radial grid spacing in Bohr radii.
    """
    if not isinstance(num_states, (int, np.integer)) or num_states < 1:
        raise ValueError("num_states must be a positive integer.")
    if num_states >= num_points:
        raise ValueError("num_states must be smaller than num_points.")

    r, dr = create_radial_grid(r_max, num_points)
    hamiltonian = build_radial_hamiltonian(r, dr, l)

    energies, eigenvectors = eigsh(
        hamiltonian,
        k=num_states,
        which="SA",
    )

    order = np.argsort(energies)
    energies = energies[order]
    eigenvectors = eigenvectors[:, order]

    # eigsh returns vectors normalised in the discrete Euclidean norm.
    # Convert them to continuum normalisation: integral |u(r)|² dr = 1.
    wavefunctions = eigenvectors.T / np.sqrt(dr)

    # Fix each arbitrary eigenvector sign for reproducible output.
    for wavefunction in wavefunctions:
        if wavefunction[np.argmax(np.abs(wavefunction))] < 0:
            wavefunction *= -1

    return r, energies, wavefunctions, dr