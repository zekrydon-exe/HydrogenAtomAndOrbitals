"""Orthogonality tests for radial eigenfunctions."""

import numpy as np

from src.hydrogen_atom.solver import solve_radial_schrodinger


def test_radial_eigenfunctions_are_orthonormal():
    _, _, wavefunctions, dr = solve_radial_schrodinger(
        r_max=50.0,
        num_points=1000,
        l=0,
        num_states=4,
    )

    overlap_matrix = (wavefunctions * dr) @ wavefunctions.T

    assert np.allclose(
        overlap_matrix,
        np.eye(4),
        atol=1e-10,
    )