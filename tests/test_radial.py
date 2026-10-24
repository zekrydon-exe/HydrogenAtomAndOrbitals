"""Tests for the radial grid and Hamiltonian."""

import numpy as np
import pytest
from scipy.sparse import issparse

from src.hydrogen_atom.radial import (
    build_radial_hamiltonian,
    create_radial_grid,
    effective_potential,
)


def test_radial_grid_spacing():
    r, dr = create_radial_grid(r_max=10.0, num_points=99)

    assert len(r) == 99
    assert dr == pytest.approx(0.1)
    assert np.allclose(np.diff(r), dr)
    assert r[0] > 0
    assert r[-1] < 10.0


def test_effective_potential_for_s_state():
    r = np.array([1.0, 2.0, 4.0])

    potential = effective_potential(r, l=0)

    assert np.allclose(potential, [-1.0, -0.5, -0.25])


def test_effective_potential_includes_centrifugal_term():
    r = np.array([1.0, 2.0])

    potential = effective_potential(r, l=1)

    expected = -1.0 / r + 1.0 / r**2

    assert np.allclose(potential, expected)


def test_hamiltonian_is_sparse_and_symmetric():
    r, dr = create_radial_grid(r_max=10.0, num_points=99)
    hamiltonian = build_radial_hamiltonian(r, dr, l=0)

    assert issparse(hamiltonian)
    assert hamiltonian.shape == (99, 99)
    assert np.allclose(
        hamiltonian.toarray(),
        hamiltonian.T.toarray(),
    )


def test_invalid_grid_is_rejected():
    with pytest.raises(ValueError):
        create_radial_grid(r_max=-1.0, num_points=99)


def test_origin_is_rejected():
    with pytest.raises(ValueError):
        effective_potential(np.array([0.0, 1.0]), l=0)