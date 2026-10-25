"""Tests for the radial grid and Hamiltonian."""

import numpy as np
import pytest
from scipy.sparse import issparse

from src.hydrogen_atom.radial import (
    build_radial_hamiltonian,
    count_radial_nodes,
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

@pytest.mark.parametrize(
    "r_max",
    [0, -1, float("nan"), float("inf"), -float("inf")],
)
def test_radial_grid_rejects_invalid_boundaries(r_max):
    with pytest.raises(ValueError):
        create_radial_grid(r_max=r_max, num_points=10)


@pytest.mark.parametrize("num_points", [1, 2.5, True, "10"])
def test_radial_grid_rejects_invalid_point_counts(num_points):
    with pytest.raises(ValueError):
        create_radial_grid(r_max=10.0, num_points=num_points)


@pytest.mark.parametrize("l", [-1, 1.5, True, float("nan")])
def test_effective_potential_rejects_invalid_angular_momentum(l):
    with pytest.raises(ValueError):
        effective_potential(np.array([1.0, 2.0]), l=l)


@pytest.mark.parametrize(
    "r",
    [
        np.array([1.0, np.nan]),
        np.array([1.0, np.inf]),
        np.array([]),
        np.array([[1.0, 2.0]]),
    ],
)
def test_effective_potential_rejects_invalid_radial_coordinates(r):
    with pytest.raises(ValueError):
        effective_potential(r, l=0)


@pytest.mark.parametrize("dr", [0, -1, float("nan"), float("inf")])
def test_hamiltonian_rejects_invalid_spacing(dr):
    r = np.array([1.0, 2.0, 3.0])

    with pytest.raises(ValueError):
        build_radial_hamiltonian(r, dr=dr, l=0)


@pytest.mark.parametrize(
    "r",
    [
        np.array([1.0, np.nan, 3.0]),
        np.array([1.0, np.inf, 3.0]),
        np.array([1.0, 2.5, 3.0]),
    ],
)
def test_hamiltonian_rejects_invalid_grids(r):
    with pytest.raises(ValueError):
        build_radial_hamiltonian(r, dr=1.0, l=0)


@pytest.mark.parametrize(
    "wavefunction",
    [
        np.array([1.0, np.nan, -1.0]),
        np.array([1.0, np.inf, -1.0]),
    ],
)
def test_node_counter_rejects_nonfinite_wavefunctions(wavefunction):
    with pytest.raises(ValueError):
        count_radial_nodes(wavefunction)


@pytest.mark.parametrize(
    "threshold",
    [float("nan"), float("inf"), -0.1, 1.0, True],
)
def test_node_counter_rejects_invalid_thresholds(threshold):
    with pytest.raises(ValueError):
        count_radial_nodes(
            np.array([1.0, -1.0, 1.0]),
            relative_threshold=threshold,
        )
