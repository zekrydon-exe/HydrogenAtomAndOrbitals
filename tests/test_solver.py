"""Tests for the hydrogen atom radial eigensolver."""

import numpy as np
import pytest

from src.hydrogen_atom.constants import GROUND_STATE_ENERGY
from src.hydrogen_atom.solver import solve_radial_schrodinger


def test_ground_state_energy():
    _, energies, _, _ = solve_radial_schrodinger(
        r_max=50.0,
        num_points=2000,
        l=0,
        num_states=1,
    )

    assert energies[0] == pytest.approx(
        GROUND_STATE_ENERGY,
        abs=2e-4,
    )


def test_wavefunction_normalisation():
    _, _, wavefunctions, dr = solve_radial_schrodinger(
        r_max=50.0,
        num_points=500,
        l=0,
        num_states=1,
    )

    norm = np.sum(np.abs(wavefunctions[0]) ** 2) * dr

    assert norm == pytest.approx(1.0, abs=1e-10)


def test_invalid_number_of_states():
    with pytest.raises(ValueError):
        solve_radial_schrodinger(
            num_points=10,
            num_states=10,
        )

@pytest.mark.parametrize(
    "kwargs",
    [
        {"r_max": 0},
        {"r_max": -1},
        {"r_max": float("nan")},
        {"r_max": float("inf")},
        {"num_points": 1},
        {"num_points": 10.5},
        {"num_points": True},
        {"l": -1},
        {"l": 1.5},
        {"l": float("nan")},
        {"num_states": 0},
        {"num_states": 10, "num_points": 10},
        {"num_states": True},
    ],
)
def test_invalid_solver_inputs(kwargs):
    arguments = {
        "r_max": 50.0,
        "num_points": 100,
        "l": 0,
        "num_states": 1,
    }
    arguments.update(kwargs)

    with pytest.raises(ValueError):
        solve_radial_schrodinger(**arguments)