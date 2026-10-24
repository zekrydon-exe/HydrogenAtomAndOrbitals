"""Validation tests for hydrogen atom excited states."""

import numpy as np
import pytest

from src.hydrogen_atom.solver import solve_radial_schrodinger


@pytest.mark.parametrize(
    ("state_index", "n", "tolerance"),
    [
        (0, 1, 2e-4),
        (1, 2, 1e-5),
        (2, 3, 2e-6),
    ],
)
def test_s_state_energies(state_index, n, tolerance):
    """Compare the first three s-state energies with analytical values."""
    _, energies, _, _ = solve_radial_schrodinger(
        r_max=50.0,
        num_points=2000,
        l=0,
        num_states=state_index + 1,
    )

    exact_energy = -1.0 / (2.0 * n**2)

    assert energies[state_index] == pytest.approx(
        exact_energy,
        abs=tolerance,
    )


def test_s_state_radial_node_counts():
    """Check that the first three s-states have n-1 radial nodes."""
    _, _, wavefunctions, _ = solve_radial_schrodinger(
        r_max=50.0,
        num_points=2000,
        l=0,
        num_states=3,
    )

    for state_index, wavefunction in enumerate(wavefunctions):
        # Ignore extremely small tail values to avoid counting numerical noise.
        threshold = 1e-7 * np.max(np.abs(wavefunction))
        significant = np.abs(wavefunction) > threshold
        signs = np.sign(wavefunction[significant])
        node_count = np.count_nonzero(signs[1:] != signs[:-1])

        expected_nodes = state_index

        assert node_count == expected_nodes

@pytest.mark.parametrize(
    ("n", "l", "tolerance"),
    [
        (2, 1, 5e-6),   # 2p
        (3, 1, 3e-6),   # 3p
        (3, 2, 1e-6),   # 3d
    ],
)
def test_energy_for_angular_momentum_states(n, l, tolerance):
    """Check energies for hydrogen states with nonzero angular momentum."""

    state_index = n - l - 1

    _, energies, _, _ = solve_radial_schrodinger(
        r_max=60.0,
        num_points=2400,
        l=l,
        num_states=state_index + 1,
    )

    expected_energy = -1.0 / (2.0 * n**2)

    assert energies[state_index] == pytest.approx(
        expected_energy,
        abs=tolerance,
    )