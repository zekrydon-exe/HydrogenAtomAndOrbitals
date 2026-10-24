"""Radial-node validation across hydrogen orbital types."""

import pytest

from src.hydrogen_atom.radial import count_radial_nodes
from src.hydrogen_atom.solver import solve_radial_schrodinger


@pytest.mark.parametrize(
    ("n", "l", "expected_nodes"),
    [
        (1, 0, 0),  # 1s
        (2, 0, 1),  # 2s
        (2, 1, 0),  # 2p
        (3, 0, 2),  # 3s
        (3, 1, 1),  # 3p
        (3, 2, 0),  # 3d
    ],
)
def test_radial_node_count(n, l, expected_nodes):
    state_index = n - l - 1

    _, _, wavefunctions, _ = solve_radial_schrodinger(
        r_max=60.0,
        num_points=2400,
        l=l,
        num_states=state_index + 1,
    )

    actual_nodes = count_radial_nodes(wavefunctions[state_index])

    assert actual_nodes == expected_nodes