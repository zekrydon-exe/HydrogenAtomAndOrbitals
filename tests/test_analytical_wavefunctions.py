"""Compare numerical radial wavefunctions with analytical solutions."""

import numpy as np
import pytest

from src.hydrogen_atom.solver import solve_radial_schrodinger


def analytical_1s_reduced_radial(r):
    """Exact reduced radial wavefunction u_1s(r) in atomic units."""
    return 2.0 * r * np.exp(-r)


def test_ground_state_wavefunction_matches_analytical_1s():
    """The numerical 1s wavefunction should match u_1s(r) = 2r exp(-r)."""
    r, energies, wavefunctions, dr = solve_radial_schrodinger(
        r_max=50.0,
        num_points=2000,
        l=0,
        num_states=1,
    )

    numerical = wavefunctions[0]
    analytical = analytical_1s_reduced_radial(r)

    # Normalise the analytical function on the same finite numerical grid.
    analytical_norm = np.sqrt(np.sum(analytical**2) * dr)
    analytical = analytical / analytical_norm

    # Eigenfunctions have an arbitrary overall sign.
    overlap = np.sum(numerical * analytical) * dr
    if overlap < 0:
        numerical = -numerical
        overlap = -overlap

    # Compare the functions using their grid-weighted L2 difference.
    l2_error = np.sqrt(np.sum((numerical - analytical)**2) * dr)

    assert energies[0] == pytest.approx(-0.5, abs=2e-4)
    assert overlap == pytest.approx(1.0, abs=1e-5)
    assert l2_error < 5e-3

def analytical_2s_reduced_radial(r):
    """Exact reduced radial wavefunction u_2s(r) in atomic units."""
    return r * (2.0 - r) * np.exp(-r / 2.0) / (2.0 * np.sqrt(2.0))


def analytical_2p_reduced_radial(r):
    """Exact reduced radial wavefunction u_2p(r) in atomic units."""
    return r**2 * np.exp(-r / 2.0) / (2.0 * np.sqrt(6.0))


def assert_wavefunction_matches(numerical, analytical, dr, tolerance=1e-2):
    """Compare normalised radial functions, accounting for arbitrary sign."""
    analytical_norm = np.sqrt(np.sum(analytical**2) * dr)
    analytical = analytical / analytical_norm

    overlap = np.sum(numerical * analytical) * dr

    # Eigenfunctions can differ by an arbitrary overall sign.
    if overlap < 0:
        numerical = -numerical
        overlap = -overlap

    l2_error = np.sqrt(np.sum((numerical - analytical)**2) * dr)

    assert overlap > 0.99
    assert l2_error < tolerance


def test_2s_wavefunction_matches_analytical_solution():
    r, energies, wavefunctions, dr = solve_radial_schrodinger(
        r_max=50.0,
        num_points=2000,
        l=0,
        num_states=2,
    )

    assert energies[1] == pytest.approx(-0.125, abs=2e-5)

    analytical = analytical_2s_reduced_radial(r)

    assert_wavefunction_matches(
        wavefunctions[1],
        analytical,
        dr,
    )


def test_2p_wavefunction_matches_analytical_solution():
    r, energies, wavefunctions, dr = solve_radial_schrodinger(
        r_max=50.0,
        num_points=2000,
        l=1,
        num_states=1,
    )

    assert energies[0] == pytest.approx(-0.125, abs=2e-5)

    analytical = analytical_2p_reduced_radial(r)

    assert_wavefunction_matches(
        wavefunctions[0],
        analytical,
        dr,
    )
