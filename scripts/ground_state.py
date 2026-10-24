"""Calculate and validate the hydrogen atom ground state."""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from src.hydrogen_atom.constants import GROUND_STATE_ENERGY
from src.hydrogen_atom.solver import solve_radial_schrodinger


# Project directories
ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "results"
FIGURES_DIR = ROOT / "figures"

RESULTS_DIR.mkdir(exist_ok=True)
FIGURES_DIR.mkdir(exist_ok=True)


def main():
    # Numerical calculation
    r, energies, wavefunctions, dr = solve_radial_schrodinger(
        r_max=50.0,
        num_points=2000,
        l=0,
        num_states=1,
    )

    energy = energies[0]
    u_numerical = wavefunctions[0]

    # Analytical hydrogen 1s reduced radial wavefunction.
    # Atomic units: u_1s(r) = 2*r*exp(-r).
    u_analytical = 2.0 * r * np.exp(-r)

    # Radial probability densities: P(r) = |u(r)|^2.
    p_numerical = np.abs(u_numerical) ** 2
    p_analytical = np.abs(u_analytical) ** 2

    # Numerical normalisation and energy error
    norm = np.sum(p_numerical) * dr
    absolute_error = abs(energy - GROUND_STATE_ENERGY)
    relative_error_percent = (
        absolute_error / abs(GROUND_STATE_ENERGY) * 100.0
    )

    print("HYDROGEN ATOM: GROUND STATE")
    print("-" * 42)
    print(f"Grid points:             {len(r)}")
    print(f"Radial boundary:         50 Bohr radii")
    print(f"Numerical energy:        {energy:.12f} Hartree")
    print(f"Analytical energy:       {GROUND_STATE_ENERGY:.12f} Hartree")
    print(f"Absolute energy error:   {absolute_error:.6e} Hartree")
    print(f"Relative energy error:   {relative_error_percent:.6f} %")
    print(f"Wavefunction norm:       {norm:.12f}")

    # Save radial data
    radial_csv = RESULTS_DIR / "ground_state_radial_data.csv"

    with radial_csv.open("w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow([
            "r_bohr",
            "u_numerical",
            "u_analytical",
            "probability_numerical",
            "probability_analytical",
        ])

        writer.writerows(zip(
            r,
            u_numerical,
            u_analytical,
            p_numerical,
            p_analytical,
        ))

    # Save energy summary
    energy_csv = RESULTS_DIR / "ground_state_energy.csv"

    with energy_csv.open("w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow([
            "energy_numerical_hartree",
            "energy_analytical_hartree",
            "absolute_error_hartree",
            "relative_error_percent",
            "wavefunction_norm",
            "num_points",
            "r_max_bohr",
        ])
        writer.writerow([
            energy,
            GROUND_STATE_ENERGY,
            absolute_error,
            relative_error_percent,
            norm,
            len(r),
            50.0,
        ])

    # Plot reduced radial wavefunction
    plt.figure(figsize=(8, 5))
    plt.plot(r, u_numerical, label="Numerical")
    plt.plot(r, u_analytical, "--", label="Analytical")
    plt.xlim(0, 10)
    plt.xlabel("Radius, r (Bohr radii)")
    plt.ylabel("Reduced radial wavefunction, u(r)")
    plt.title("Hydrogen 1s Ground-State Wavefunction")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "ground_state_wavefunction.png", dpi=300)
    plt.close()

    # Plot radial probability density
    plt.figure(figsize=(8, 5))
    plt.plot(r, p_numerical, label="Numerical")
    plt.plot(r, p_analytical, "--", label="Analytical")
    plt.xlim(0, 10)
    plt.xlabel("Radius, r (Bohr radii)")
    plt.ylabel("Radial probability density, |u(r)|²")
    plt.title("Hydrogen 1s Radial Probability Density")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "ground_state_probability.png", dpi=300)
    plt.close()

    print("\nFiles saved:")
    print(f"  {radial_csv.relative_to(ROOT)}")
    print(f"  {energy_csv.relative_to(ROOT)}")
    print("  figures/ground_state_wavefunction.png")
    print("  figures/ground_state_probability.png")


if __name__ == "__main__":
    main()