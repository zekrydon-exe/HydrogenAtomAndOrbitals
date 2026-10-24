"""Calculate hydrogen s-state energies and compare with analytical values."""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from src.hydrogen_atom.solver import solve_radial_schrodinger


ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "results"
FIGURES_DIR = ROOT / "figures"

RESULTS_DIR.mkdir(exist_ok=True)
FIGURES_DIR.mkdir(exist_ok=True)


def main():
    # Calculate the four lowest l=0 states.
    r, energies, wavefunctions, dr = solve_radial_schrodinger(
        r_max=50.0,
        num_points=2000,
        l=0,
        num_states=4,
    )

    # For l=0, the bound states correspond to n=1, 2, 3, ...
    n_values = np.arange(1, 5)
    analytical_energies = -1.0 / (2.0 * n_values**2)

    absolute_errors = np.abs(energies - analytical_energies)
    relative_errors = (
        absolute_errors / np.abs(analytical_energies) * 100.0
    )

    print("HYDROGEN ATOM: LOWEST FOUR s-STATES")
    print("-" * 82)
    print(
        f"{'n':>3} {'Numerical (Ha)':>19} "
        f"{'Analytical (Ha)':>19} {'Abs. error (Ha)':>19} {'Rel. error (%)':>15}"
    )

    for i, n in enumerate(n_values):
        print(
            f"{n:3d} {energies[i]:19.10f} "
            f"{analytical_energies[i]:19.10f} "
            f"{absolute_errors[i]:19.6e} "
            f"{relative_errors[i]:15.6f}"
        )

    # Save energy comparison.
    csv_path = RESULTS_DIR / "excited_state_energies.csv"

    with csv_path.open("w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow([
            "n",
            "l",
            "numerical_energy_hartree",
            "analytical_energy_hartree",
            "absolute_error_hartree",
            "relative_error_percent",
        ])

        for i, n in enumerate(n_values):
            writer.writerow([
                n,
                0,
                energies[i],
                analytical_energies[i],
                absolute_errors[i],
                relative_errors[i],
            ])

    # Plot numerical and analytical energy levels.
    fig, ax = plt.subplots(figsize=(8, 5))

    for i, n in enumerate(n_values):
        ax.hlines(
            energies[i],
            n - 0.3,
            n + 0.3,
            linewidth=2,
            label="Numerical" if i == 0 else None,
        )
        ax.hlines(
            analytical_energies[i],
            n - 0.3,
            n + 0.3,
            linestyle="--",
            linewidth=2,
            label="Analytical" if i == 0 else None,
        )

    ax.set_xticks(n_values)
    ax.set_xlabel("Principal quantum number, n")
    ax.set_ylabel("Energy (Hartree)")
    ax.set_title("Hydrogen Atom: Lowest Four s-State Energies")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "excited_state_energies.png", dpi=300)
    plt.close(fig)

    # Plot reduced radial wavefunctions.
    fig, ax = plt.subplots(figsize=(9, 5))

    for i, n in enumerate(n_values):
        ax.plot(r, wavefunctions[i], label=f"Numerical {n}s")

    ax.set_xlim(0, 35)
    ax.set_xlabel("Radius, r (Bohr radii)")
    ax.set_ylabel("Reduced radial wavefunction, u(r)")
    ax.set_title("Hydrogen Atom: Lowest Four s-State Wavefunctions")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "excited_state_wavefunctions.png", dpi=300)
    plt.close(fig)

    print(f"\nSaved energy data: {csv_path.relative_to(ROOT)}")
    print("Saved figures:")
    print("  figures/excited_state_energies.png")
    print("  figures/excited_state_wavefunctions.png")


if __name__ == "__main__":
    main()