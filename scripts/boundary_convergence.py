"""Study radial-boundary convergence for the hydrogen 4s state."""

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
    # Keep dr = 0.025 Bohr approximately constant while increasing r_max.
    r_max_values = [40.0, 50.0, 70.0, 100.0]
    dr_target = 0.025

    # For the hydrogen 4s state: n = 4, l = 0.
    # The analytical energy is E_n = -1/(2n^2) Hartree.
    n = 4
    l = 0
    exact_energy = -1.0 / (2.0 * n**2)

    rows = []

    print("HYDROGEN 4s RADIAL-BOUNDARY CONVERGENCE")
    print("-" * 84)
    print(
        f"{'r_max (a0)':>12} {'Points':>9} {'dr (a0)':>12} "
        f"{'Energy (Ha)':>18} {'Abs. error (Ha)':>20}"
    )

    for r_max in r_max_values:
        # The grid excludes both endpoints, so dr = r_max/(num_points + 1).
        num_points = round(r_max / dr_target) - 1

        r, energies, wavefunctions, dr = solve_radial_schrodinger(
            r_max=r_max,
            num_points=num_points,
            l=l,
            num_states=n,
        )

        # The fourth eigenvalue is the 4s state for l = 0.
        energy = float(energies[n - 1])
        absolute_error = abs(energy - exact_energy)

        # Check the numerical eigenfunction's discrete radial normalisation.
        norm = float(np.sum(wavefunctions[n - 1] ** 2) * dr)

        print(
            f"{r_max:12.2f} {num_points:9d} {dr:12.6f} "
            f"{energy:18.10f} {absolute_error:20.8e}"
        )

        rows.append({
            "principal_quantum_number": n,
            "angular_momentum_quantum_number": l,
            "r_max_bohr": r_max,
            "num_points": num_points,
            "dr_bohr": dr,
            "energy_hartree": energy,
            "analytical_energy_hartree": exact_energy,
            "absolute_error_hartree": absolute_error,
            "radial_norm": norm,
        })

    # Save numerical results.
    csv_path = RESULTS_DIR / "boundary_convergence.csv"

    with csv_path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    # Plot energy error against the maximum radial distance.
    r_max_values = np.array([row["r_max_bohr"] for row in rows])
    errors = np.array([row["absolute_error_hartree"] for row in rows])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.semilogy(r_max_values, errors, "o-")
    ax.set_xlabel(r"Maximum radial distance, $r_{\max}$ (Bohr)")
    ax.set_ylabel("Absolute energy error (Hartree)")
    ax.set_title(r"Hydrogen 4s Energy: Radial-Boundary Convergence")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()

    figure_path = FIGURES_DIR / "boundary_convergence.png"
    fig.savefig(figure_path, dpi=300)
    plt.close(fig)

    print(f"\nSaved data: {csv_path.relative_to(ROOT)}")
    print(f"Saved plot: {figure_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
