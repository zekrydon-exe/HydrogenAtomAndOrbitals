"""Calculate hydrogen energies for different angular momentum states."""

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
    # Allowed hydrogen states: (principal quantum number n, angular momentum l).
    states = [
        (1, 0),  # 1s
        (2, 0),  # 2s
        (2, 1),  # 2p
        (3, 0),  # 3s
        (3, 1),  # 3p
        (3, 2),  # 3d
    ]

    rows = []

    print("HYDROGEN ATOM: s-, p- AND d-STATES")
    print("-" * 78)
    print(
        f"{'State':>6} {'n':>3} {'l':>3} "
        f"{'Numerical (Ha)':>19} {'Analytical (Ha)':>19} {'Abs. error':>12}"
    )

    for n, l in states:
        # The first bound state for a given l has n = l + 1.
        # Its index among the l-sector eigenvalues is n - l - 1.
        state_index = n - l - 1

        _, energies, _, _ = solve_radial_schrodinger(
            r_max=60.0,
            num_points=2400,
            l=l,
            num_states=state_index + 1,
        )

        numerical_energy = energies[state_index]
        analytical_energy = -1.0 / (2.0 * n**2)
        absolute_error = abs(numerical_energy - analytical_energy)

        state_name = f"{n}{'spdf'[l]}"

        print(
            f"{state_name:>6} {n:3d} {l:3d} "
            f"{numerical_energy:19.10f} "
            f"{analytical_energy:19.10f} "
            f"{absolute_error:12.4e}"
        )

        rows.append({
            "state": state_name,
            "n": n,
            "l": l,
            "numerical_energy_hartree": numerical_energy,
            "analytical_energy_hartree": analytical_energy,
            "absolute_error_hartree": absolute_error,
        })

    # Save the comparison table.
    csv_path = RESULTS_DIR / "orbital_state_energies.csv"

    with csv_path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    # Plot the calculated energy levels.
    fig, ax = plt.subplots(figsize=(9, 6))

    for row in rows:
        n = row["n"]
        energy = row["numerical_energy_hartree"]
        ax.hlines(energy, n - 0.35, n + 0.35, linewidth=2)
        ax.text(n + 0.38, energy, row["state"], va="center")

    ax.set_xlabel("Principal quantum number, n")
    ax.set_ylabel("Energy (Hartree)")
    ax.set_title("Hydrogen Atom Energy Levels by Orbital Type")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "orbital_state_energies.png", dpi=300)
    plt.close(fig)

    print(f"\nSaved data: {csv_path.relative_to(ROOT)}")
    print("Saved plot: figures/orbital_state_energies.png")


if __name__ == "__main__":
    main()