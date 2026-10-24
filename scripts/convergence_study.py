"""Study radial-grid convergence for hydrogen atom energy levels."""

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
    # Hold r_max fixed to investigate grid-spacing effects.
    r_max = 50.0
    grid_sizes = [250, 500, 1000, 2000, 4000]

    # Ground-state analytical energy in Hartree.
    exact_energy = -0.5

    rows = []

    print("GROUND-STATE GRID CONVERGENCE")
    print("-" * 72)
    print(
        f"{'Points':>8} {'dr (a0)':>12} {'Energy (Ha)':>18} "
        f"{'Abs. error (Ha)':>18} {'Error ratio':>13}"
    )

    previous_error = None

    for num_points in grid_sizes:
        r, energies, _, dr = solve_radial_schrodinger(
            r_max=r_max,
            num_points=num_points,
            l=0,
            num_states=1,
        )

        energy = energies[0]
        error = abs(energy - exact_energy)

        error_ratio = (
            previous_error / error
            if previous_error is not None
            else float("nan")
        )

        print(
            f"{num_points:8d} {dr:12.6f} {energy:18.10f} "
            f"{error:18.8e} {error_ratio:13.5f}"
        )

        rows.append({
            "num_points": num_points,
            "r_max_bohr": r_max,
            "dr_bohr": dr,
            "energy_hartree": energy,
            "absolute_error_hartree": error,
            "error_ratio": error_ratio,
        })

        previous_error = error

    # Save results to CSV.
    csv_path = RESULTS_DIR / "ground_state_convergence.csv"

    with csv_path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    # Plot absolute error against grid spacing.
    dr_values = np.array([row["dr_bohr"] for row in rows])
    errors = np.array([row["absolute_error_hartree"] for row in rows])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.loglog(dr_values, errors, "o-", label="Numerical error")

    # Reference slope for second-order convergence.
    reference = errors[-1] * (dr_values / dr_values[-1]) ** 2
    ax.loglog(dr_values, reference, "--", label=r"Reference slope: $\Delta r^2$")

    ax.set_xlabel("Grid spacing, Δr (Bohr radii)")
    ax.set_ylabel("Absolute energy error (Hartree)")
    ax.set_title("Hydrogen Ground-State Grid Convergence")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "ground_state_convergence.png", dpi=300)
    plt.close(fig)

    print(f"\nSaved data: {csv_path.relative_to(ROOT)}")
    print("Saved plot: figures/ground_state_convergence.png")


if __name__ == "__main__":
    main()