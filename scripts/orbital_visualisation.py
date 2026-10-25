"""Generate 3D isosurfaces of hydrogen orbital probability densities."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import sph_harm_y
from skimage.measure import marching_cubes

from src.hydrogen_atom.solver import solve_radial_schrodinger

ROOT = Path(__file__).resolve().parent.parent
FIGURES_DIR = ROOT / "figures"
FIGURES_DIR.mkdir(exist_ok=True)


def calculate_orbital(n, l, m, grid_size=100, extent=12.0):
    """Calculate |psi|^2 on a Cartesian grid for a hydrogen orbital."""
    if not (n >= 1 and 0 <= l < n and -l <= m <= l):
        raise ValueError("Invalid hydrogen quantum numbers (n, l, m).")

    r_max = max(40.0, extent * 2)
    num_points = 2400
    state_index = n - l - 1

    r, energies, wavefunctions, _ = solve_radial_schrodinger(
        r_max=r_max,
        num_points=num_points,
        l=l,
        num_states=state_index + 1,
    )

    # R(r) = u(r)/r
    radial_wavefunction = wavefunctions[state_index] / r

    axis = np.linspace(-extent, extent, grid_size)
    spacing = axis[1] - axis[0]
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")

    radius = np.sqrt(x**2 + y**2 + z**2)
    theta = np.zeros_like(radius)
    phi = np.arctan2(y, x)

    nonzero = radius > 0
    theta[nonzero] = np.arccos(
        np.clip(z[nonzero] / radius[nonzero], -1.0, 1.0)
    )

    radial_grid = np.concatenate(([0.0], r))

    # R(0) is finite for s orbitals and vanishes for l > 0.
    # Estimate the s-orbital value at the origin by extrapolating
    # the first few numerical radial values.
    if l == 0:
        radial_at_origin = float(
            np.polyfit(r[:3], radial_wavefunction[:3], deg=2)[-1]
        )
    else:
        radial_at_origin = 0.0

    radial_values = np.concatenate(
        ([radial_at_origin], radial_wavefunction)
    )

    R = np.interp(
        radius.ravel(),
        radial_grid,
        radial_values,
        left=0.0,
        right=0.0,
    ).reshape(radius.shape)

    Y = sph_harm_y(l, m, theta, phi)
    density = np.abs(R * Y) ** 2

    return x, y, z, density, energies[state_index], spacing


def plot_orbital(n, l, m, grid_size=100, extent=12.0):
    """Extract and plot a probability-density isosurface."""
    x, y, z, density, energy, spacing = calculate_orbital(
        n, l, m, grid_size=grid_size, extent=extent
    )

    maximum = float(np.max(density))
    if not np.isfinite(maximum) or maximum <= 0:
        raise RuntimeError(
            f"Invalid probability density for orbital ({n}, {l}, {m})."
        )

    # A fraction of the maximum density defines the isosurface.
    level = 0.15 * maximum

    vertices, faces, _, _ = marching_cubes(
        density,
        level=level,
        spacing=(spacing, spacing, spacing),
    )

    # marching_cubes returns coordinates relative to the grid origin.
    vertices += np.array([-extent, -extent, -extent])

    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")

    ax.plot_trisurf(
        vertices[:, 0],
        vertices[:, 1],
        vertices[:, 2],
        triangles=faces,
        linewidth=0.05,
        antialiased=True,
        alpha=0.85,
    )

    ax.set_xlim(-extent, extent)
    ax.set_ylim(-extent, extent)
    ax.set_zlim(-extent, extent)
    ax.set_box_aspect((1, 1, 1))

    ax.set_xlabel("x (Bohr radii)")
    ax.set_ylabel("y (Bohr radii)")
    ax.set_zlabel("z (Bohr radii)")
    ax.set_title(
        f"Hydrogen orbital (n={n}, l={l}, m={m})\n"
        f"Energy = {energy:.6f} Hartree"
    )

    output = FIGURES_DIR / f"orbital_n{n}_l{l}_m{m}.png"
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)

    print(
        f"Saved {output.relative_to(ROOT)} "
        f"(isosurface level={level:.3e}, vertices={len(vertices)})"
    )


def main():
    plot_orbital(n=1, l=0, m=0, grid_size=100, extent=8.0)
    plot_orbital(n=2, l=1, m=0, grid_size=100, extent=12.0)
    plot_orbital(n=3, l=2, m=0, grid_size=100, extent=15.0)


if __name__ == "__main__":
    main()