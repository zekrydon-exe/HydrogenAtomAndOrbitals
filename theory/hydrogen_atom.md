# Hydrogen Atom: Theory and Numerical Method

## 1. Objective

The hydrogen atom is one of the most important exactly solvable systems in
non-relativistic quantum mechanics. It consists of an electron interacting
with a proton through the Coulomb potential.

This project solves the radial Schrödinger equation numerically and compares
the resulting energies and radial wavefunctions with analytical results.

The calculation uses atomic units and assumes a fixed nucleus. Relativistic
effects, spin-dependent interactions, and other corrections are neglected.

## 2. Atomic Units

In atomic units, the following quantities are set equal to unity:

- Reduced Planck constant: hbar = 1
- Electron mass: m_e = 1
- Elementary charge magnitude: e = 1
- Coulomb constant factor: 1/(4 pi epsilon_0) = 1

Distances are expressed in Bohr radii, and energies are expressed in Hartree.

Useful conversions are approximately:

- 1 Bohr radius = 0.529177 Angstrom
- 1 Hartree = 27.2114 electronvolts

The implementation uses the electron mass as unity and neglects the finite
proton mass. Consequently, its reference energies are the infinite-nuclear-mass
hydrogen energies.

## 3. The Schrödinger Equation

The time-independent Schrödinger equation is

    H psi = E psi

where H is the Hamiltonian, psi is the spatial wavefunction, and E is the
energy eigenvalue.

For an electron in the Coulomb field of a fixed proton, the Hamiltonian in
atomic units is

    H = -1/2 nabla^2 - 1/r

The potential energy is therefore

    V(r) = -1/r

Because this potential depends only on the distance from the nucleus, it is
spherically symmetric. The wavefunction can be separated into radial and
angular parts:

    psi(r, theta, phi) = R(r) Y_l^m(theta, phi)

Here:

- r is the distance from the nucleus.
- theta and phi are the polar and azimuthal angles.
- R(r) is the radial wavefunction.
- Y_l^m(theta, phi) is a spherical harmonic.
- l is the orbital angular momentum quantum number.
- m is the magnetic quantum number, with -l <= m <= l.

The angular functions describe the angular dependence, while the radial
equation determines the allowed radial functions and energies for a given l.

## 4. The Radial Schrödinger Equation

The radial equation for R(r) is

    -1/(2 r^2) d/dr [r^2 dR/dr]
    + l(l+1)/(2 r^2) R
    - 1/r R
    = E R

It is convenient to define the reduced radial wavefunction

    u(r) = r R(r)

Substituting this definition gives the one-dimensional reduced radial equation:

    -1/2 d^2u/dr^2 + V_eff(r) u(r) = E u(r)

where the effective potential is

    V_eff(r) = -1/r + l(l+1)/(2 r^2)

The effective potential contains two contributions:

1. Coulomb attraction: -1/r
2. Centrifugal term: l(l+1)/(2 r^2)

For l = 0, the centrifugal term vanishes. For l > 0, it contributes a
positive term that strongly affects the wavefunction near the origin.

The implementation evaluates this effective potential in
`src/hydrogen_atom/radial.py`.

## 5. Why Use the Reduced Radial Wavefunction?

The ordinary radial wavefunction R(r) and reduced radial wavefunction u(r)
are related by

    u(r) = r R(r)

Their normalisation conditions are equivalent:

    integral from 0 to infinity of |R(r)|^2 r^2 dr = 1

    integral from 0 to infinity of |u(r)|^2 dr = 1

The factor r^2 in the first expression comes from the volume element in
spherical coordinates.

For a spherically symmetric probability distribution, the radial probability
density is

    P(r) = |u(r)|^2 = r^2 |R(r)|^2

Thus, P(r) dr gives the probability of finding the electron between r and
r + dr after integrating over all angles.

The project solves for u(r), whose equation has the form of a one-dimensional
Schrödinger eigenvalue problem. This makes a finite-difference discretisation
straightforward.

The reduced radial function also has a useful boundary condition at the
origin: for the physical bound states considered here, u(0) = 0.

## 6. Boundary Conditions

For bound states, the radial wavefunction must remain physically acceptable.

At the origin:

    u(0) = 0

At large distance:

    u(r) -> 0 as r -> infinity

A numerical calculation cannot extend to infinity. Instead, the radial
domain is truncated at a finite maximum distance r_max, and the calculation
uses a finite grid of interior points.

The wavefunction is taken to vanish at the grid endpoints. Increasing r_max
allows more of the outer wavefunction to be represented, but requires more
grid points if the grid spacing is to remain fixed.

A boundary that is too close can distort the energy, especially for excited
states whose radial distributions extend farther from the nucleus.

The project investigates this effect in
`scripts/boundary_convergence.py`.

## 7. Analytical Hydrogen Energy Levels

For a hydrogen atom with a fixed nucleus, the bound-state energies are

    E_n = -1/(2 n^2) Hartree

where n is the principal quantum number:

    n = 1, 2, 3, ...

The energy depends on n but not on l or m in this model.

The allowed quantum numbers satisfy

    l = 0, 1, ..., n - 1

    m = -l, ..., 0, ..., +l

The conventional orbital labels include:

- l = 0: s
- l = 1: p
- l = 2: d
- l = 3: f

Examples of analytical energies are:

| State | n | l | Energy (Hartree) |
|---|---:|---:|---:|
| 1s | 1 | 0 | -0.500000 |
| 2s | 2 | 0 | -0.125000 |
| 2p | 2 | 1 | -0.125000 |
| 3s | 3 | 0 | -0.055556 |
| 3p | 3 | 1 | -0.055556 |
| 3d | 3 | 2 | -0.055556 |
| 4s | 4 | 0 | -0.031250 |

The equality of energies for different l values with the same n is a
characteristic degeneracy of the ideal Coulomb problem.

A numerical grid can introduce small differences between these energies.
Such differences are discretisation errors, not physical splitting in the
model used here.

## 8. Radial Nodes

A radial node is a zero of the radial wavefunction at a positive radius,
excluding the origin.

For a hydrogenic bound state, the number of radial nodes is

    number of radial nodes = n - l - 1

Examples:

| State | Expected radial nodes |
|---|---:|
| 1s | 0 |
| 2s | 1 |
| 2p | 0 |
| 3s | 2 |
| 3p | 1 |
| 3d | 0 |

Node counting provides a useful check that a computed eigenvector corresponds
to the intended radial state.

For a fixed l, the eigenvalues are ordered by increasing energy. The state
with zero-based eigenvector index j corresponds to

    n = l + j + 1

provided the requested eigenstates have been correctly ordered.

## 9. Finite-Difference Discretisation

The reduced radial equation is

    -1/2 d^2u/dr^2 + V_eff(r) u(r) = E u(r)

Choose a uniform radial grid with spacing dr. At an interior point r_i, the
second derivative is approximated by the central finite difference:

    d^2u/dr^2 at r_i
    approximately [u_(i+1) - 2u_i + u_(i-1)] / dr^2

Substitution into the radial equation produces a matrix eigenvalue problem:

    H u = E u

The resulting tridiagonal Hamiltonian has diagonal elements

    H_(i,i) = 1/dr^2 + V_eff(r_i)

and nearest-neighbour off-diagonal elements

    H_(i,i+1) = H_(i+1,i) = -1/(2 dr^2)

The matrix is real and symmetric. Its eigenvalues approximate the allowed
energies, and its eigenvectors approximate the reduced radial wavefunctions
on the grid.

The implementation constructs this sparse matrix in
`src/hydrogen_atom/radial.py` and solves for its lowest eigenvalues in
`src/hydrogen_atom/solver.py`.

Sparse eigensolvers avoid computing the full eigenspectrum when only a few
low-energy states are required.

## 10. Numerical Normalisation

The continuum normalisation condition is

    integral |u(r)|^2 dr = 1

On a uniform grid, the integral can be approximated by a numerical quadrature.
For the interior-grid eigenvectors used by the solver, the wavefunctions are
scaled using the grid spacing so that their discrete norm is approximately

    sum_i |u_i|^2 dr = 1

This normalisation is important when calculating probability densities and
comparing numerical wavefunctions with analytical solutions.

The exact accuracy of the integral depends on the discretisation and
quadrature convention.

## 11. Convergence and Error Analysis

Two distinct numerical effects should be examined.

### 11.1 Grid-spacing convergence

The grid spacing is reduced while the radial boundary is held fixed.
For a central finite-difference approximation to the second derivative,
the leading discretisation error is generally second order in dr.

In an asymptotic convergence regime, this means the energy error is expected
to scale approximately as

    error proportional to dr^2

Halving dr should therefore reduce the leading error by approximately a
factor of four, provided boundary error and other effects are sufficiently
small.

The project studies this behaviour in
`scripts/convergence_study.py`.

### 11.2 Radial-boundary convergence

The maximum radial distance r_max is increased while the grid spacing is
kept approximately constant.

This tests how strongly the finite outer boundary affects a computed state.
Boundary and discretisation errors can coexist, so increasing r_max does not
guarantee a perfectly monotonic decrease in total energy error.

The project studies this behaviour in
`scripts/boundary_convergence.py`.

### 11.3 Comparison with analytical energies

For each state, compare the numerical energy E_num with the exact energy
E_exact using the absolute error

    absolute error = |E_num - E_exact|

The relative percentage error, when E_exact is nonzero, is

    relative error (%) = 100 |E_num - E_exact| / |E_exact|

Energy agreement alone is not a complete validation of a wavefunction.
Additional checks include normalisation, orthogonality, radial-node counts,
and comparison with analytical radial functions.

## 12. Orbital Visualisation

The full spatial wavefunction is assembled from its radial and angular parts:

    psi(r, theta, phi) = R(r) Y_l^m(theta, phi)

Since R(r) = u(r)/r, the radial eigenfunction computed by the solver can be
converted to R(r) away from the origin. The angular dependence is represented
using spherical harmonics.

The probability density is

    rho(r, theta, phi) = |psi(r, theta, phi)|^2

Orbital plots can visualise the spatial distribution of this density.
The displayed shape depends on the selected quantum numbers and on the
density threshold used to construct the surface.

The probability density is not the same as the radial probability density:
rho is a probability density per unit volume, whereas |u(r)|^2 describes
probability per unit radial distance after angular integration.

## 13. Model Limitations

The present model assumes:

- A non-relativistic electron.
- A fixed proton and infinite nuclear mass.
- A pure Coulomb interaction.
- No spin-orbit coupling or other relativistic corrections.
- No external electric or magnetic fields.
- A finite numerical domain and finite grid spacing.

The calculated energies should therefore be compared with the analytical
energies of the same idealised model, rather than with precision experimental
hydrogen energies that include additional physical effects.

## 14. Reproducibility

The main calculations can be run from the project root using:

```bash
python -m scripts.ground_state
python -m scripts.excited_states
python -m scripts.orbital_states
python -m scripts.convergence_study
python -m scripts.boundary_convergence
python -m scripts.orbital_visualisation