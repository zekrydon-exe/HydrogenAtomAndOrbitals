# Hydrogen Atom and Orbitals: Numerical Solution of the Radial Schrödinger Equation

## Abstract

This project investigates the numerical solution of the hydrogen atom's radial Schrödinger equation using the finite-difference method. The radial Hamiltonian is discretised on a uniform grid, and its lowest eigenvalues and eigenvectors are calculated using sparse matrix techniques implemented in Python.

The numerical energy levels are compared with the analytical hydrogen spectrum. Additional validation includes radial-wavefunction normalisation, orthogonality, radial-node counts, grid-convergence analysis and an investigation of finite radial-boundary effects. Three-dimensional probability-density visualisations are also generated for selected hydrogen orbitals using numerical radial functions and analytical spherical harmonics.

The ground-state calculations demonstrate second-order convergence with grid refinement. A separate boundary-convergence study of the 4s state shows that the computational domain must be sufficiently large to obtain accurate excited-state energies. The project combines numerical implementation, physical validation and reproducible scientific outputs.

---
## 1. Introduction

The hydrogen atom is a fundamental system in quantum mechanics. Because it consists of a single electron interacting with a proton, its non-relativistic Schrödinger equation admits analytical solutions. These solutions provide an excellent reference for testing numerical methods used in quantum chemistry.

The hydrogen atom is particularly suitable for computational investigation because its energy spectrum, radial wavefunctions and angular wavefunctions are known analytically. Numerical calculations can therefore be evaluated against exact reference results.

The present project focuses on solving the radial Schrödinger equation using finite differences. The resulting energy eigenvalues and radial eigenfunctions are examined to assess the accuracy, stability and convergence of the numerical implementation.

### 1.1 Objectives

The principal objectives are:

1. Formulate the hydrogen atom's radial Schrödinger equation in atomic units.

2. Construct a discrete radial Hamiltonian using central finite differences.

3. Calculate bound-state energies using a sparse eigensolver.

4. Compare numerical energy levels with analytical hydrogen energies.

5. Validate the normalisation and orthogonality of radial eigenfunctions.

6. Verify the expected radial-node structure of selected states.

7. Investigate the effect of grid refinement on numerical accuracy.

8. Investigate the effect of the outer radial boundary on excited-state energies.

9. Generate three-dimensional probability-density visualisations for selected orbitals.

10. Organise the implementation, tests, numerical results and figures into a reproducible project.

---
## 2. Theoretical Background

### 2.1 The hydrogen atom Hamiltonian

In atomic units, the non-relativistic Hamiltonian for the hydrogen atom is
$$
\hat{H}=-\frac{1}{2}\nabla^2-\frac{1}{r},
$$
where:

- $\hat{H}$ is the Hamiltonian operator;

- $\nabla^2$ is the Laplacian;

- $r$ is the electron–proton separation;

- the Coulomb potential is $-1/r$.

Atomic units simplify the equations by setting the electron mass, reduced Planck constant, elementary charge and Coulomb constant to unity under the usual atomic-unit convention.

Distances are expressed in Bohr radii, denoted $a_0$, and energies are expressed in Hartree.

### 2.2 Separation of variables

The Coulomb potential depends only on the radial coordinate. Consequently, the wavefunction can be separated into radial and angular components:
$$
\psi\_{nlm}(r,\theta,\phi)

\=

R\_{nl}(r)Y_l^m(\theta,\phi),
$$
where:

- $n$ is the principal quantum number;

- $l$ is the orbital angular momentum quantum number;

- $m$ is the magnetic quantum number;

- $R\_{nl}(r)$ is the radial wavefunction;

- $Y_l^m(\theta,\phi)$ is a spherical harmonic.

The allowed quantum numbers satisfy
$$
n=1,2,3,\ldots,
$$$$
l=0,1,\ldots,n-1,
$$$$
m=-l,-l+1,\ldots,l.
$$
The values of $l=0,1,2$ correspond to s, p and d orbitals, respectively.

### 2.3 Reduced radial wavefunction

The reduced radial wavefunction is defined as
$$
u\_{nl}(r)=rR\_{nl}(r).
$$
This transformation simplifies the radial equation by removing the first-derivative term.

The reduced radial Schrödinger equation becomes
$$
\left[

-\frac{1}{2}\frac{d^2}{dr^2}

+V\_{\mathrm{eff}}(r)

\right]u\_{nl}(r)

\=

E\_{nl}u\_{nl}(r),
$$
where the effective potential is
$$
V\_{\mathrm{eff}}(r)

\=

-\frac{1}{r}

\+

\frac{l(l+1)}{2r^2}.
$$
The effective potential contains two contributions:

1. The attractive Coulomb potential, $-1/r$.

2. The centrifugal term, $l(l+1)/(2r^2)$.

The centrifugal contribution depends on the angular momentum quantum number. It vanishes for s-states but is nonzero for p-, d- and higher-angular-momentum states.

### 2.4 Radial normalisation

The full wavefunction is normalised according to
$$
\int |\psi\_{nlm}|^2\,d\tau=1.
$$
With normalised spherical harmonics, the radial normalisation condition becomes
$$
\int_0^\infty |R\_{nl}(r)|^2r^2\,dr=1.
$$
Since $u\_{nl}(r)=rR\_{nl}(r)$, this is equivalent to
$$
\int_0^\infty |u\_{nl}(r)|^2\,dr=1.
$$
The radial probability density is therefore
$$
P\_{nl}(r)=|u\_{nl}(r)|^2.
$$
The probability of finding the electron between $r$ and $r+dr$ is approximately
$$
P\_{nl}(r)\,dr.
$$
This quantity must be distinguished from the three-dimensional probability density,
$$
\rho(r,\theta,\phi)=|\psi\_{nlm}(r,\theta,\phi)|^2.
$$
### 2.5 Analytical energy spectrum

The exact non-relativistic bound-state energies of hydrogen in atomic units are
$$
\boxed{E_n=-\frac{1}{2n^2}\ \text{Hartree}}.
$$
The analytical energy depends on the principal quantum number $n$, but not on $l$ or $m$.

Consequently, states with the same $n$ have the same analytical energy in the non-relativistic Coulomb model. For example,
$$
E\_{2s}=E\_{2p}=-0.125\ \text{Hartree},
$$
and
$$
E\_{3s}=E\_{3p}=E\_{3d}

\=-0.0555555556\ \text{Hartree}.
$$
This degeneracy provides a useful test of the numerical implementation.

### 2.6 Radial nodes

For a hydrogenic bound state, the number of radial nodes is
$$
\boxed{N\_{\mathrm{radial}}=n-l-1}.
$$
A radial node is a non-origin radial position at which the radial wavefunction vanishes and changes sign.

Examples include:

\| State | $n$ | $l$ | Expected radial nodes |

\|---|---:|---:|---:|

\| 1s | 1 | 0 | 0 |

\| 2s | 2 | 0 | 1 |

\| 2p | 2 | 1 | 0 |

\| 3s | 3 | 0 | 2 |

\| 3p | 3 | 1 | 1 |

\| 3d | 3 | 2 | 0 |

These values provide an additional physical validation criterion beyond energy agreement.

---
## 3. Numerical Method

### 3.1 Radial grid

The infinite radial domain is approximated by a finite computational interval,
$$
0\<r\<r\_{\max}.
$$
The interior radial grid contains evenly spaced points. The grid spacing is determined by the outer boundary and the number of interior points.

The origin is excluded from the interior grid because the effective potential contains singular terms at $r=0$.

### 3.2 Boundary conditions

The reduced radial wavefunction satisfies the regular-origin condition
$$
u(0)=0.
$$
At the finite outer boundary, the calculation imposes
$$
u(r\_{\max})=0.
$$
The outer boundary approximates the decay of a bound-state wavefunction at large distance. Its adequacy depends on the state being calculated and the chosen boundary radius.

Spatially extended excited states can be more sensitive to the outer boundary than the ground state.

### 3.3 Finite-difference approximation

The second derivative is approximated by the central finite-difference expression
$$
\frac{d^2u}{dr^2}

\approx

\frac{u\_{i+1}-2u_i+u\_{i-1}}{(\Delta r)^2}.
$$
Substituting this approximation into the radial Schrödinger equation gives
$$
-\frac{1}{2}

\frac{u\_{i+1}-2u_i+u\_{i-1}}{(\Delta r)^2}

\+

V\_{\mathrm{eff}}(r_i)u_i

\=

Eu_i.
$$
This can be expressed as a matrix eigenvalue problem,
$$
\mathbf{H}\mathbf{u}=E\mathbf{u},
$$
where $\mathbf{H}$ is the discrete Hamiltonian matrix.

### 3.4 Hamiltonian matrix

The finite-difference kinetic-energy operator produces a tridiagonal matrix. The potential contributes to the diagonal elements.

The resulting Hamiltonian is sparse and symmetric. These properties allow the use of efficient sparse symmetric eigensolvers rather than dense matrix diagonalisation.

The lowest eigenvalues correspond to the lowest-energy radial states for a selected value of $l$.

### 3.5 Wavefunction normalisation

The computed eigenvectors are normalised using a discrete approximation to the radial integral:
$$
\sum_i |u_i|^2\Delta r\approx 1.
$$
This ensures that the numerical wavefunctions follow the radial normalisation convention used throughout the project.

The numerical eigenvectors have arbitrary overall signs. A wavefunction and its negative represent the same physical state. When comparing numerical and analytical radial functions, their signs must therefore be aligned before calculating a difference.

---
## 4. Software Implementation

The project is implemented in Python using NumPy, SciPy, Matplotlib, pytest and scikit-image.

### 4.1 Core modules

The source package is organised under `src/hydrogen_atom/`.

\| File | Purpose |

\|---|---|

\| `constants.py` | Stores constants and analytical reference values. |

\| `radial.py` | Constructs the radial grid, effective potential and sparse Hamiltonian; includes radial-node counting. |

\| `solver.py` | Solves the radial eigenvalue problem and returns energies and radial eigenfunctions. |

\| `\_\_init\_\_.py` | Defines the Python package. |

### 4.2 Calculation scripts

The `scripts/` directory contains independent programs for the main calculations.

\| Script | Purpose |

\|---|---|

\| `ground_state.py` | Calculates the ground-state energy, evaluates its error and generates radial plots. |

\| `excited_states.py` | Calculates the lowest four s-state energies. |

\| `orbital_states.py` | Compares selected s-, p- and d-state energies with analytical values. |

\| `convergence_study.py` | Investigates ground-state energy convergence with grid refinement. |

\| `boundary_convergence.py` | Investigates 4s-state energy convergence as the radial boundary increases. |

\| `orbital_visualisation.py` | Generates three-dimensional orbital probability-density isosurfaces. |

The calculations save numerical data as CSV files and figures as PNG files under `results/` and `figures/`.

### 4.3 Automated tests

The test suite checks:

- Radial-grid construction and spacing.

- Effective-potential values for different angular momenta.

- Inclusion of the centrifugal term.

- Hamiltonian sparsity and symmetry.

- Rejection of invalid grid parameters.

- Ground-state energy accuracy.

- Wavefunction normalisation.

- Selected excited-state energy values.

- Expected radial-node counts.

- Orthogonality of radial eigenfunctions.

- Input validation for invalid parameters and non-finite values.

- Numerical-versus-analytical radial wavefunction comparisons for the 1s, 2s and 2p states.

The latest confirmed test run reported **70 passed tests**.

---
## 5. Ground-State Energy

### 5.1 Calculation parameters

The ground-state calculation used:

\| Parameter | Value |

\|---|---:|

\| Principal quantum number, $n$ | 1 |

\| Angular momentum, $l$ | 0 |

\| Outer radial boundary | $50\,a_0$ |

\| Interior grid points | 2,000 |

\| Analytical energy | $-0.5$ Hartree |

### 5.2 Numerical result

The reported numerical energy was
$$
E\_{\mathrm{numerical}}

\approx -0.4999219774\ \text{Hartree}.
$$
The analytical energy is
$$
E\_{\mathrm{exact}}=-0.5\ \text{Hartree}.
$$
The absolute energy error is
$$
\Delta E

\=

|E\_{\mathrm{numerical}}-E\_{\mathrm{exact}}|.
$$
Using the reported numerical value,
$$
\boxed{\Delta E\approx7.80\times10^{-5}\ \text{Hartree}}.
$$
The relative error is
$$
\epsilon\_{\mathrm{rel}}

\=

\frac{|E\_{\mathrm{numerical}}-E\_{\mathrm{exact}}|}

{|E\_{\mathrm{exact}}|}

\times100\%.
$$
The reported relative error is approximately
$$
\boxed{\epsilon\_{\mathrm{rel}}\approx0.0156\%}.
$$
The numerical result is close to the analytical ground-state energy.

### 5.3 Ground-state wavefunction

The ground-state script generates plots of the reduced radial wavefunction and radial probability density.

For the 1s state, the reduced radial wavefunction has no radial nodes. The radial probability density is non-negative and is proportional to $|u(r)|^2$.

The reported wavefunction norm was
$$
\sum_i |u_i|^2\Delta r\approx1.
$$
The ground-state tests also check that the numerical energy lies within a specified tolerance of the analytical result.

---
## 6. Excited-State Energy Levels

### 6.1 s-state spectrum

The lowest four s-state energies were calculated numerically and compared with the analytical hydrogen energies.

The exact reference values are
$$
E_1=-0.5,

\qquad

E_2=-0.125,
$$$$
E_3=-0.0555555556,

\qquad

E_4=-0.03125,
$$
with all energies expressed in Hartree.

The reported calculations show close agreement for the first three s-states. The 4s state displayed a larger relative error in the initial calculation, motivating a separate investigation of radial-boundary convergence.

The numerical data are saved in `results/excited_state_energies.csv`, and the associated figures are saved under `figures/`.

### 6.2 s-, p- and d-state comparisons

Selected states were calculated using different angular momentum quantum numbers.

\| State | $n$ | $l$ | Numerical energy (Hartree) | Analytical energy (Hartree) | Absolute error (Hartree) |

\|---|---:|---:|---:|---:|---:|

\| 1s | 1 | 0 | −0.4999219644 | −0.5000000000 | $7.8036\times10^{-5}$ |

\| 2s | 2 | 0 | −0.1249951216 | −0.1250000000 | $4.8784\times10^{-6}$ |

\| 2p | 2 | 1 | −0.1250016264 | −0.1250000000 | $1.6264\times10^{-6}$ |

\| 3s | 3 | 0 | −0.0555545919 | −0.0555555556 | $9.6370\times10^{-7}$ |

\| 3p | 3 | 1 | −0.0555563051 | −0.0555555556 | $7.4959\times10^{-7}$ |

\| 3d | 3 | 2 | −0.0555556198 | −0.0555555556 | $6.4242\times10^{-8}$ |

These results were obtained using a calculation with an outer boundary of $60\,a_0$ and 2,400 interior points.

The small energy differences between states with the same principal quantum number are numerical errors. The analytical model predicts exact degeneracy with respect to $l$ and $m$.

The data are saved in `results/orbital_state_energies.csv`, with the energy-level figure saved as `figures/orbital_state_energies.png`.

---
## 7. Grid-Convergence Study

### 7.1 Objective

The purpose of the grid-convergence study is to determine how the ground-state energy error changes with radial-grid resolution.

The outer boundary was held fixed at $50\,a_0$, while the number of interior grid points was increased.

This isolates the effect of grid refinement as far as possible, although the finite outer boundary remains present.

### 7.2 Results

The ground-state calculation produced the following values.

\| Grid points | Grid spacing ($a_0$) | Numerical energy (Hartree) | Absolute error (Hartree) | Error ratio |

\|---:|---:|---:|---:|---:|

\| 250 | 0.199203 | −0.4951358022 | $4.8642\times10^{-3}$ | — |

\| 500 | 0.099800 | −0.4987611470 | $1.2389\times10^{-3}$ | 3.92637 |

\| 1,000 | 0.049950 | −0.4996885125 | $3.1149\times10^{-4}$ | 3.97722 |

\| 2,000 | 0.024988 | −0.4999219774 | $7.8023\times10^{-5}$ | 3.99227 |

\| 4,000 | 0.012497 | −0.4999804800 | $1.9520\times10^{-5}$ | 3.99707 |

The analytical ground-state energy is $-0.5$ Hartree.

### 7.3 Interpretation

The central finite-difference approximation to the second derivative has a leading truncation error proportional to $(\Delta r)^2$.

In the grid-converged regime, the energy error is therefore expected to behave approximately as
$$
|E(\Delta r)-E\_{\mathrm{exact}}|

\propto(\Delta r)^2.
$$
If the grid spacing is halved, the leading error should decrease by approximately a factor of four.

The observed error ratios approach four as the grid is refined. The final reported ratio is approximately 3.997, which is consistent with second-order convergence.

This result supports the numerical consistency of the finite-difference implementation for the ground state.

The numerical results are saved in `results/ground_state_convergence.csv`, and the convergence plot is saved in `figures/ground_state_convergence.png`.

---
## 8. Radial-Boundary Convergence of the 4s State

### 8.1 Motivation

The hydrogen atom extends over an infinite radial domain, but numerical calculations use a finite outer boundary.

A boundary that is adequate for the ground state may be insufficient for an excited state because excited-state wavefunctions can extend farther from the nucleus.

The initial 4s calculation exhibited a larger relative error than the lower s-states. A boundary-convergence diagnostic was therefore performed.

### 8.2 Calculation procedure

The analytical reference energy is
$$
E_4=-\frac{1}{2(4)^2}

\=-0.03125\ \text{Hartree}.
$$
The outer boundary was varied from $40\,a_0$ to $100\,a_0$. The number of interior grid points was adjusted to maintain a grid spacing of approximately $0.025\,a_0$.

### 8.3 Results

\| Outer boundary ($a_0$) | Interior grid points | Numerical energy (Hartree) | Absolute error (Hartree) |

\|---:|---:|---:|---:|

\| 40 | 1,599 | −0.0305514661 | $6.9853\times10^{-4}$ |

\| 50 | 1,999 | −0.0312040275 | $4.5972\times10^{-5}$ |

\| 70 | 2,799 | −0.0312496415 | $3.5848\times10^{-7}$ |

\| 100 | 3,999 | −0.0312496948 | $3.0517\times10^{-7}$ |

### 8.4 Interpretation

The absolute error decreases substantially as the boundary increases from $40\,a_0$ to $70\,a_0$.

The results at $70\,a_0$ and $100\,a_0$ are close, indicating that the finite-boundary effect has become small for this state at the chosen grid resolution.

The remaining error includes discretisation error and any residual boundary effect. Consequently, the results do not imply that the numerical energy is exact or that the boundary error is identically zero.

This study demonstrates the importance of assessing the outer boundary separately from grid-spacing convergence.

---
## 9. Wavefunction Validation

Energy agreement is necessary but not sufficient to establish the correctness of a numerical eigenfunction. Additional physical properties are used to validate the computed radial wavefunctions.

### 9.1 Normalisation

The reduced radial wavefunction must satisfy
$$
\int_0^\infty |u\_{nl}(r)|^2\,dr=1.
$$
On the discrete grid, the integral is approximated by
$$
\sum_i |u_i|^2\Delta r.
$$
The automated test suite includes a normalisation test for the computed ground-state wavefunction. The reported test passed within the specified tolerance.

### 9.2 Orthogonality

For a fixed angular momentum quantum number $l$, distinct radial eigenstates should be orthogonal:
$$
\int_0^\infty

u\_{nl}(r)u\_{ml}(r)\,dr

\=

\delta\_{nm}.
$$
The numerical overlap matrix is constructed using the discrete grid integration convention:
$$
S\_{nm}

\=

\sum_i u\_{ni}u\_{mi}\Delta r.
$$
For normalised, mutually orthogonal states, the overlap matrix should approximate the identity matrix.

An automated test checks the orthonormality of the first four computed s-state eigenfunctions.

### 9.3 Radial-node counts

The expected number of radial nodes is
$$
N\_{\mathrm{radial}}=n-l-1.
$$
The implementation includes a node-counting function that detects sign changes in the reduced radial wavefunction while excluding sufficiently small tail values to reduce false counts caused by numerical noise.

Automated tests cover the 1s, 2s, 2p, 3s, 3p and 3d states.

The reported tests passed for these selected states.

Node counting is a numerical diagnostic, and its reliability depends on the grid resolution and the treatment of very small amplitudes.

### 9.4 Energy validation

The automated tests compare selected numerical energies with the analytical hydrogen spectrum.

The tested states include s-states and selected states with nonzero angular momentum.

These comparisons validate the implementation across multiple angular momentum sectors rather than relying exclusively on the ground-state result.

### 9.5 Analytical radial-wavefunction comparisons

The automated test suite compares numerical reduced radial wavefunctions with analytical reference functions for the 1s, 2s and 2p states.

Before calculating the difference, both functions are normalised on the numerical grid and their arbitrary overall signs are aligned. The comparison uses a discrete weighted L2 error, alongside an energy check, so that the assessment considers the shape of the wavefunction as well as the eigenvalue.

These tests provide a direct check of the radial eigenfunctions against known hydrogenic solutions. Extending the same comparison to additional states would provide broader coverage of excited-state behaviour.

---
## 10. Three-Dimensional Orbital Visualisation

### 10.1 Wavefunction construction

The three-dimensional orbital wavefunction is expressed as
$$
\psi\_{nlm}(r,\theta,\phi)

\=

R\_{nl}(r)Y_l^m(\theta,\phi).
$$
The radial component is obtained from the numerical solver through
$$
R\_{nl}(r)=\frac{u\_{nl}(r)}{r}.
$$
The angular component is calculated using spherical harmonics.

The probability density is
$$
\rho(r,\theta,\phi)

\=

|\psi\_{nlm}(r,\theta,\phi)|^2.
$$
### 10.2 Visualisation method

The orbital visualisation script evaluates the probability density on a three-dimensional Cartesian grid.

The density is then converted into an isosurface using the marching-cubes algorithm implemented in scikit-image.

An isosurface connects points having the same probability-density value. It provides a three-dimensional representation of the spatial distribution of the orbital density.

The current visualisation script generates figures for:

- 1s: $n=1,\ l=0,\ m=0$;

- 2p: $n=2,\ l=1,\ m=0$;

- 3d: $n=3,\ l=2,\ m=0$.

### 10.3 Expected orbital structures

The 1s probability density is spherically symmetric.

The 2p orbital with $m=0$ has angular dependence proportional to $\cos\theta$. Its probability density has two lobes along the $z$-axis, separated by a nodal plane.

The 3d orbital with $m=0$ has angular dependence proportional to
$$
3\cos^2\theta-1.
$$
Its probability density contains two lobes along the $z$-axis and an equatorial toroidal region.

The generated figures provide a spatial representation of the computed radial functions combined with the analytical angular functions.

These are probability-density visualisations, not direct plots of the complex wavefunction phase.

---
## 11. Software Testing and Reproducibility

### 11.1 Test suite

The project uses pytest to automate numerical and physical validation.

The latest confirmed complete test run reported:
$$
\boxed{70\text{ tests passed}}.
$$
The tests cover radial-grid construction, effective-potential values, Hamiltonian properties, input validation, energy comparisons, wavefunction normalisation, orthogonality, radial-node counts and numerical-versus-analytical wavefunction comparisons for the 1s, 2s and 2p states.

The tests provide regression protection against changes to the solver and radial Hamiltonian.

### 11.2 Reproducing the calculations

The calculations are organised as independent Python modules and can be executed from the project root.

The main calculation commands are:

```bash

python -m scripts.ground_state

python -m scripts.excited_states

python -m scripts.orbital_states

python -m scripts.convergence_study

python -m scripts.boundary_convergence

python -m scripts.orbital_visualisation

```

The test suite can be run using:

```bash

PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest -v

```

Plugin autoloading is disabled in this command because the development environment previously encountered an unrelated pytest plugin import problem involving ROS 2 and the missing `yaml` module.

### 11.3 Generated outputs

The project retains numerical results as CSV files and visualisations as PNG files.

Important outputs include:

- `results/ground_state_radial_data.csv`

- `results/ground_state_energy.csv`

- `results/excited_state_energies.csv`

- `results/orbital_state_energies.csv`

- `results/ground_state_convergence.csv`

- `results/boundary_convergence.csv`

- `figures/ground_state_wavefunction.png`

- `figures/ground_state_probability.png`

- `figures/excited_state_energies.png`

- `figures/excited_state_wavefunctions.png`

- `figures/orbital_state_energies.png`

- `figures/ground_state_convergence.png`

- `figures/boundary_convergence.png`

- `figures/orbital_n1_l0_m0.png`

- `figures/orbital_n2_l1_m0.png`

- `figures/orbital_n3_l2_m0.png`

The precise set of files present depends on which scripts have been executed.

---
## 12. Limitations

Although the numerical results are encouraging, several limitations remain.

### 12.1 Finite computational domain

The exact radial domain extends to infinity, whereas the numerical domain ends at a finite $r\_{\max}$.

The 4s boundary study demonstrates that insufficient radial extent can introduce substantial energy error. Boundary convergence should therefore be assessed for additional excited states before generalising the observed behaviour to the entire spectrum.

### 12.2 Finite grid spacing

The finite-difference approximation introduces discretisation error.

The ground-state convergence study demonstrates second-order convergence over the tested grid resolutions. This does not, by itself, establish the same error behaviour for every excited state or angular momentum sector.

### 12.3 Radial wavefunction comparison

The energy spectrum, normalisation, orthogonality and selected radial-node counts have been tested.

Automated tests now compare numerical radial wavefunctions against analytical solutions for the 1s, 2s and 2p states. A broader comparison with tabulated wavefunction errors across additional states remains a possible extension.

### 12.4 Three-dimensional visualisation

The orbital visualisations are generated on finite Cartesian grids and use a selected isovalue.

Their appearance therefore depends on grid resolution, spatial extent and the chosen density threshold. An isosurface does not represent a boundary beyond which the electron cannot be found; it represents a selected probability-density level.

### 12.5 Physical model

The project uses the non-relativistic hydrogen atom model. Relativistic corrections, spin–orbit interactions, finite nuclear mass corrections beyond the adopted approximation, and external-field effects are not included in the current calculation.

---
## 13. Future Work

Potential extensions of the project include:

1. Systematic grid-convergence studies for multiple excited states and angular momentum sectors.

2. Automated radial-boundary convergence tests.

3. Numerical-versus-analytical radial wavefunction comparisons for a wider range of states.

4. Additional orthogonality tests across s-, p- and d-state sectors.

5. Quantitative residual checks of the discrete eigenvalue equation.

6. Improved orbital visualisation and comparison of different magnetic quantum numbers.

7. Investigation of hydrogenic ions with nuclear charge $Z>1$.

8. Extension towards multi-electron systems and molecular electronic-structure methods.

These extensions would provide a pathway from a single-electron numerical solver towards more advanced computational quantum chemistry.

---
## 14. Conclusion

This project implements a finite-difference solution of the hydrogen atom's radial Schrödinger equation using a sparse Hamiltonian and numerical eigenvalue methods.

The calculated ground-state energy agrees closely with the analytical result. The grid-convergence study demonstrates the expected second-order error reduction as the radial grid is refined.

The computed s-, p- and d-state energies reproduce the analytical hydrogen spectrum to within small numerical errors. Additional tests examine wavefunction normalisation, orthogonality and radial-node counts.

The 4s radial-boundary study demonstrates that the computational domain is an important source of numerical error for excited states. Increasing the outer boundary substantially improves the calculated energy until the remaining error becomes small at the chosen grid resolution.

Finally, the project combines numerical radial solutions with spherical harmonics to generate three-dimensional orbital probability-density visualisations.

Overall, the work establishes a tested and reproducible computational foundation for solving central-potential Schrödinger equations and provides a basis for further development in computational quantum chemistry.