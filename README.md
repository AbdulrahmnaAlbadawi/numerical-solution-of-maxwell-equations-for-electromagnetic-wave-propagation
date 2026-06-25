# Numerical Solution of Maxwell's Equations for Electromagnetic Wave Propagation

An introduction to time-domain methods for solving Maxwell's equations in 1D, comparing the **explicit Finite-Difference Time-Domain (FDTD)** method against the **implicit Crank-Nicolson** method. The project simulates a Gaussian pulse propagating through 1D free space and visualizes the fundamental trade-offs between the two approaches.

This was developed as a course project for **MAT232E**.

## Overview

Maxwell's equations govern electromagnetic phenomena, but most real-world problems are too complex for analytical solutions. This project demonstrates two computational engines for solving the 1D time-domain Maxwell's equations and highlights the engineering trade-offs between them.

Starting from the two coupled curl equations in a source-free, lossless medium:

```
μ₀ ∂Hy/∂t = ∂Ex/∂z
ε₀ ∂Ex/∂t = ∂Hy/∂z
```

the continuous derivatives are discretized into finite differences and solved on a computational grid.

## The Two Methods

| Feature | FDTD (Explicit) | Crank-Nicolson (Implicit) |
|---|---|---|
| Method type | Explicit | Implicit |
| Core operation | Simple "leapfrog" arithmetic updates | Solve a linear system `A·x = b` each step |
| Stability | Conditionally stable (CFL: `c·Δt ≤ Δz`) | Unconditionally stable |
| Cost per step | Low, `O(N)` | High, matrix solve |
| Accuracy | 2nd-order in space and time | 2nd-order in time |
| Key weakness | Stability constraint forces tiny time steps on fine grids | Heavier, more complex steps |
| Typical use | Broadband pulse propagation, large simple geometries | Problems needing large time steps relative to grid size |

**FDTD** uses the staggered Yee grid with leapfrog time-stepping: the future field is computed directly from past values, with no matrix inversion. Fast per step, but bound by the Courant-Friedrichs-Lewy (CFL) condition.

**Crank-Nicolson** averages the spatial derivative over the current and future time steps, producing a tridiagonal system that must be solved at every step. This breaks the CFL "speed limit" and remains stable for any `Δt` — in this simulation it uses a time step 4× larger than the FDTD limit.

## Repository Structure

```
.
├── matlab/        # MATLAB simulation source
├── docs/          # Presentation slides and supporting material
└── README.md
```

## Running the Simulation

Requires **MATLAB** (no additional toolboxes needed).

1. Open MATLAB and navigate to the `matlab/` directory.
2. Run the main script:
   ```matlab
   MAT232E_Project_CODE_G7
   ```
3. An animated figure displays both solutions side by side:
   - **Blue (solid):** FDTD, CFL-constrained
   - **Red (dashed):** Crank-Nicolson, unconditionally stable, running at 4× the time step

### Key Parameters

| Parameter | Value | Description |
|---|---|---|
| `c` | 3×10⁸ m/s | Speed of light |
| `z_limit` | 1.0 m | Domain length |
| `dz` | 0.005 m | Spatial step |
| `dt_fdtd` | `dz/c` | FDTD time step (CFL limit) |
| `dt_cn` | `4·dt_fdtd` | Crank-Nicolson time step |

The source is a Gaussian pulse injected at one-third of the domain length, with **Perfect Electric Conductor (PEC)** boundaries.

## Numerical Dispersion

A key artifact demonstrated by the simulation is **numerical dispersion**: in a vacuum all frequencies travel at exactly `c`, but on a discrete grid, frequencies with wavelengths approaching `Δz` are artificially slowed. This causes the pulse to spread and distort as it propagates. A finer grid (more points per wavelength) reduces this error.

## Key Takeaway

There is no single "best" method — only the right tool for the job. FDTD is fast, efficient, and intuitive for broadband problems as long as you respect its CFL stability limit. Implicit methods like Crank-Nicolson offer unconditional stability, a critical advantage when fine grids would otherwise force prohibitively small time steps. The core engineering skill is understanding these trade-offs to select, implement, and interpret the right method.

## Authors

- **Abdulrahman Albadawi** (031021047) — Researcher + Presenter
- **Abdisalam Hersi** (031022064) — Code + Presenter

## References

1. A Review of Progress in FDTD Maxwell's Equations Modeling of Impulsive Subionospheric Propagation Below 300 kHz
2. FDTD Maxwell's Equations Models for Nonlinear Electrodynamics and Optics
3. FDTD for Hydrodynamic Electron Fluid Maxwell Equations
4. High-Accuracy FDTD Solution of the Absorbing Wave Equation, and Conducting Maxwell's Equations Based on a Nonstandard Finite-Difference Model
5. Investigation of Numerical Errors of the Two-Dimensional ADI-FDTD Method for Maxwell's Equations Solution
6. Stable FEM-FDTD Hybrid Method for Maxwell's Equations
7. The Finite-Difference Time-Domain (FDTD) and the Finite-Volume Time-Domain (FVTD) Methods in Solving Maxwell's Equations
