# molecular-dynamics
# Simple 2D Molecular Dynamics Simulation

This repository contains C++ codes for a simple 2D Molecular Dynamics (MD) simulation of a two-component fluid interacting via the Lennard-Jones potential. Two different ensembles are provided: a Microcanonical (NVE) ensemble and a Canonical (NVT) ensemble using the Nosé-Hoover thermostat.

This code is intended as an educational example for learning the fundamentals of MD simulations.

---

## Features

-   **System**: 2D, 80-particle, two-component mixture.
-   **Potential**: Shifted Lennard-Jones (LJ) potential with a cutoff radius.
-   **Integrator**: Velocity Verlet algorithm.
-   **Thermostat**: Nosé-Hoover method for NVT simulations.
-   **Boundary Conditions**: Periodic Boundary Conditions (PBC) with the Minimum Image Convention (MIC).
-   **Initialization**: Generates a random initial configuration without particle overlaps and assigns initial velocities from a Maxwell-Boltzmann distribution (for NVT).

---

## Code Description

This repository includes two main source files:

-   `md_NVE.cpp`: Simulates the system in the **Microcanonical (NVE) ensemble**. In this simulation, the number of particles (N), volume (V), and total energy (E) are conserved.
-   `md_NVT.cpp`: Simulates the system in the **Canonical (NVT) ensemble**. The number of particles (N), volume (V), and temperature (T) are kept constant. Temperature control is achieved using the Nosé-Hoover thermostat.

---

