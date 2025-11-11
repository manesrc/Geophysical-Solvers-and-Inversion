# Physically-Consistent Solvers for Geophysical Inversion

This repository contains the core numerical methods developed during my Ph.D., "Numerical techniques for the solution of thermal problems in the context of geophysical inversions."

This repository is submitted as a request from my Ph.D. jury asking for a Manual on the developed methods.

---

## Project & Code Overview

The full research project involved developing a complete framework to solve the geophysical inverse problem of identifying the Lithosphere-Asthenosphere Boundary (LAB) from Surface 
Heat Flux (SHF) data, while ensuring physical consistency.

This framework consists of three main components:

1.  **[Solver 1 - UPLOADED] Subdivision solver (using Nitsche's method):** A physically-consistent solver that splits the domain (Lithosphere/Asthenosphere) and solves two separate 
problems, enforcing the continuity of temperature and flux at the interface -> supporting journal paper DOI: 10.1108/hff-10-2023-0649.
2.  **[Solver 2 - UPLOADED] Convection Solver:** A novel, single-domain solver where convection velocities are determined to naturally adjust the isotherm value at 
the interface -> supporting journal paper in preparation.
3.  **[Inverse Framework - To be uploaded] MCMC Inversion:** The "englobing" code that uses the two solvers above as an engine to solve the inverse problem, using Bayesian MCMC to 
quantify the uncertainty of the LAB's geometry  -> supporting journal paper in preparation

### Repository Status

This repository currently contains the **full MATLAB code for the first and second solvers (Solver 1: Nitsche's Method and Solver 2: Convection Solver)**.
* The core functions and modules are located in the uploaded folders each describing specific tasks developed during the code execution.
* The main script to run the first solver is 'ComputeTemp_NitscheMethod.m' 
* The main script to run the second solver is 'main_uAdjustment.m'
        * within the second solver there exist an iterative procedure to obtain 'equilibrium LAB configurations': those LAB geometries that produce a Stokes problem solution that when placed
                                                                                                             in a thermal problem their result comply the isotherm condition at the geometry of the LAB
        * also there exist different codes to produce velocity basis (two kinds available, with and without distance-related attenuation of viscosity)

The code for tthe MCMC inversion framework is currently being cleaned and prepared for publication and will be uploaded in the near future.

### How to Run (Solver 1)

* **Language:** MATLAB R2024a 
* **Core Requirements:** `OptimizationToolbox` 

---

## Skills and knowledge needed for these tools 

* **Developing novel methods:** The solvers were built from conservation equations to solve the problem of physical inconsistency in existing models.
* **Solving inverse problems:** The code is the "forward model" engine designed for use within a full geophysical inversion framework.
* **Quantifying uncertainties:** The full framework uses this solver within a Bayesian MCMC loop.
* **Validate rigorousity:** The methods were validated against analytically-derived benchmark problems and real-world data (SHF).
* **Programming in MATLAB:** Advanced proficiency in MATLAB for numerical simulation.
* **Writing and documenting, reproducible scientific software:** The code shows the ability to work among colleagues to produce software.
  
---

## License

This project is licensed under the MIT License. Please see the `LICENSE` file for details.
