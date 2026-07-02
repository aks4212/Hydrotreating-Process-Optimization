# Hydrotreating Process Optimization

This repository contains the rigorous thermodynamic modeling, validation, and multi-variable optimization of a heavy diesel Hydrotreating Process featuring an H2 recycle loop. The project leverages both DWSIM and MATLAB to evaluate sulfur conversion efficiency against CAPEX/OPEX trade-offs.

---

## Project Structure

* **`DWSIM_Simulation/`**: Contains the core flowsheet models simulating 8 distinct unit operations, closed mass/energy balances, and reaction kinetics across 11 chemical species.
* **`MATLAB_Validation/`**: Computational scripts utilizing the Chao-Seader thermodynamic property package to independently validate the DWSIM flowsheet results.
* **`Optimization_Data/`**: Datasets and plotting scripts for the Pareto optimization, mapping process variables against economic and efficiency constraints.

---

## Core Concepts & Technologies

* **Process Simulation (DWSIM & MATLAB):** Engineered a highly constrained simulation environment, ensuring strict mathematical closure with 0 Degrees of Freedom (DOF=0).
* **Thermodynamic Modeling:** Applied Chao-Seader thermodynamics to accurately predict vapor-liquid equilibrium (VLE) in complex hydrocarbon and hydrogen mixtures.
* **Reaction Kinetics & Recycle Loops:** Modeled complex desulfurization kinetics within a continuous H2 recycle system to maximize resource efficiency.
* **Multi-Objective Optimization:** Conducted Pareto optimization across 4 distinct process variables to balance theoretical chemical yield with real-world energy economics.

---

## Execution & Setup

1. Clone the repository:
   ```bash
   git clone [https://github.com/ask4212/Diesel-Hydrotreating-Optimization.git](https://github.com/aks4212/Diesel-Hydrotreating-Optimization.git)
