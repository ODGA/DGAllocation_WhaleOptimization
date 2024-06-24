# Optimal Allocation of Multiple Distributed Generation in IEEE 33-Bus Radial Distribution System
## Overview
This project focuses on optimizing the location and size/dispatch of multiple distributed generation (DG) units in the IEEE 33-Bus Radial Distribution System. The optimization is performed using Whale Optimization Algorithm (WOA) and Enhanced Whale Optimization Algorithm (E-WOA).

<image src="https://github.com/ODGA/DGAllocation_WhaleOptimization/assets/42317843/43400151-2e58-4fec-b8e7-fce55cbab384" width=40% height=40%> <image src="https://github.com/ODGA/DGAllocation_WhaleOptimization/assets/42317843/61d32ebb-aab3-4e7c-aa6d-0ea36c5e21ad" width=41% height=41%>


## Features
- Optimization of DG location and size/dispatch.
- Supports IEEE 33-Bus and 5-Bus Radial Distribution Systems.
- Customizable simulation parameters.
- Progress counter for simulation status.
- Intermediate results saving and continuation of interrupted simulations.
## Getting Started
### Prerequisites
- MATLAB R2024a or earlier.
### Hardware Used
- Intel i3 8100.
- 8GB RAM.
### Installation
1. Download the repository from GitHub.
2. Extract the folder to your computer.
3. Ensure all functions and saved data matrices are inside the same folder.
### Running Simulations

**Whale Optimization Algorithm (WOA)**

<image src="https://github.com/ODGA/DGAllocation_WhaleOptimization/assets/42317843/034c8776-40e3-4df8-9216-c2e522f92f2c" width=30% height=30%>


To run the simulation using WOA:

1. Open MATLAB.
2. Navigate to the project folder.
3. Run the _Run_WOA.m_ file.
   
**Enhanced Whale Optimization Algorithm (E-WOA)**

<image src="https://github.com/ODGA/DGAllocation_WhaleOptimization/assets/42317843/4dba5a99-0233-4e55-adc4-026b899b8ef5" width=30% height=30%>

To run the simulation using E-WOA:

1. Open MATLAB.
2. Navigate to the project folder.
3. Run the _Run_EWOA.m_ file.

### Customization

Before running the main files, you can modify the following parameters:

- Test system to use (5-bus or 33-bus).
- Number of distributed generation units to install (numDGs).
- Number of hours of load demand (numHours).
- Minimum and maximum capacity of each DG.
- Population size.
- Maximum iterations.
- Total runs.

Note: Using the 5-bus test system sets numDGs and numHours to 1 by default.

### Simulation Progress
A progress counter is displayed during the simulations. Upon reaching 100%, the simulations are complete.

### Results
After running the main file, the load flow for the best solutions is solved, and the optimal solutions and objective function values among the runs are obtained. These results, along with the simulation data, are stored in the matrix file:

_[ALGORITHM]_Results_[No.]Bus_Hour[No.].mat_

Where:

- _[ALGORITHM]_ could be WOA or EWOA.
- _[No.]_ Bus could be 5 or 33.
- Hour _[No.]_ could be 1, 3, 24, or 96.

### Additional Notes

For longer simulations, it is safe to stop the simulation. The main files allow for the continuation of the simulations. The results for each run are saved in the matrix file.
