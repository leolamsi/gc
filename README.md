# DPGC Molecular Dynamics Simulation and Analysis Scripts

This repository contains the scripts necessary to perform molecular dynamics simulations and analyze the resulting trajectory data. The workflow involves generating simulation inputs, running the LAMMPS simulation, and then using MATLAB scripts for post-processing and analysis.


## 1.  **Generate Input Files:**

Use the python script `gen_init_ghost_nouniform.py` to generate the input file required for the LAMMPS simulation. This script creates the initialization file, defining particle pair interactions and particle positions.

## 2.  **Run LAMMPS Simulation:**

Execute the LAMMPS simulation using the script `run_sim.lmp`. This script is the primary LAMMPS input file that drives the simulation.

## 3.  **Analyze Simulation Data:**    
The following MATLAB scripts are used to analyze the simulation trajectory data:
     *`trajread_leo.m`: This script reads the LAMMPS trajectory dump file, unwraps periodic boundary conditions, and calculates the mean squared displacement (MSD), relaxation time (τ), and diffusion coefficient (D). This script calls all the other MATLAB scripts. This is the only script that has to be run on MATLAB manually.
     *`msdcal_leo.m`: This script calculates the MSD as a function of time.
     *`D_from_msd.m`: This script calculates the diffusion coefficient (D) and relaxation time (τ) from the MSD data. 
     *`fscal.m`: This script calculates the self-intermediate scattering function (SISF) and overlap. 



###  Analysis details: `trajread_leo.m`

This MATLAB function reads the LAMMPS trajectory dump file, unwraps periodic boundary conditions, and calls all the other analysis functions.

* **Input:**
    * `filename`:  Trajectory filename.
    * `T`: Temperature.
* **Output files:**
    * A file ending with `.msd` containing the time and MSD.
    * A file ending with `.tauD` containing the temperature (T), relaxation time (τ), and diffusion coefficient (D). (This data is not used for our manuscript as it is simply and estimation)
    * A file ending with `.fs` containing the time, SISF, and overlap.

