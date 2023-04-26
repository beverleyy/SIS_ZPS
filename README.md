# Code for E2K session at Zhangde Primary

2 lessons on very-introductory Python coding and running a computational fluid dynamics simulation were conducted over 2 weeks at Zhangde Primary's E2K program as part of MOE's Scientists-in-Schools program, in collaboration with A*STAR IHPC. This repository contains the solver file, the Jupyter notebook provided to the students (`cavity_student_copy.ipynb`), and the instructor notebook (`cavity.ipynb`) which generates the animated lid-driven cavity flow solution.

## How to use

The codes in this repository were meant to be as IT-friendly as possible for use in a primary school environment. However, before running the codes, you must first install Jupyter on your computer. I recommend [Anaconda](https://www.anaconda.com/download) if user is not familiar with using a terminal, although you can also install standalone Python and use a package manager like pip to install Jupyter.

If you have Git, you can clone the repository: 
```
git clone https://github.com/beverleyy/SIS_ZPS
```
Otherwise you can download the entire repository as a ZIP archive and extract it to your drive.

If you have not installed `numpy` and `matplotlib` (both are included with standard Anaconda installations), install those two packages now.

Start Jupyter notebook in the repository folder, then open the relevant notebook and run it.

If you are running the Instructor version of the code, make a new directory to save the transient flow solutions in: `mkdir result` (MacOS/Linux/Windows PowerShell) or just make a new folder

## Implementation

As this exercise was created for a target audience of Primary 5 students, all the CFD code has been set aside in the file `solvers.py`. The code solves the incompressible Navier-Stokes equations in two dimensions using the finite volume method. Central differencing is applied for diffusion terms and upwind differencing for convective terms. Number of inner iterations for pressure solver can be controlled using `nit` variable.

The Jupyter notebook imports the `solve_cavity()` function from `solvers.py` and runs it per timestep. The student version will require input for the domain boundaries, grid size, timestep size, number of timesteps, and the lid velocity. 
