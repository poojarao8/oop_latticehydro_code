# oop_latticehydro_code
This repo contains the code for solving the lattice fluids equations for incompressible fluids using C++ and parallelized using MPI.The pressure Poisson equation is solved using the HYPRE library. Currently, the solve used for pressure equation is pre-conditioned conjugate gradient, but it will be changed to a multigrid solver (Boomer AMG). 
The next significant task is to incorporate hybrid parallelization using GPUs.


