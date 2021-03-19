# GPU_Inhom_solver

This code is written by Kamalnath Kadivel and modified by Taiwu Yu. The inhomogeneous solver can be derived by 

Shen, Yao, et al. "An improvement on the three-dimensional phase-field microelasticity theory for elastically and structurally inhomogeneous solids." Scripta Materialia 60.10 (2009): 901-904.

The code can be downloaded through the following command:

git clone https://github.com/TaiwuYu/GPU_Inhom_solver.git

The main files include:

# constant.h

regulates system size and all the major variables

# time_evolution.c

programmes the governing equation of phase field model and calculates chemical potentials and elastic driving force

# inhom_v3.c

the inhomogeneous modulus solver aims to calculate the local stress field 



