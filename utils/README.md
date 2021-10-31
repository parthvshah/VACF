# Utilities

This folder contains the utilities used for the results obtained. Following is a brief description of the same:

- clean.py - Cleans the output file of LAMMPS and structures data for read by our programs
- diffusion.py - Calculates the coefficient of diffusion by intergrating the Velocity Auto-Correlation Function
- particles.py - Dissociates the LAMMPS output file into multiple files, each containing velocities of a particular particle
- plot.py - Plots the Velocity Auto-Correlation Function with `matplotlib`
- sub.py - A script to automatically attempt a submission of a scipt to PBS at a given interval