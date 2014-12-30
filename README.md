FCC-simulation
==============

This repository contains fundamental simulation means and methods in Molecular Dynamics, which are given as examplary programs of FCC crystal under hypothesis of Lennard-Jones Potential and Longevin thermal condition.

==============
Update log 20141218

now the radius distribution function can be calculated for the whole time steps, i.e., able to calculate all the 1000 time steps in the example.


==============
Update log 20141217

bug fixed for pair.cpp

==============
Update log 20141216

new feature add to this program: calculating one frame of atomcfg.xyz to get the radius distribution function via pair.cpp and pair.h.


==============
Update log 20141203

Till this version, this simulation program serves as the mean of simulating:

1. Simulate the FCC crystal viberation in real world.

2. The interactive force is determined according to Lennard-Jones Potential

3. The adjustment of the force is given by Longevin thermal condition

4. Use a customized parameter file, parameters.txt to customize the parameters needed for calculation.

The realization of FCC simulation is made up of several parts:

1. a.out: The executable program that is compiled by G++;

2. atomcfg.xyz: The output of the position of each atom at each time step, which can be used for further study, such as visualization in vmd.

3. log.dat: The thermal information for each time step, such as potential energy, knetic enery and temperature.

4. memory.h & memory.cpp: The class for creating high dimension array.

5. random.h & random.cpp: The class for creating random number with different distributions.

6. ThreeD.h & ThreeD.cpp: The class for initializing a FCC crystal, simluating its movements according to Lennard-Jones Potential and ouputing the result in standardlized .xyz and .dat files.

7. main.cpp: The main function that realizing the simulation.

8. parameters.txt: The file enablint customized data input.
