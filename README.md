FCC-simulation
==============

This repository contains fundamental simulation program of FCC crystal under hypothesis of Lennard-Jones Potential and Longevin thermal condition in Molecular Dynamics.

Update log 20141203
Till this version, this simulation program serve as the mean of simulating:

1. Simulate the FCC crystal viberation in real world.

2. The interactive force is determined according to Lennard-Jones Potential

3. The adjustment of the force is given by Longevin thermal condition

4. Use a customized parameter file, parameters.txt to customize the parameters needed for calculation.

The realization of FCC simulation is made up of several parts:

1. a.out: The executable program that is compiled by G++;

2. atomcfg.xyz: The output of the position of each atom at each time step, which can be used for further study, such as visualization in vmd.

3. log.dat: The thermal information for each time step, such as potential energy, knetic nenery and temperature.

4. memory.h & memory.cpp

5. random.h & random.cpp

6. ThreeD.h & ThreeD.cpp

7. main.cpp

8. parameters.txt:
