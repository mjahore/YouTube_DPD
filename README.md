# Let's Write a DPD Simulation in C: A YouTube Miniseries

This code is the result of a series of YouTube videos aimed at briefly introducing the dissipative particle dynamics (DPD) simulation framework, and guiding the viewer through the simple process of writing a fairly efficient simulation from scratch. Although there are several great applications out there that can do DPD simulations (e.g., LAMMPS, HOOMD-blue, Materials Studio, and others), it can be very instructive to write one's own simulation. My mantra is that if I can't write it myself, I don't truly understand what is going on under the hood -- and I don't like to trust black boxes in science! For those who are interested, my research group maintains its own DPD code written in Fortran 77 that supports parallel/distributed computing using the Message Passing Interface (MPI). PD<sup>2</sup> is available at https://www.github.com/mjahore/PD2.

The YouTube miniseries can be viewed at: https://www.youtube.com/playlist?list=PLx9F_EFL1lVaAqnF-k091LjhjGvC0pLOq

This source is laid out to mirror the progress that is made in each video. The final, full code is contained in the Final_Code/ subdirectory of this repository. I've made every attempt to make the code as readable as possible with ample comments. I also did a little bit of optimization on the code in Advanced_Code/ to make it a little bit more general purpose. 

## Video 1:

In this video, I introduce the necessary fundamentals of DPD for writing a simulation, discuss the parameters needed for the simulation, and begin to write the code. I detail how we read the parameters from a file, create a random initial configuration, and write the particle coordinates to an .xyz file for viewing in an application like VMD or Ovito.

## Video 2:

In this video, I correct 3 minor mistakes from video 1, and discuss the use of a cell list to make the processing of our forces more efficient. In the code, I implement the cell list, harmonic force, and calculate the root mean-squared radius of gyration.

## Video 3:

In this video, I implement the cell list to efficiently calculate the conservative, random, and dissipative forces between all particle pairs separated by a distance of r<sub>c</sub> or less. I discuss the velocity-Verlet algorithm, and implement it along with two functions to calculate the average momentum and temperature in the system. I run a quick test and compare the output to values in the literature, showing that this code produces results that agree with published results with similar parameters but generated with PD<sup>2</sup>. PD<sup>2</sup> is a parallel implementation of DPD written in Fortran77 that utilizes MPI. It is open source and available at https://www.github.com/mjahore/PD2.

## Video 4:

In this video, I compare simulations of 5% polymer solutions (N = 10, 20, 30, 50, and 100) to polymer melts with identical chain lengths to observe how R<sub>g</sub> scales with (N - 1). We observe excellent quantitative agreement between the DPD results and classical polymer physics theories. I also demonstrate that the code can simulate phase separation of a binary mixture of immiscible fluids.
