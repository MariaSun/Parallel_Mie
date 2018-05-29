# Parallel_Mie
The package performs Mie scattering routine in parallel using MPI interface.

To compile and run the package one would need to install OpenMP (say, https://github.com/firemodels/fds/wiki/Installing-Open-MPI-on-OS-X)

The paclage uses pre-excisting subroutines. The links and references are provided in the comments in first lines of the .c-files

Incoming plane-wave photon beam is generated into 1X1 mum ring. The scattered photons' position, weight and amplitude are calculated in photon.c (main-function). Library function mesh.c redistributes the photons into a mesh, adding up the wieghts of the photons arriving at the same sells of the detector.

MMA notebook gives the example of the output of the code. The corresponding pront statement is in mesh.c.
