#!/bin/sh

nProc=1
Nx=8
Ny=8
maxIter=1000
nThreads=1
ngp=1
casename=0

mpirun -n $nProc ./poissonSolver.o $Nx $Ny $maxIter $ngp $nThreads $casename -v
