#!/bin/sh

nProc=1
Nx=8
Ny=8
maxIter=1000
nThreads=1
ngp=1
casename=1
tol=1.e-07
printFreq=20

# All these runs work
# mpirun -n $nProc ./poissonSolver.o $Nx $Ny $maxIter $ngp $nThreads $casename -t $tol -v
# mpirun -n $nProc ./poissonSolver.o $Nx $Ny $maxIter $ngp $nThreads $casename -t $tol -v $printFreq
# mpirun -n $nProc ./poissonSolver.o $Nx $Ny $maxIter $ngp $nThreads $casename -v $printFreq -t $tol 
# mpirun -n $nProc ./poissonSolver.o $Nx $Ny $maxIter $ngp $nThreads $casename -v -t $tol 

mpirun -n $nProc ./poissonSolver.o $Nx $Ny $maxIter $ngp $nThreads $casename -v $printFreq -t $tol 
