#!/bin/sh

### test all usage outcomes and error messages
nProc=1
Nx=8
Ny=8
maxIter=1000
nThreads=1
ngp=1
casename=0

echo "// --- incorrect number of inputs --- //"
mpirun -n $nProc ./poissonSolver.o $Nx $Ny $maxIter $ngp $nThreads 
echo " "

echo "// --- bad grid size --- //"
mpirun -n $nProc ./poissonSolver.o -1 $Ny $maxIter $ngp $nThreads $casename
echo " "

echo "// --- bad grid size --- //"
mpirun -n $nProc ./poissonSolver.o $Nx -1 $maxIter $ngp $nThreads $casename
echo " "

echo "// --- bad nthreads --- //"
mpirun -n $nProc ./poissonSolver.o $Nx $Ny $maxIter $ngp -1 $casename
echo " "

echo "// --- bad ngp --- //"
mpirun -n $nProc ./poissonSolver.o $Nx $Ny $maxIter -1 $nThreads $casename
echo " "

echo "// --- bad casename --- //"
mpirun -n $nProc ./poissonSolver.o $Nx $Ny $maxIter $ngp $nThreads -1
echo " "

echo "// --- bad tol input --- //"
mpirun -n $nProc ./poissonSolver.o $Nx $Ny $maxIter $ngp $nThreads 0 -t
echo " "

echo "// --- bad tol input --- //"
mpirun -n $nProc ./poissonSolver.o $Nx $Ny $maxIter $ngp $nThreads 0 -v -t
echo " "
