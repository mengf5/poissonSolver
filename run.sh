#!/bin/sh

# Usage: ./poissonSolver.o Nx Ny maxIter ngp nThreads casename [-t tol] [-v printFreq] 
#        Nx, Ny    = grid size in x and y 
#        maxIter   = maximum number of iterations 
#        ngp       = number of ghost points 
#        nThreads  = number of threads 
#        casename  = number corresponding to case 
#        -t        = optional flag for tolerance specification 
#        tol       = input argument for tolerance 
#        -v        = optional argument for verbose output 
#        printFreq = frequency to print iteration progress. Set negative for no iteration progress 
#                    Set to -2 for compact output

nProc=4
Nx=40
Ny=40
maxIter=100000
nThreads=1
ngp=11
casename=2
tol=1.e-10
printFreq=-2

# All these runs work
# mpirun -n $nProc ./poissonSolver.o $Nx $Ny $maxIter $ngp $nThreads $casename -t $tol -v
# mpirun -n $nProc ./poissonSolver.o $Nx $Ny $maxIter $ngp $nThreads $casename -t $tol -v $printFreq
# mpirun -n $nProc ./poissonSolver.o $Nx $Ny $maxIter $ngp $nThreads $casename -v $printFreq -t $tol 
#mpirun
#mpirun -n $nProc ./poissonSolver.o $Nx $Ny $maxIter $ngp $nThreads $casename -v $printFreq -t $tol 
mpirun -n $nProc ./poissonSolver.o $Nx $Ny $maxIter $ngp $nThreads $casename -v $printFreq -t $tol 

# for N in 80 160 320
# do
#     for iNgp in 1 2 4 8
#     do
# 	Ngp="$(($iNgp*$N/80))"

# 	Nt=1
# 	Np=8
# 	mpirun -n $Np ./poissonSolver.o $N $N $maxIter $Ngp $Nt $casename -v $printFreq -t $tol 

# 	Nt=2
# 	Np=4
# 	mpirun -n $Np ./poissonSolver.o $N $N $maxIter $Ngp $Nt $casename -v $printFreq -t $tol
	
# 	Nt=4
# 	Np=2
# 	mpirun -n $Np ./poissonSolver.o $N $N $maxIter $Ngp $Nt $casename -v $printFreq -t $tol 

# 	Nt=8
# 	Np=1
# 	mpirun -n $Np ./poissonSolver.o $N $N $maxIter $Ngp $Nt $casename -v $printFreq -t $tol 	
#     done    
# done
    
    



# for N in 20 40 80 160
# do
#     mpirun -n $nProc ./poissonSolver.o $N $N $maxIter $ngp $nThreads $casename -v -2 -t $tol 
# done
    
