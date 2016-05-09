#!/bin/bash

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

# parse input args
nProc=$1
Nx=$2
Ny=$2
nIterTot=$3
nThreads=$4
partition=$5

# keep these constant
casename=2
tol=1.e-10
printFreq=-2 # this option prints our results neatly

# do 4 runs for each ghost point config
for iNgp in 1 2 4 8
do

    # choose number of ghost points
    Ngp="$(($iNgp*$N/80))"
    
    # total number of iterations
    maxIter="$(($nIterTot/$Ngp))"

    # run job
    srun --partition $partition \
	--ntasks $nProc --overcommit \
	./poissonSolver.o \
	$Nx $Ny $maxIter \
	$nThreads \
	$casename \
	-v $printFreq \
	-t $tol \
	>> results.txt
done    



