#!/bin/bash

# biggest case in this file is
# Nproc = 32*64 = 2048
N=4096
# Jacobi converges with N^2 steps
# however, that will take forever. 
# instead just take a small portion of
# the steps. This takes into account
# sub iterations
nIterTot=$N

# submit batch jobs for each node config
# total batch jobs = 3*4 = 12
for nodes in 8 16 32
do
    # partition we are using
    partition=small
    
    # run for different ranks-thread configs

    # calculate total ranks
    nProc="$((64*$nodes))"
    nThreads=0
    
    # submit
    sbatch --partition $partition \
	--nodes $nodes \
	--time 60 \
	./bgq_single_run_p1.sh \
	$nProc \
	$N \
	$nIterTot \
	$nThreads \
	$partition

    nProc="$((16*$nodes))"
    nThreads=4
    
    sbatch --partition $partition \
	--nodes $nodes \
	--time 60 \
	./bgq_single_run_p1.sh \
	$nProc \
	$N \
	$nIterTot \
	$nThreads \
	$partition

    nProc="$((4*$nodes))"
    nThreads=16
    
    sbatch --partition $partition \
	--nodes $nodes \
	--time 60 \
	./bgq_single_run_p1.sh \
	$nProc \
	$N \
	$nIterTot \
	$nThreads \
	$partition

    nProc="$((1*$nodes))"
    nThreads=64
    
    sbatch --partition $partition \
	--nodes $nodes \
	--time 60 \
	./bgq_single_run_p1.sh \
	$nProc \
	$N \
	$nIterTot \
	$nThreads \
	$partition
	
done

