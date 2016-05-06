make:
	mpicc poissonSolverWithoutRandom.c poissonProblemInputs.c -lm -I. -pthread -g -o poissonSolver.o
