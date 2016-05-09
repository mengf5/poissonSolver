make:
	mpicc poissonSolverWithoutRandom.c poissonProblemInputs.c -lm -I. -lpthread -g -o poissonSolver.o
