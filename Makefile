make:
	mpicc poissonSolverWithoutRandom.c poissonProblemInputs.c -lm -lpthread -g -o poissonSolver.o
