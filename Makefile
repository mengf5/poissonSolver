make:
	mpicc poissonSolverWithoutRandom.c poissonProblemInputs.c -lm -lpthread -o poissonSolver.o
