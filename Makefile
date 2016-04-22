make:
	mpicc poissonSolverWithoutRandom.c -lm -lpthread -o poissonSolver.o
