5/9/16
Instructions for BGQ:
 - change the code to have
   #define BGQ 1
   instead of 
   #define BGQ 0
 - change makefile from 
	mpicc poissonSolverWithoutRandom.c poissonProblemInputs.c -lm -I. -lpthread -g -o poissonSolver.o	
   to
	mpixlc poissonSolverWithoutRandom.c poissonProblemInputs.c -lm -I. -pthread -O3 -o poissonSolver.o	
 - run the command 
   chmod +x *.sh
   so you can run the shell files


Read me first:

There are several ways to get rid of the need for extern in the extra function files. I have done so in what I believe to be the simplest way, although it's not the cleanest. Here's how it works:

(1) "Public" variables: any global variables that will be used in both files must be declared in the exact same way in both files.
(2) External functions: functions that are designed in poissonProblemInputs must be declared in poissonSolverWithoutRandom before they are called. For example, the function initializeProblemInputs is designed in poissonProblemInputs, but it is delcared in poissonSolverWithoutRandom before main.

I have run several tests with this file and also with smaller test files, and everything works fine. If you are having trouble with this, please let me know.