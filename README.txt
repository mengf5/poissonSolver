%-----------------------------------------------------------------------------------
This is a description for how to make and run our parallel Poisson equation solver in BG/Q
%-----------------------------------------------------------------------------------

a) Make executable::

$ module load xl
$ make 

b) Run the command

The usage information is read in to the executable from::
<bgq_single_run_gp.sh>
A detailed introduction of the input variables is in the beginning of this file.

We specify these input variables in(If you want to change any input, dive in this)::
<bgq_run_cases_new.sh>

To run the case, you need to,

$ chmod +x *.sh
$ bgq_run_cases_new.sh

%-----------------------------------------------------------------------------------
Changes for Kratos or your local machine
%-----------------------------------------------------------------------------------
 - change the code <poissonSolverWithoutRandom.c> to have
   #define BGQ 0
   instead of 
   #define BGQ 1
 - change makefile from 
	mpixlc poissonSolverWithoutRandom.c poissonProblemInputs.c -lm -I. -pthread -O3 -o poissonSolver.o
   to::
	mpicc poissonSolverWithoutRandom.c poissonProblemInputs.c -lm -I. -lpthread -g -o poissonSolver.o	
 


%-----------------------------------------------------------------------------------
Group members and contributions
%-----------------------------------------------------------------------------------
Total work of this project is separated in the following sections. We put the numbers we worked on after our name. 

  Code::
1. implemented the algorithm,
2. implemented the parallel part of the code,
3. implemented the supported part in the code(usage, error checking, etc.)

  Testing::
4. Checked the convergence rate of the code and make sure the code works on our local machine.

Runs and Results::
5. Ruined performance tests on BG\Q,
6. Parsed output data, made presentable tables and plots.

Paper::
7. Problem specification
8. Code Implementation
9. Results
10. Summary and future work
11. References


Fanlong Meng, 1,2,3,4,5,6,9
Edward Rusu, 1,3,7,8,9,10,11
Daniel Serino, 1,2,3,4,5*,6,10,11


