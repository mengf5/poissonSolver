// Run with mpirun -n p ./a.out n maxIter thrds g

// --- Acitivity Log --- //
// 4/21/16 ER
// - Initial Set up of MPI and data parsing.
// - Setup the Maybe Function

#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<pthread.h>
#include<time.h>

// --- Global Variables --- //
// MPI
MPI_Comm network; // This is MPI_COMM_WORLD
int nProc; // Number of processors
int myRank; // My rank
int nThreads; // Number of threads
int threadSupport; // level of desired thread support

// Grid
int Nx, Ny; // Grid cells in x and y
int N_matrix; // length of matrix

// Solver
int maxIter; // Maximum number of iterations before killing the program
int ngp; // Number of ghost cell layers to communicate

// Flags
int verboseFlag;

// --- Usage --- //
// error handling
void usage(char* name, int error) {
  printf("Error %d in %s: \n", error, name);
  switch(error) {
  case 1:
    printf("Incorrect number of input arguments.\n");
    printf("Usage: %s Nx Ny maxIter ngp nThreads [-v] \n",name);
    printf("       Nx, Ny   = grid size in x and y \n");
    printf("       maxIter  = maximum number of iterations \n");
    printf("       ngp      = number of ghost points \n");
    printf("       nThreads = number of threads \n");
    printf("       -v       = optional argument for verbose output \n");
    break;
  case 2:
    printf("   Unable to parse size of the grid.\n");
    printf("   The size of the grid should be a positive integer.\n");
    break;
  case 3:
    printf("   Unable to parse the maximum number of iterations.\n");
    printf("   Maximum iterations should be a positive integer.\n");
    break;
  case 4:
    printf("   Unable to parse the number of threads.\n");
    printf("   The number of threads should be a positive integer.\n");
    break;
  case 5:
    printf("   Unable to parse the number of ghost layers.\n");
    printf("   The number of ghost layers should be a positive integer.\n");
    break;
  default:
    printf("   It is literally impossible to get here in the program.\n");
    printf("   This must be the work of the divine Maybe function!\n");
    break;
  }
}

// input handling
int inputHandler(int argc, char* argv[]) {
  // verbose option
  verboseFlag = 0;

  // Incorrect number of input arguments
  printf("%d",argc);
  if (argc == 7) {
    if (strcmp(argv[6],"-v") == 0) {
      verboseFlag = 1;
      if (myRank == 0) {
	printf("%s: Verbose output enabled \n",argv[0]);
      }
    } else {
      return 1;
    }
  } else if (argc != 6) {
    return 1;
  }

  char* ptr;
  // Parse the size of the grid
  Nx = strtol(argv[1], &ptr, 10); 
  if (*ptr != '\0' || Nx < 1) {
    return 2;
  }

  Ny = strtol(argv[2], &ptr, 10); 
  if (*ptr != '\0' || Ny < 1) {
    return 2;
  }
  
  // Parse the number of iterations before death
  maxIter = strtol(argv[3], &ptr, 10); 
  if (*ptr != '\0' || maxIter < 1) {
    return 3;
  }
  // Parse the number of ghost cell layers to communicate
  ngp = strtol(argv[4], &ptr, 10); 
  if (*ptr != '\0' || ngp < 1) {
    return 5;
  }
  
  // Parse the nubmer of threads
  nThreads = strtol(argv[5], &ptr, 10); 
  if (*ptr != '\0' || nThreads < 1) {
    return 4;
  }

  // print thread support for verbose option
  if (myRank == 0) {
    if (verboseFlag) {
      switch(threadSupport) {
      case 0:
	printf("Only one thread will execute.\n");
	break;
      case 1:
	printf("Only the main thread will make MPI calls.\n");
	break;
      case 2:
	printf("Only one thread will make MPI calls at a time.\n");
	break;
      case 3:
	printf("All threads may make MPI calls without restriction.\n");
	break;
      default:
	printf("Program error.\n");
	return 6;
      }
    }
  }
  
  return 0;
}

// initialize MPI
void initializeMPI(int argc, char* argv[]) {
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &threadSupport);
  network = MPI_COMM_WORLD;
  MPI_Comm_size(network, &nProc);
  MPI_Comm_rank(network, &myRank);
}

// --- Main --- //
int main(int argc, char* argv[]) {

  // --- MPI Initialization --- //
  initializeMPI(argc,argv);
  
  // --- Parse input --- //
  int inputFlag = inputHandler(argc,argv);
  if (inputFlag != 0) {
    if (myRank == 0) {
      usage(argv[0],inputFlag);
    }
    MPI_Finalize();
    return 0;
  }
  
  // --- Matrix Initialization --- //
  /* int nLayer = 0; */
  /* int i,j; */
  /* for(int i = 0; i < maxIter; ++i){ */
  /*   for(int j = 0; j < nLayer; ++j){  */
  /*     //caculation in each mpirank */
  /*     //(jacobi/CG/..) */
  /*   } */

  /*   MPI_Barrier(network); */

  /*   // if(err < criterion) */
  /*   //     break;   */

  /*   //communication(send and recieve parallel ghost lines) */


  /*   MPI_Barrier(network); */


  /* } */
    
    
    
    
    
    
    
    
    
    
  MPI_Finalize();
  return 0;
}


// --- Functions --- //

