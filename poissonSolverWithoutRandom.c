// Run with mpirun -n p ./a.out n maxIter thrds g

// --- Acitivity Log --- //
// 4/21/16 ER
// - Initial Set up of MPI and data parsing.
// - Setup the Maybe Function

// 4/22/16 DS
// - implemented Jacobi, advance grid
// - poissonProblemInputs.c

// 4/24/16 ER
// - Removed the need for extern in poissonProblemInputs.c. See READ_ME.txt

// 4/26/16 FM
// - confirmed y0 is not usable due to some wired internal function
// - couldnot figure out inline so I define forcingFunction before initilizationProblemInputs function 
// - change MPI_THREAD_MULTIPLE to MPI_THREAD_FUNNELED, do yours works with mutiple?
// - change ia,ib and add ja,jb
// - change forcingFunction and (b,t,l,r)bc 

// 4/26/16 DS
// - implemented applyBCs()
// - implemented residual error check
// - ran test with 1 rank, it converges!

// --- TODO --- //
// - initial guess needed
// - maxIter should account for sub-iterations
// - clean up the (1/dx)'s in the scheme
// - more testing but not constant solution
// - input tolerance as optional parameter?
// - communication
// - threads

#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<pthread.h>
#include<time.h>
#include<math.h>
#include<string.h>

// --- Type Defs --- //
// grid containers
typedef double** GRID;
typedef double* SUBGRID;

// --- Global Variables --- //
// MPI
MPI_Comm network; // This is MPI_COMM_WORLD
int nProc; // Number of processors
int myRank; // My rank
int nThreads; // Number of threads
int threadSupport; // level of desired thread support
int partition = 1;// partition scheme, 1 for square 2 for slender rect. // for now 1 is not available. 

// Grid
int Nx, Ny; // Global grid cells in x and y
int Npx, Npy; // Number of processors in x and y dimensions
int M; // Local grid cells in box
int Mx, My; // Local grid cells in rect. in x and y 
int N_matrix; // length of matrix

int ia,ib; // local index of i=0 to i=Mx
int ja,jb; // local index of j=0 to j=My
int I, J; // indices

// Note: y0 and y1 are taken due to some internal c reason
double x0_, x1_, y0_, y1_; // local axis limits
double X0_, X1_, Y0_, Y1_; // global axis limits
double dx, dy; // local grid spacing
double dX, dY; // global grid spacing

GRID Un; // solution at iteration
GRID Unp1; // solution at next iteration
GRID F; // forcing function
// boundary conditions
SUBGRID LBC;
SUBGRID RBC;
SUBGRID BBC;
SUBGRID TBC;

// Solver
int maxIter; // Maximum number of iterations before killing the program
int ngp; // Number of ghost cell layers to communicate
double tol = 1.e-06;
int casenumber;

// Flags
int verboseFlag;
char* verboseTag;
int printFrequency=10; // print iteration status every 10 iterations

// --- Function declarations --- //


// --- Usage --- //
// error handling
void usage(char* name, int error) {
  printf("Error %d in %s: \n", error, name);
  switch(error) {
  case 1:
    printf("Incorrect number of input arguments.\n");
    printf("Usage: %s Nx Ny maxIter ngp nThreads casename ",name);
    printf("[-t tol] [-v] \n");
    printf("       Nx, Ny   = grid size in x and y \n");
    printf("       maxIter  = maximum number of iterations \n");
    printf("       ngp      = number of ghost points \n");
    printf("       nThreads = number of threads \n");
    printf("       casename = number corresponding to case \n");
    printf("       -t       = optional flag for tolerance ");
    printf("specification \n");
    printf("       tol      = input argument for tolerance \n");
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
  case 6:
    printf("   Invalid casename selected.\n");
    break;
  case 7:
    printf("   [-t] was an input argument, but [tol] was not specified. \n");
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
  // max number of inputs = 10
  // min number of inputs = 7
  // [-t tol] and [-v] are optional 
  if (argc > 7) {
    if (strcmp(argv[7],"-v") == 0) {
      // first optional argument is [-v]
      if (myRank == 0) {
	verboseFlag = 1;
	verboseTag = "v>>";
	printf("%s Verbose output enabled \n",verboseTag);
      }

      // check for a second optional argument
      if ((argc > 8) && (strcmp(argv[8],"-t") == 0)) {
	// [-t] is an argument
	if (argc == 10) {
	  // set value of tol
	  tol = atof(argv[9]);
	} else if (argc > 10) {
	  return 1;
	} else {
	  // tol not specified
	  return 7;
	}
      }

    } else if (strcmp(argv[7],"-t") == 0) {
      // first optional argument is [-t]
      if (argc > 9) {
	// set value of tol
	tol = atof(argv[9]);
      } else {
	// tol not specified
	return 7;
      }

      // check for second optional argument
      if (argc == 10) {
	// [-v] is inputted
	if (myRank == 0) {
	  verboseFlag = 1;
	  verboseTag = "v>>";
	  printf("%s Verbose output enabled \n",verboseTag);
	}
      } else {
	return 1;
      }

    } else {
      return 1;
    }
  } else if (argc != 7) {
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

  // make sure number of processors
  // divides number of elements,
  // and the resulting square blocks
  // divide the number of rows and columns
  if (partition == 1){
    M = (int) sqrt(Nx*Ny/nProc);}
  else if(partition == 2){
    Mx = Nx;
    My = Ny/nProc;
  }

  // todo change this check
  // make sure this error check works for our situation
  /* if ((Nx*Ny-M*M*nProc != 0) || (Nx%M != 0) || (Ny%M != 0)) { */
  /*   return 2; */
  /* } */
  
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

  // processors in each dimension
  if (partition == 1){
    Npx = Nx/M;
    Npy = Ny/M;
  }
  else if (partition == 2) {
    Npx = 1;
    Npy = nProc;
  }

  // casenumber 
  casenumber = strtol(argv[6], &ptr, 10); 

  // todo print inputs

  // Let keep it this way assuming we nondimensionalize 
  // every problem to unit square

  // Global Axis limits
  // todo, define in input file
  X0_ = 0.;
  X1_ = 1.;
  Y0_ = 0.;
  Y1_ = 1.;

  // Grid spacing
  dx = (X1_-X0_)/(1.*Nx);
  dy = (Y1_-Y0_)/(1.*Ny);
  dX = (X1_-X0_)/(1.*Npx);
  dY = (Y1_-Y0_)/(1.*Npy);

  // Local Axis limits
  if (partition == 1){
    x0_ = X0_ + dX*(myRank%Npx);
    y0_ = Y0_ + dY*(myRank%Npx);}
  else if (partition == 1){
    x0_ = X0_;
    y0_ = Y0_ + dY*(myRank%nProc);
  }
  
  x1_ = x0_ + dX;
  y1_ = y0_ + dY;

  // print thread support for verbose option
  if (myRank == 0) {
    if (verboseFlag) {
      printf("%s ",verboseTag);
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
  //    int request = MPI_THREAD_MULTIPLE;
  // does yours support the multiple? 
    int request = MPI_THREAD_FUNNELED;
    MPI_Init_thread(&argc, &argv, request, &threadSupport);
    if (request != threadSupport) {
        if (myRank == 0) {
            printf("You requested level %d support but only have level %d support.\n", request, threadSupport);
        }
    }
  network = MPI_COMM_WORLD;
  MPI_Comm_size(network, &nProc);
  MPI_Comm_rank(network, &myRank);
}

// initialize Grid
void initializeGrid() {

  // calculate ia and ib
  if (partition == 1){
    ia = ngp;
    ib = M+ngp;
    ja = ngp;
    jb = M+ngp;
  }
  else if (partition == 2){
    ia = ngp;
    ib = Mx+ngp;
    ja = ngp;
    jb = My+ngp;
  }
  

  // todo, be wary of sizes of all arrays
  // 2D grid for solution
  Unp1 = (GRID) malloc(sizeof(SUBGRID)*(M+1+2*ngp));
  Un   = (GRID) malloc(sizeof(SUBGRID)*(M+1+2*ngp));
  // 2D grid for forcing function
  F    = (GRID) malloc(sizeof(SUBGRID)*(M+1+2*ngp));

  // 1D grid for boundaries
  LBC = (SUBGRID) malloc(sizeof(double)*(M+1+2*ngp));
  RBC = (SUBGRID) malloc(sizeof(double)*(M+1+2*ngp));
  BBC = (SUBGRID) malloc(sizeof(double)*(M+1+2*ngp));
  TBC = (SUBGRID) malloc(sizeof(double)*(M+1+2*ngp));

  // Allocate 
  for (I = ia-ngp; I <= ib+ngp; ++I) {
    Unp1[I] = (SUBGRID) malloc(sizeof(double)*(M+1+2*ngp));
    Un[I]   = (SUBGRID) malloc(sizeof(double)*(M+1+2*ngp));
    F[I]    = (SUBGRID) malloc(sizeof(double)*(M+1+2*ngp));

    LBC[I]  = 0.;
    RBC[I]  = 0.;
    BBC[I]  = 0.;
    TBC[I]  = 0.;

    // initialize all cells
    for (J = ja-ngp; J <= jb+ngp; ++J) {
      Unp1[I][J] = 0.;
      Un[I][J]   = 0.;
      F[I][J]    = 0.;
    }
  }
}

// Jacobi
void jacobiStep() {
  // Scheme is
  // 1/dx^2 * (U_{i+1,j}^n - 2 U_{i,j}^{n+1} + U_{i-1,j}^n) + ...
  // 1/dy^2 * (U_{i,j+1}^n - 2 U_{i,j}^{n+1} + U_{i,j-1}^n) = f_{i,j}

  // Todo, make it look nicer  
  double preMultiplier = (.5/(1./(dx*dx)+1./(dy*dy)));
  // Todo, change to ja and jb
  for (I = ia-ngp+1; I <= ib+ngp-1; ++I) {
    for (J = ja-ngp+1; J <= jb+ngp-1; ++J) {
      Unp1[I][J] = preMultiplier*
	( (1./(dx*dx))*(Un[I+1][J]+Un[I-1][J]) +
	  (1./(dy*dy))*(Un[I][J+1]+Un[I][J-1]) -
	  F[I][J] );
    }
  }
}

void advanceGrid() {
  // Copy grid Unp1 to Un
  // Todo, copy everything?
  for (I = ia-ngp; I <= ib+ngp; ++I) {
    for (J = ja-ngp; J <= jb+ngp; ++J) {
      Un[I][J] = Unp1[I][J];
    }
  }  
}

void updateBCs() {
  // Update boundary cells if applicable

  // left is a boundary
  if (x0_ == X0_) {
    I = ia;
    for (J = ja-ngp; J <= jb+ngp; ++J) {
      Unp1[I][J] = LBC[J];
    }
  }
  // bottom is a boundary
  if (y0_ == Y0_) {
    J = ja;
    for (I = ia-ngp; J <= ib+ngp; ++J) {
      Unp1[I][J] = BBC[J];
    }    
  }
  // right is a boundary
  if (x1_ == X1_) {
    I = ib;
    for (J = ja-ngp; J <= jb+ngp; ++J) {
      Unp1[I][J] = RBC[J];
    }
  }
  // top is a boundary
  if (y1_ == Y1_) {
    J = jb;
    for (I = ia-ngp; J <= ib+ngp; ++J) {
      Unp1[I][J] = TBC[J];
    }    
  }

}

double calculateMaximumResidualError() {

  double err = 0.;
  double check;

  for (I = ia+1; I <= ib-1; ++I) {
    for (J = ja+1; J <= jb-1; ++J) {
      // use residual of equation
      check = fabs( (1./(dx*dx))*(Unp1[I+1][J]-2.*Unp1[I][J]+Unp1[I-1][J]) +
		   (1./(dy*dy))*(Unp1[I][J+1]-2.*Unp1[I][J]+Unp1[I][J-1]) -
		   F[I][J] );
      if (check > err) {
	err = check;
      }
    }
  }

  return err;
}



// --- External Functions --- //
int initializeProblemInputs(int caseNumber );



// --- Main --- //
int main(int argc, char* argv[]) {

  // --- MPI Initialization --- //
  initializeMPI(argc, argv);
  
  // --- Parse input --- //
  int inputFlag = inputHandler(argc,argv);
  if (inputFlag != 0) {
    if (myRank == 0) {
      usage(argv[0],inputFlag);
    }
    MPI_Finalize();
    return 0;
  }

  if (verboseFlag) {
    printf("%s MPI initialization is done \n",verboseTag);
  }

  // --- Grid Initialization --- //
  initializeGrid();
  if (verboseFlag) {
    printf("%s Grid initialization is done \n",verboseTag);
  }

  // --- Problem Inputs --- //
  inputFlag = initializeProblemInputs(casenumber);
  if (inputFlag == 1) {
    if (myRank == 0) {
      usage(argv[0],6);
    }
    MPI_Finalize();
    return;
  }

  if (verboseFlag) {
    printf("%s \n",verboseTag);
    printf("%s Begin iterations: \n",verboseTag);
  }

  // --- Iteration Loop --- //
  int n, k; // TODO, better names for these vars
  double err;

  // perform iterations
  for (n = 0; n < maxIter; ++n){

    // perform sub iterations for the amount of
    // ghost points
    for (k = 0; k < ngp; ++k) {
      // --- Jacobi --- //
      jacobiStep();

      // --- Update BCs --- //
      updateBCs();
      
      // --- Advance --- //
      advanceGrid();
    }

    // todo figure this out
    MPI_Barrier(network);

    // --- Communication --- //
    // TODO, do this

    // todo implement this

    err = calculateMaximumResidualError();
    if (verboseFlag) {
      if (n%printFrequency == 0) {
	// todo, make a table, don't repeat n and err
	printf("    n = %4d: err = %6.3e \n",n,err);
      }
    }

    // todo, need to consider multiple ranks 
    if (err < tol) {
        break;
    }

    // communication(send and recieve parallel ghost lines)
    MPI_Barrier(network);
  }

  if (verboseFlag) {
    if (n == maxIter) {
      printf("%s Solver exited because the maximum ",verboseTag);
      printf("number of iterations were exceeded \n");
      printf("    maxIter = %6d \n",maxIter);
      printf("    err     = %6.3e \n",err);
      printf("    tol     = %6.3e \n",tol);
    } else {
      printf("%s Solver exited because the residual ",verboseTag);
      printf("error is less than the tolerance \n");
      printf("    nIter = %6d \n",n);
      printf("    err   = %6.3e \n",err); 
      printf("    tol   = %6.3e \n",tol); 
    }
  }
    
  MPI_Finalize();
  return 0;
}


