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
// - implemented applyBCs() (and debugged)
// - implemented residual/global error check
// - ran tests with 1 rank, it converges!
// - see testResults.txt
// - rotated the slenderness b/c U[I][J], we want to send messages with blocks in y

// 5/5/16 all
// - debugged comm and ghost points

// --- TODO --- //
// - initial guess needed
// - maxIter should account for sub-iterations
// - clean up the (1/dx)'s in the scheme
// - more testing but not constant solution
// - input tolerance as optional parameter?
// - communication
// - threads
// - may need slender in other direction

// 4/28/16 FM
// - add array TGL/BGL for communication ghost lines
// - implemented MPI_Isend and Irecv for send/recv ghost lines
// - implemented updateGhosts to send Ghost line to Un
// - didn't mess up the serial code.
// 4/29/16 FM
// - change the comm to slender block
// - form a vector for all ghost lines and send them together.

// --- TODO --- //
// - run some test on this
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
MPI_Request requestSend1, requestRecv1;
MPI_Request requestSend2, requestRecv2;
MPI_Status status1;
MPI_Status status2;

int nProc; // Number of processors
int myRank; // My rank
int nThreads; // Number of threads
int threadSupport; // level of desired thread support
int partition = 2; // partition scheme, 1 for square 2 for slender rect.
// for now 1 is not available. 

// Grid
int Nx, Ny; // Global grid cells in x and y
int Npx, Npy; // Number of processors in x and y dimensions
int Mx, My; // Local grid cells in rect. in x and y 
int N_matrix; // length of matrix
int Mx_t; // number of grid cells per thread

int ia,ib; // local index of i=0 to i=Mx
int ja,jb; // local index of j=0 to j=My
int I, J, J1, G; // indices


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
// communication conditions
SUBGRID LGL;
SUBGRID RGL;
SUBGRID LGLp1;
SUBGRID RGLp1;


// Solver
int maxIter; // Maximum number of iterations before killing the program
int ngp; // Number of ghost cell layers to communicate
double tol = 1.e-06;
int casenumber;
int checkSolution;

// Flags
int verboseFlag;
char* verboseTag;
int printFrequency=-1; // print iteration status every 10 iterations

// --- Function declarations --- //


// --- Usage --- //
// error handling
void usage(char* name, int error) {
  printf("Error %d in %s: \n", error, name);
  switch(error) {
  case 1:
    printf("Incorrect number of input arguments.\n");
    printf("Usage: %s Nx Ny maxIter ngp nThreads casename ",name);
    printf("[-t tol] [-v printFreq] \n");
    printf("       Nx, Ny    = grid size in x and y \n");
    printf("       maxIter   = maximum number of iterations \n");
    printf("       ngp       = number of ghost points \n");
    printf("       nThreads  = number of threads \n");
    printf("       casename  = number corresponding to case \n");
    printf("       -t        = optional flag for tolerance ");
    printf("specification \n");
    printf("       tol       = input argument for tolerance \n");
    printf("       -v        = optional argument for verbose output \n");
    printf("       printFreq = frequency to print iteration progress. ");
    printf("Set negative for no iteration \n");
    printf("                   progress. Set to -2 for compact output \n");
    break;
  case 2:
    printf("Unable to parse size of the grid.\n");
    printf("The size of the grid should be a positive integer.\n");
    break;
  case 3:
    printf("Unable to parse the maximum number of iterations.\n");
    printf("Maximum iterations should be a positive integer.\n");
    break;
  case 4:
    printf("Unable to parse the number of threads.\n");
    printf("The number of threads should be a positive integer.\n");
    break;
  case 5:
    printf("Unable to parse the number of ghost layers.\n");
    printf("The number of ghost layers should be a positive integer.\n");
    break;
  case 6:
    printf("Invalid casename selected.\n");
    break;
  case 7:
    printf("[-t] was an input argument, but [tol] was not specified. \n");
    break;
  case 8:
    printf("partition = 1, the grid must be able to be divided into square ");
    printf("patches for each rank. \n");
    printf("The squares must be fitted in x and y. \n");
    break;
  case 9:
    printf("partition = 2, the grid must be able to be divided into slender ");
    printf("rectangles. \n");
    printf("The squares must be fitted in y. \n");
    break;
  case 10:
    printf("Invalid parition choice, must choose 1 for squares or 2 for slender rectangles\n");
    break;
  default:
    printf("It is literally impossible to get here in the program.\n");
    printf("This must be the work of the divine Maybe function!\n");
    break;
  }
}

// input handling
int inputHandler(int argc, char* argv[]) {
  // verbose option
  verboseFlag = 0;

  // Incorrect number of input arguments

  // max number of inputs = 11
  int maxIn = 11;
  // min number of inputs = 7
  int minIn = 7;
  // [-t tol] and [-v printFrequency] are optional 
  if (argc > minIn) {
    if (strcmp(argv[7],"-v") == 0) {
      // first optional argument is [-v]
      // check for a second optional argument
      if ((argc > 8) && (strcmp(argv[8],"-t") == 0)) {
	// [-t] is the next argument and 
	// printFrequency is not specified, therefore
	// max # of inputs is now (maxIn-1)
	if (argc == maxIn-1) { 
	  // set value of tol
	  tol = atof(argv[maxIn-1-1]);
	} else if (argc > maxIn-1) {
	  // wrong number of input args
	  return 1;
	} else {
	  // tol not specified
	  return 7;
	}
      } else if (argc > 8) {
	// printFrequency is the next argument
	printFrequency = atoi(argv[8]);
	
	// [-t] is the next argument and 
	// printFrequency is specified, therefore
	// max # of inputs is (maxIn)
	if (argc == maxIn) { 
	  // set value of tol
	  tol = atof(argv[maxIn-1]);
	} else if (argc > maxIn) {
	  // wrong number of input args
	  return 1;
	} else {
	  // tol not specified
	  return 7;
	}
      }

      // enable verbose output
      if ((myRank == 0) && (printFrequency != -2)) {
	verboseFlag = 1;
	verboseTag = "v>>";
	printf("%s Verbose output enabled \n",verboseTag);
      }

    } else if (strcmp(argv[7],"-t") == 0) {
      // first optional argument is [-t]
      if ((argc > 8) && (strcmp(argv[8],"-v") != 0)) {
	// set value of tol
	tol = atof(argv[8]);
      } else {
	// tol not specified
	return 7;
      }

      // check for second optional argument
      if (argc > 9) {
	// [-v] is inputted
	if (argc == maxIn) {
	  printFrequency = atoi(argv[maxIn-1]);
	} else if (argc > maxIn) {
	  return 1;
	}
	// check for printFrequency input

	// enable verbose output
	if ((myRank == 0) && (printFrequency != -2)) {
	  verboseFlag = 1;
	  verboseTag = "v>>";
	  printf("%s Verbose output enabled \n",verboseTag);
	}
      } else if (argc != 9) {
	// incorrect number of inputs
	return 1;
      }

    } else {
      // incorrect number of inputs
      return 1;
    }
  } else if (argc != minIn) {
    // incorrect number of inputs
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
    Mx = (int) sqrt(Nx*Ny/nProc);
    My = Mx;
  } else if (partition == 2){
    Mx = Nx/nProc;
    My = Ny;
  }

  // todo change this check
  // make sure this error check works for our situation
  if (partition == 1) {
    if ((Nx*Ny-Mx*My*nProc != 0) || (Nx%Mx != 0) || (Ny%My != 0)) {
      return 8;
    }
  } else if (partition == 2) {
    if (Nx%Mx != 0) {
      return 9;
    }
  } else {
    return 10;
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
  if (nThreads == 0) {
    Mx_t = Mx;
  } else {
    Mx_t = Mx/nThreads;
  }
  if ((*ptr != '\0') || (nThreads < 0)) {
    return 4;
  }
  if ((nThreads != 0) && (nThreads*Mx_t != Mx)) {
    return 4;
  }

  // processors in each dimension
  if (partition == 1){
    Npx = Nx/Mx;
    Npy = Ny/My;
  } else if (partition == 2) {
    Npx = nProc;
    Npy = 1;
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
  else if (partition == 2) {
    x0_ = X0_ + dX*(myRank%nProc);
    y0_ = Y0_;
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

  ia = ngp;
  ib = Mx+ngp;
  ja = ngp;
  jb = My+ngp;

  // todo, be wary of sizes of all arrays
  // 2D grid for solution
  Unp1 = (GRID) malloc(sizeof(SUBGRID)*(Mx+1+2*ngp));
  Un   = (GRID) malloc(sizeof(SUBGRID)*(Mx+1+2*ngp));
  // 2D grid for forcing function
  F    = (GRID) malloc(sizeof(SUBGRID)*(Mx+1+2*ngp));

  // 1D grid for boundaries
  LBC = (SUBGRID) malloc(sizeof(double)*(My+1+2*ngp));
  RBC = (SUBGRID) malloc(sizeof(double)*(My+1+2*ngp));
  BBC = (SUBGRID) malloc(sizeof(double)*(Mx+1+2*ngp));
  TBC = (SUBGRID) malloc(sizeof(double)*(Mx+1+2*ngp));
  // 1D grid for communication ghost line
  // p1 is the recieving bay, the others will be sent
  LGL = (SUBGRID) malloc(sizeof(double)*(My+1+2*ngp)*ngp);
  RGL = (SUBGRID) malloc(sizeof(double)*(My+1+2*ngp)*ngp);
  LGLp1 = (SUBGRID) malloc(sizeof(double)*(My+1+2*ngp)*ngp);
  RGLp1 = (SUBGRID) malloc(sizeof(double)*(My+1+2*ngp)*ngp);

  // Allocate 
  for (I = ia-ngp; I <= ib+ngp; ++I) {
    Unp1[I] = (SUBGRID) malloc(sizeof(double)*(My+1+2*ngp));
    Un[I]   = (SUBGRID) malloc(sizeof(double)*(My+1+2*ngp));
    F[I]    = (SUBGRID) malloc(sizeof(double)*(My+1+2*ngp));

    BBC[I]  = 0.;
    TBC[I]  = 0.;

    // initialize all cells
    for (J = ja-ngp; J <= jb+ngp; ++J) {
      Unp1[I][J] = 0.;
      Un[I][J]   = 0.;
      F[I][J]    = 0.;      
    }
  }

  for (J = ja-ngp; J <= jb+ngp; ++J) {
    LBC[J]  = 0.;
    RBC[J]  = 0.;
  }

  for (G = 1; G<= ngp; ++G){
    for (J = ja-ngp; J <= jb+ngp; ++J) {
      J1 = J + (My+1+2*ngp)*(G-1);
      LGL[J1] = Un[ia-G][J];
      RGL[J1] = Un[ib+G][J];
      LGLp1[J1] = Un[ia-G][J];
      RGLp1[J1] = Un[ib+G][J];
    }
  }
}

// Jacobi
void jacobiStep(int iStart,int iEnd) {
  // Scheme is
  // 1/dx^2 * (U_{i+1,j}^n - 2 U_{i,j}^{n+1} + U_{i-1,j}^n) + ...
  // 1/dy^2 * (U_{i,j+1}^n - 2 U_{i,j}^{n+1} + U_{i,j-1}^n) = f_{i,j}

  // dy^2 * (U_{i+1,j}^n + U_{i-1,j}^n) + ...
  // dx^2 * (U_{i,j+1}^n + U_{i,j-1}^n) - ...
  // dx^2 dy^2 * (f_{i,j}) = ...
  // 2 (dy^2 + dx^2) * (U_{i,j})

  /* printf("rank %d: in jacobi \n",myRank); */
  // Todo, make it look nicer  
  double A_ = .5*(dy*dy)/(dy*dy+dx*dx);
  double B_ = .5*(dx*dx)/(dy*dy+dx*dx);
  double C_ = -.5*(dx*dx*dy*dy)/(dy*dy+dx*dx);
  // Todo, change to ja and jb
  /* for (I = ia-ngp+1; I <= ib+ngp-1; ++I) { */
  for (I = iStart; I <= iEnd; ++I) {
    for (J = ja-ngp+1; J <= jb+ngp-1; ++J) {
      /* Unp1[I][J] = preMultiplier* */
      /* 	( (1./(dx*dx))*(Un[I+1][J]+Un[I-1][J]) + */
      /* 	  (1./(dy*dy))*(Un[I][J+1]+Un[I][J-1]) - */
      /* 	  F[I][J] ); */
      /* printf("rank %d: in loop \n",myRank); */
      /* printf("i+1 %f \n",Un[I+1][J]); */
      /* printf("i-1 %f\n",Un[I-1][J]); */
      /* printf("j+1 %f\n",Un[I][J+1]); */
      /* printf("j-1 %f\n",Un[I][J-1]); */
      Unp1[I][J] = A_*(Un[I+1][J]+Un[I-1][J]) +
	B_*(Un[I][J+1]+Un[I][J-1]) +
	C_*F[I][J];
      /* printf("rank %d: (%d,%d)\n",myRank,I,J); */
    }
  }
}

void advanceGrid(int iStart, int iEnd) {
  // Copy grid Unp1 to Un
  // Todo, copy everything?
  /* for (I = ia-ngp; I <= ib+ngp; ++I) { */
  for (I = iStart; I <= iEnd; ++I) {
    for (J = ja-ngp; J <= jb+ngp; ++J) {
      Un[I][J] = Unp1[I][J];
    }
  }
}

void updateBCs(int iStart, int iEnd) {
  // Update boundary cells if applicable

  // left is a boundary, and thread handles BC
  /* printf("istart = %d, ia-ngp = %d \n",iStart,ia-ngp); */
  if ((x0_ == X0_) && (iStart == ia-ngp)) {
    I = ia;
    /* printf("fuckfuckfuckfuckfuckfuck \n"); */
    for (J = ja-ngp; J <= jb+ngp; ++J) {
      Unp1[I][J] = LBC[J];
      /* printf("u(%d,%d) = %9.3e, or %9.3e \n",I-ia,J-ja,Unp1[I][J],LBC[J]); */
    }
  }
  // bottom is a boundary
  if (y0_ == Y0_) {
    J = ja;
    /* for (I = ia-ngp; I <= ib+ngp; ++I) { */
    for (I = iStart; I <= iEnd; ++I) {
      Unp1[I][J] = BBC[I];
    }    
  }
  // right is a boundary, and thread handles BC
  /* printf("iend = %d, ib+ngp = %d \n",iEnd,ib+ngp); */
  if ((x1_ == X1_) && (iEnd == ib+ngp)) {
    I = ib;
    for (J = ja-ngp; J <= jb+ngp; ++J) {
      Unp1[I][J] = RBC[J];
    }
  }
  // top is a boundary
  if (y1_ == Y1_) {
    J = jb;
    /* for (I = ia-ngp; I <= ib+ngp; ++I) { */
    for (I = iStart; I <= iEnd; ++I) {
      Unp1[I][J] = TBC[I];
    }    
  }


}

void collectGhosts() {
  // Update boundary cells if applicable


  // form left ghost array
  if (x0_ != X0_) {
    for (G = 1; G<= ngp; ++G){
      for (J = ja-ngp; J <= jb+ngp; ++J) {
	J1 = J + (My+1+2*ngp)*(G-1);
	LGL[J1] = Un[ia+G][J];
	/* printf("LGL[%d] = %f \n",J1,LGL[J1]); */
      }
    }
  }

  // form right ghost array
  if (x1_ != X1_) {
    for (G = 1; G<= ngp; ++G){
      for (J = ja-ngp; J <= jb+ngp; ++J) {
	J1 = J + (My+1+2*ngp)*(G-1);
	RGL[J1] = Un[ib-G][J];
	/* printf("RGL[%d] = %f \n",J1,RGL[J1]); */
      }
    }
  }

}

void updateGhosts() {
  // Update boundary cells if applicable
  /* printf("rank %d: \n",myRank); */
  // form left ghost array
  if (x0_ != X0_) {
    for (G = 1; G<= ngp; ++G){
      for (J = ja-ngp; J <= jb+ngp; ++J) {
	J1 = J + (My+1+2*ngp)*(G-1);
	/* printf("rank %d: LGLp1[%d] = %f \n",myRank,J1,LGLp1[J1]); */
	Un[ia-G][J]  = LGLp1[J1];
	/* printf("rank %d: (I,J) = (%d,%d)",myRank,ia-G,J); */

      }
    }
  }

  // form right ghost array
  if (x1_ != X1_) {
    for (G = 1; G<= ngp; ++G){
      for (J = ja-ngp; J <= jb+ngp; ++J) {
	J1 = J + (My+1+2*ngp)*(G-1);
	Un[ib+G][J] = RGLp1[J1];;
	/* printf("rank %d: RGLp1[%d] = %f \n",myRank,J1,RGLp1[J1]); */
      }
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

void printInputs() {
  // grid ngp, maxIter tol
  // ranks, , , nThreads, partition type

  printf("%s Input summary: \n",verboseTag);
  printf("    %5s %5s %4s %8s %10s \n","Nx","Ny","ngp","maxIter","tol");
  printf("    %5d %5d %4d %8d %10.3e \n",Nx,Ny,ngp,maxIter,tol);
  printf("    %5s %9s %10s \n","ranks","nThreads","partition");
  printf("    %5d %9d ",nProc,nThreads);
  if (partition == 1) {
    printf("%10s \n","square");
  } else {
    printf("%10s \n","slender");
  }
  printf(" \n");
}

void *doOneStep(void *i_thread) {
  // thread number

  // get thread id
  long thread_id = (long) i_thread;

  // get index start and end
  int iStart0 = ((int)(thread_id))*Mx_t+ia;
  int iEnd0 = iStart0+Mx_t;
  int iStart = iStart0;
  int iEnd = iEnd0;
  int k;
  /* printf("tid = %d \n",(int) thread_id); */
  /* printf("Mx_t = %d, ia = %d \n",Mx_t,ia); */
  /* printf("t = %d, beginning of thread. iStart = %d, iEnd = %d \n", */
  /* 	 thread_id,iStart0,iEnd0); */

  for (k = 0; k < ngp; ++k) {
    // handle corner cases
    if (thread_id == 0) {
      iStart = iStart0-ngp+1   -1;
    } 
    if (thread_id == nThreads-1) {
      iEnd = iEnd0+ngp-1;
    }

    // --- Jacobi --- //
    /* printf("t = %d, before jacobi \n",thread_id); */
    jacobiStep(iStart+1,iEnd);
    /* printf("t = %d, after jacobi \n",thread_id); */

    // handle corner cases again
    if (thread_id == 0) {
      iStart = iStart   -1;
    } 
    if (thread_id == nThreads-1) {
      iEnd = iEnd     +1;
    }

    /* printf("t = %d, beginning of thread. iStart = %d, iEnd = %d \n", */
    /* 	   thread_id,iStart,iEnd); */

    // --- Update BCs --- //
    updateBCs(iStart+1,iEnd);
    /* printf("t = %d, after update \n",thread_id); */


    // --- Advance --- //
    advanceGrid(iStart+1,iEnd);
    /* printf("t = %d, after advance \n",thread_id); */
  }

  pthread_exit(NULL);
}


// --- External Functions --- //
int initializeProblemInputs(int caseNumber );
double checkError(int caseNumber );

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
    printInputs();
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
    return 0;
  }

  if (verboseFlag) {
    printf("%s Problem initialization is done \n",verboseTag);
    printf("\n");
    printf("%s Begin iterations: \n",verboseTag);
  }

  // --- Iteration Loop --- //
  int n, k; // TODO, better names for these vars
  double err, maxErr;

  // thread vars  
  int rc;
  long i_thread;
  pthread_t threads[nThreads];

  // perform iterations
  for (n = 0; n < maxIter; ++n){
    /* printf("-------- n = %d ---------\n",n); */

    // perform sub iterations for the amount of
    // ghost points
    if (nThreads == 0) {
      // subiterations
      for (k = 0; k < ngp; ++k) {
	// --- Jacobi --- //
	jacobiStep(ia-ngp+1,ib+ngp-1);

	// --- Update BCs --- //
	updateBCs(ia-ngp,ib+ngp);
      
	// --- Advance --- //
	advanceGrid(ia-ngp,ib+ngp);
      
      } 
    } else {
      /* printf("rank %d: before create \n",myRank); */
      // launch threads and do subiterations
      for (i_thread = 0; i_thread < nThreads; ++i_thread) {
	rc = pthread_create(&threads[i_thread],NULL,doOneStep,(void*)i_thread);
      }
      /* printf("rank %d: after create fuck \n",myRank); */

      // handle errors
      if (rc) {
	if (myRank == 0) {
	  usage(argv[0],6);
	}
	MPI_Finalize();
	exit(0);
      }
      /* printf("rank %d: after check fuck \n",myRank); */

      // join
      for (i_thread = 0; i_thread < nThreads; ++i_thread) {
	pthread_join(threads[i_thread],NULL);
      }

      /* printf("rank %d: after join fuck \n",myRank); */
    }


    /* printf("rank %d: before collect ghost\n ",myRank); */
    // --- Collect ghost line value of Un --- //
    collectGhosts();
    /* printf("rank %d: after collect ghost\n",myRank); */
    // --- Communication --- //
    // --- exchange BBC/TBC --- //

    if( nProc > 1){

      if (myRank == 0) {
	MPI_Irecv(&(RGLp1[0]), (My+1+2*ngp)*ngp, MPI_DOUBLE, myRank+1, 1234, network,&requestRecv2); 
	MPI_Isend(&(RGL[0])  , (My+1+2*ngp)*ngp, MPI_DOUBLE, myRank+1, 1234, network,&requestSend2);
      }
      else if (myRank==nProc-1){
	MPI_Irecv(&(LGLp1[0]), (My+1+2*ngp)*ngp, MPI_DOUBLE, myRank-1, 1234, network,&requestRecv1);
	MPI_Isend(&(LGL[0])  , (My+1+2*ngp)*ngp, MPI_DOUBLE, myRank-1, 1234, network,&requestSend1);
      }
      else{
	MPI_Irecv(&(LGLp1[0]), (My+1+2*ngp)*ngp, MPI_DOUBLE, myRank-1, 1234, network,&requestRecv1);
	MPI_Irecv(&(RGLp1[0]), (My+1+2*ngp)*ngp, MPI_DOUBLE, myRank+1, 1234, network,&requestRecv2); 
	MPI_Isend(&(LGL[0])  , (My+1+2*ngp)*ngp, MPI_DOUBLE, myRank-1, 1234, network,&requestSend1);
	MPI_Isend(&(RGL[0])  , (My+1+2*ngp)*ngp, MPI_DOUBLE, myRank+1, 1234, network,&requestSend2);
      }

      if (myRank != 0) {
	// todo, check if sending is needed
	MPI_Wait(&requestSend1, &status1);
	MPI_Wait(&requestRecv1, &status1);
      }
      if (myRank != nProc-1) {
	MPI_Wait(&requestSend2, &status2);
	MPI_Wait(&requestRecv2, &status2);
      }
    }
      
    /* printf("rank %d: before update ghost \n",myRank); */
    // --- Update Ghost line value of Un --- //
    updateGhosts();
    /* printf("rank %d: after update ghost \n",myRank); */

    /* printf("rank %d: hi \n",myRank); */
    // todo figure this out
    /* MPI_Barrier(network); */

    // --- Communication --- //
    // TODO, do this

    // todo implement this

    err = calculateMaximumResidualError();
    if (verboseFlag) {
      if (printFrequency < 0) {

      } else if (n%printFrequency == 0) {
    	// todo, make a table, don't repeat n and err
    	printf("    n = %4d: err = %6.3e \n",n,err);
      }
    }

    // todo, need to consider multiple ranks 
    // barrier
    MPI_Allreduce(&err,&maxErr,1,MPI_DOUBLE,MPI_MAX,network);
    if (maxErr < tol) {
      
      break;
    }

    // communication(send and recieve parallel ghost lines)
    /* MPI_Barrier(network); */
  }

  /* if (verboseFlag) { */
  /*   printf("\n"); */
  /*   if (n == maxIter) { */
  /*     printf("%s Solver exited because the maximum ",verboseTag); */
  /*     printf("number of iterations were exceeded \n"); */
  /*     printf("    maxIter = %6d \n",maxIter); */
  /*     printf("    err     = %6.3e \n",err); */
  /*     printf("    tol     = %6.3e \n",tol); */
  /*   } else { */
  /*     printf("%s Solver exited because the residual ",verboseTag); */
  /*     printf("error is less than the tolerance \n"); */
  /*     printf("    nIter = %9d \n",n); */
  /*     printf("    err   = %9.3e \n",err);  */
  /*     printf("    tol   = %9.3e \n",tol);  */
  /*   } */
  /* } */

  double globalError, maxGlobalError;
  globalError = checkError(casenumber);
  MPI_Allreduce(&globalError,&maxGlobalError,1,MPI_DOUBLE,MPI_MAX,network);
  if (verboseFlag) {
    printf("%s globalError = %9.3e \n",verboseTag,maxGlobalError);
  } else if ((printFrequency == -2) && (myRank == 0)){

    printf("%5s %5s %4s %6s %6s %9s %10s %10s %5s \n",
  	   "Nx","Ny","ngp","nIter","nRanks","nThreads","res","err","case");
    printf("%5d %5d %4d %6d %6d %9d %10.3e %10.3e %5d \n",
  	   Nx,Ny,ngp,n,nProc,nThreads,err,maxGlobalError,casenumber);
  }
  
  MPI_Finalize();
  return 0;
}


