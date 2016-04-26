// --- TODO --- //
// - handle more inputs (Axis limits) in casename
// - 

#include<math.h>
#include<string.h>
#include<stdio.h>

typedef double** GRID;
typedef double* SUBGRID;


// --- globals --- //
// MPI
int nProc; // Number of processors
int myRank; // My rank
int nThreads; // Number of threads

// Grid
int Npx, Npy; // Number of processors in x and y dimensions
int ia,ib; // local index of i=0 to i=Mx
int ja,jb; // local index of j=0 to j=My


double x0_, x1_, y0_, y1_; // local axis limits
double X0_, X1_, Y0_, Y1_; // global axis limits
double dx, dy; // local grid spacing
double dX, dY; // global grid spacing

GRID F; // forcing function
// boundary conditions
SUBGRID LBC;
SUBGRID RBC;
SUBGRID BBC;
SUBGRID TBC;

// indices
int I, J; 

// Solver
int ngp; // Number of ghost cell layers to communicate

// Flags
int verboseFlag;

// --- Function Decs --- //
// forcing used in L u = f
double forcingFunction(double x, double y, int caseNumber);
// left, right, bottom, and top boundary conditions
double lbc(double y, int caseNumber);
double rbc(double y, int caseNumber);
double bbc(double x, int caseNumber);
double tbc(double x, int caseNumber);

double forcingFunction(double x, double y, int caseNumber) {
  if( caseNumber == 0) {
    return 0.;}
  else if( caseNumber == 1){
    return 1.;}
}  
  

double lbc(double y, int caseNumber) {
  if( caseNumber == 0) {
    return 0.;}
  else if( caseNumber == 1){
    return 1.;}
}

double rbc(double y, int caseNumber) {
  if( caseNumber == 0) {
    return 0.;}
  else if( caseNumber == 1){
    return 1.;}
}

double tbc(double y, int caseNumber) {
  if( caseNumber == 0) {
    return 0.;}
  else if( caseNumber == 1){
    return 1.;}
}

double bbc(double y, int caseNumber) {
  if( caseNumber == 0) {
    return 0.;}
  else if( caseNumber == 1){
    return 1.;}
}

// --- Main Implementation --- //
/*int initializeProblemInputs(char* casename, 
			    GRID F, 
			    SUBGRID LBC, SUBGRID RBC,
			    SUBGRID BBC, SUBGRID TBC) {*/
int initializeProblemInputs(int caseNumber) {
  // Choose the correct problem inputs based on the casename
  // and fill in forcing function grid, boundary conditions

  /* inline double forcingFunction(double x, double y) { */
  /*   return 0.; */
  /* } */

  // --- Declare Inlines --- //
  /* if (strcmp(casename,"constant") == 0) { */
    // --- Constant Solution --- //
    // solution is u(x,y) = 1.
  /* inline double forcingFunction(double x, double y) { */
  /*   return 0.; */
  /* } */
  /* inline double lbc(double y) { */
  /*   return 1.; */
  /* } */
  /* inline double rbc(double y) { */
  /*   return 1.; */
  /* } */
  /* inline double bbc(double y) { */
  /*   return 1.; */
  /* } */
  /* inline double tbc(double y) { */
  /*   return 1.; */
  /* } */
  
  /* } else { */
  // --- Invalid Casename --- //
  // casename does not match any preset
  /*   return 1; */
  /* } */
  
  
  if( caseNumber == 0) {
    printf("can we do any better than this \n");
  }
  else if( caseNumber == 1){
    printf("doesnot work! go back to constant! \n");
    
  }  

  
  // --- Fill In Grids --- //
  // Forcing function      
  double x, y;
  for (I = ia-ngp; I <= ib+ngp; ++I) {
    for (J = ja-ngp; J <= jb+ngp; ++J) {
      // Find x and y, note that I and J are c 
      // indices that start at 0. Offset them by ia
      // to find the math index
      x = x0_ + dx*(I-ia);
      y = y0_ + dy*(J-ja);

      F[I][J] = forcingFunction(x,y,caseNumber);
    }
  }

  // BCs (left and right)
  for (J = ja-ngp; J <= jb+ngp; ++J) {
    y = y0_ + dy*(J-ja);
    // TODO, got a seg fault here
    LBC[J] = lbc(y,caseNumber);
    RBC[J] = rbc(y,caseNumber);
  }

  // BCs (bottom and top)
  for (I = ia-ngp; I <= ib+ngp; ++I) {
    x = x0_ + dx*(I-ia);
    BBC[I] = bbc(x,caseNumber);
    TBC[I] = tbc(x,caseNumber);
  }

  return 0;
}


