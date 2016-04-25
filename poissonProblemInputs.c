// --- TODO --- //
// - handle more inputs (Axis limits) in casename
// - 

#include<math.h>

typedef double** GRID;
typedef double* SUBGRID;


// --- globals --- //
// MPI
int nProc; // Number of processors
int myRank; // My rank
int nThreads; // Number of threads

// Grid
int Npx, Npy; // Number of processors in x and y dimensions
int ia; // local index of i=0
int ib; // local index of i=M

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
double forcingFunction(double x, double y);
// left, right, bottom, and top boundary conditions
double lbc(double y);
double rbc(double y);
double bbc(double x);
double tbc(double x);


// --- Main Implementation --- //
/*int initializeProblemInputs(char* casename, 
			    GRID F, 
			    SUBGRID LBC, SUBGRID RBC,
			    SUBGRID BBC, SUBGRID TBC) {*/
int initializeProblemInputs(char* casename) {
  // Choose the correct problem inputs based on the casename
  // and fill in forcing function grid, boundary conditions

  /*  inline double forcingFunction(double x, double y) { */
  /*   return 0.; */
  /* } */

  // --- Declare Inlines --- //
  /* if (strcmp(casename,"constant") == 0) { */
    // --- Constant Solution --- //
    // solution is u(x,y) = 1.
    inline double forcingFunction(double x, double y) {
      return 0.;
    }
    inline double lbc(double y) {
      return 1.;
    }
    inline double rbc(double y) {
      return 1.;
    }
    inline double bbc(double y) {
      return 1.;
    }
    inline double tbc(double y) {
      return 1.;
    }

  /* } else { */
    // --- Invalid Casename --- //
    // casename does not match any preset
  /*   return 1; */
  /* } */

  // --- Fill In Grids --- //
  // Forcing function      
  double x, y;
  for (I = ia-ngp; I <= ib+ngp; ++I) {
    for (J = ia-ngp; J <= ib+ngp; ++J) {
      // Find x and y, note that I and J are c 
      // indices that start at 0. Offset them by ia
      // to find the math index
      x = x0_ + dx*(I-ia);
      y = y0_ + dy*(J-ia);
      
      F[I][J] = forcingFunction(x,y);
    }
  }

  // BCs (left and right)
  // todo, replace with ja and jb
  for (J = ia-ngp; J <= ib+ngp; ++J) {
    y = y0_ + dy*(J-ia);
    // TODO, got a seg fault here
    /* LBC[2] = lbc(y); */
    /* LBC[2] = 1.; */
    /* RBC[2] = rbc(y); */
  }

  // BCs (bottom and top)
  /* for (I = ia-ngp; I <= ib+ngp; ++I) { */
  /*   x = x0_ + dx*(I-ia); */
  /*   BBC[I] = bbc(x); */
  /*   TBC[I] = tbc(x); */
  /* } */

  return 0;
}


