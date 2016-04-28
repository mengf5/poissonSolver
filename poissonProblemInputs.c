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
GRID Un; 
GRID Unp1;
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
char* verboseTag;

// Case specific constants
double a2 = 1.3;
double b2 = 1.9;

// --- Function Decs --- //
// forcing used in L u = f
double forcingFunction(double x, double y, int caseNumber) {
  if( caseNumber == 0) {
    return 0.;
  } else if( caseNumber == 1) {
    return 0.;
  } else if (caseNumber == 2) {
    return -(a2*a2+b2*b2)*cos(a2*x)*cos(b2*y);
  }
}  
  
// left, right, bottom, and top boundary conditions
double lbc(double y, int caseNumber) {
  if( caseNumber == 0) {
    return 1.;
  } else if( caseNumber == 1) {
    return exp(2*x0_)*sin(2*y);
  } else if (caseNumber == 2) {
    return cos(a2*x0_)*cos(b2*y);
  }

}

double rbc(double y, int caseNumber) {
  if( caseNumber == 0) {
    return 1.;
  } else if( caseNumber == 1) {
    return exp(2*x1_)*sin(2*y);
  } else if (caseNumber == 2) {
    return cos(a2*x1_)*cos(b2*y);
  }
}

double tbc(double x, int caseNumber) {
  if( caseNumber == 0) {
    return 1.;
  } else if( caseNumber == 1) {
    return exp(2*x)*sin(2*y1_);
  } else if (caseNumber == 2) {
    return cos(a2*x)*cos(b2*y1_);
  }
}

double bbc(double x, int caseNumber) {
  if( caseNumber == 0) {
    return 1.;
  } else if( caseNumber == 1) {
    return exp(2*x)*sin(2*y0_);
  } else if (caseNumber == 2) {
    return cos(a2*x)*cos(b2*y0_);
  }
}

double uTrue(double x, double y, int caseNumber) {
  if (caseNumber == 0) {
    return 1.;
  } else if (caseNumber == 1) {
    return exp(2*x)*sin(2*y);
  } else if (caseNumber == 2) {
    return cos(a2*x)*cos(b2*y);
  }
}

// --- Check Error --- ///
double checkError(int caseNumber) {
  // if available, check the max error
  
  // casenumber does not have an 'exact' solution
  if ((caseNumber < 0) || (caseNumber > 2)) {
    return -1.;
  }

  double err = 0.;
  double check, x, y;

  for (I = ia; I <= ib; ++I) {
    for (J = ja; J <= jb; ++J) {
      x = x0_ + dx*(I-ia);
      y = y0_ + dy*(J-ja);

      // use residual of equation
      check = fabs(Unp1[I][J]-uTrue(x,y,caseNumber));
      /* printf("(%d,%d): | %9.3e - %9.3e | = %9.3e \n",I-ia,J-ja, */
      /* 	     Unp1[I][J],uTrue(x,y,caseNumber),check); */
      if (check > err) {
	err = check;
      }
    }
  }

  return err;
}

// --- Main Implementation --- //
int initializeProblemInputs(int caseNumber) {
  // Choose the correct problem inputs based on the casename
  // and fill in forcing function grid, boundary conditions

  // casenumber is invalid
  if ((caseNumber < 0) || (caseNumber > 2)) {
    return 1;
  }
  // print casename
  if (verboseFlag) {
    printf("%s ",verboseTag);
    if (caseNumber == 0) {
      printf("casename: constant \n");
    } else if (caseNumber == 1) {
      printf("casename: harmonic function \n");
    } else if (caseNumber == 2) {
      printf("casename: trig solution \n");
    } else {
      return 1;
    }  
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



