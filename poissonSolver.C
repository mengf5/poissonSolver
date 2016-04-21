//********* Revision Summary *************

//********* Includes *********************
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>
#include<mpi.h>
#include<clcg4.h>
#include<mpi.h>
#include<pthread.h>

//*************** Usage *****************
void usage(char* name, int errNo) {
  printf("Error %d in %s \n",errNo,name);
  switch (errNo) {
  case 0:

    break;
  default:

    break;
  }
}

//*************** MAIN ******************
int main(int argc, char** argv){

  //initialization
  // since the calculation is grid free, i guess we just need initial guess and rhs

  for(int i=0; i<nIteration; ++i){
    for(int j=0; j<nLayer;++j){ 
      //caculation in each mpirank
      //(jacobi/CG/..)
    }

    MPI_Barrier 

      if(err < criterion)
	break;  

    //communication(send and recieve parallel ghost lines)


    MPI_Barrier 


      }

  return 0
    };
