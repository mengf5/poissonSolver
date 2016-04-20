//This is my understanding on what we need to do based on today's meeting
//look very similar to the Game Of Life hw to me.

#include whatever.h

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
