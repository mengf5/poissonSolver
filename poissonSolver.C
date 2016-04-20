#include 

int main(int argc, char** argv){

for(int i=0; i<nIteration; ++i){
  
for(int j=0; j<nLayer;++j){ 
%%caculation in each mpirank
(jacobi/CG/..)
}

MPI_Barrier 

%%communication(send and recieve parallel ghost lines)


MPI_Barrier 
}

return 0
};
