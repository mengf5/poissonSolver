// Compile with mpicc poissonSolver.c clcg4.o -lpthread
// Run with mpirun -n p ./a.out n maxIter thrds g


// --- Acitivity Log --- //
// 4/21/16 ER
// - Initial Set up of MPI and data parsing.
// - Setup the Maybe Function

#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<pthread.h>
#include "clcg4.h"
#include <time.h>

// --- Global Variables --- //
MPI_Comm ntwrk; // This is MPI_COMM_WORLD
int p; // Number of processors
int rnk; // My rank
int n; // Size of the matrix
int g; // Number of ghost cell layers to communicate
int thrds; // Number of threads
int maxIter; // Maximum number of iterations before killing the program




// --- Function Declarations --- //
void Usage(char* name, int error); // Ouput program info to the user
void Maybe(char* q, double rand); // The maybe function; the void output reminds the user that this function is too good for him






int main(int argc, char* argv[]) {

    // --- MPI Initialization --- //
    int thrd_sprt;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thrd_sprt);
    ntwrk = MPI_COMM_WORLD;
    MPI_Comm_size(ntwrk, &p);
    MPI_Comm_rank(ntwrk, &rnk);
    
    
    clock_t start = clock(), diff; // for maybe function

    
    if(rnk == 0) {
        switch(thrd_sprt) {
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
                return 0;
        }
    }

    // --- Parse input --- //
    if (argc != 5) { // Incorrect number of input arguments
        if (rnk == 0) {
            Usage(argv[0], 0);
        }
        MPI_Finalize();
        return 0;
    }

    char* ptr;
    n = strtol(argv[1], &ptr, 10); // Parse the size of the matrix
    if (*ptr != '\0' || n < 1) {
        if (rnk == 0) {
            Usage(argv[0], 1);
        }
        MPI_Finalize();
        return 0;
    }
    
    maxIter = strtol(argv[2], &ptr, 10); // Parse the number of iterations before death
    if (*ptr != '\0' || maxIter < 1) {
        if (rnk == 0) {
            Usage(argv[0], 2);
        }
        MPI_Finalize();
        return 0;
    }
    
    thrds = strtol(argv[3], &ptr, 10); // Parse the nubmer of threads
    if (*ptr != '\0' || thrds < 1) {
        if (rnk == 0) {
            Usage(argv[0],3);
        }
        MPI_Finalize();
        return 0;
    }
    
    g = strtol(argv[4], &ptr, 10); // Parse the number of ghost cell layers to communicate
    if (*ptr != '\0' || g < 1) {
        if (rnk == 0) {
            Usage(argv[0],4);
        }
        MPI_Finalize();
        return 0;
    }
    

    
    
    
    
    
    
    
    
    
    
    
    
    // --- Matrix Initialization --- //
    
    
    
    
    
    int nLayer = 0;
    int i,j;
    for(int i = 0; i < maxIter; ++i){
        for(int j = 0; j < nLayer; ++j){ 
            //caculation in each mpirank
            //(jacobi/CG/..)
        }

        MPI_Barrier(ntwrk);

           // if(err < criterion)
           //     break;  

        //communication(send and recieve parallel ghost lines)


        MPI_Barrier(ntwrk);


    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


    
    
    
    
    
    // --- For maybe function --- //
    diff = clock() - start;
    InitDefault(); // Initialize the RNG
    // --- Maybe? --- //
    if (rnk == 0) {
        // printf("%f\n%f\n", GenVal(rnk), GenVal(rnk));
        // long int guess = GenVal(rnk);
        double guess = GenVal(diff);
        //printf("%f\n" guess);
        Maybe("Will I ever be as cool as Fan Long?\n", guess);
    }



    MPI_Finalize();
    return 0;
}







// --- Functions --- //
void Usage(char* name, int error) {
    printf("Error %d in %s: \n", error, name);
    switch(error) {
        case 0:
            printf("   Incorrect number of input arguments.\n");
            printf("   Arguments should be int n, int maxIter, int thrds, int g \n");
            break;
        case 1:
            printf("   Unable to parse size of the matrix.\n");
            printf("   The size of the matrix should be a positive integer.\n");
            break;
        case 2:
            printf("   Unable to parse the maximum number of iterations.\n");
            printf("   Maximum iterations should be a positive integer.\n");
            break;
        case 3:
            printf("   Unable to parse the number of trheads.\n");
            printf("   The number of threads should be a positive integer.\n");
            break;
        case 4:
            printf("   Unable to parse the number of ghost layers.\n");
            printf("   The number of ghost layers should be a positive integer.\n");
            break;
        default:
            printf("   It is literally impossible to get here in the program.\n");
            printf("   This must be the work of the divine Maybe function!\n");
            break;
    }
}

void Maybe(char* q, double rand) {
    int out = (int) (10*rand);

    printf("%f\n", rand);
    printf("%d\n", out);

    printf("\n%s\n", q);
    switch(out) {
        case 0:
            printf("   Yes!\n");
            break;
        case 1:
            printf("   No.\n");
            break;
        case 2:
            printf("   You suck at running this program!\n");
            break;
        case 3:
            printf("   This is a dumb question.\n");
            break;
        case 4:
            printf("   Ben and Paulina are destined to be together.\n");
            break;
        case 5:
            printf("   What does the fox say?\n");
            break;
        case 6:
            printf("   Alice.\n");
            break;
        case 7:
            printf("   Why?\n");
            break;
        default:
            printf("   Maybe?\n");
            break;
    }
}