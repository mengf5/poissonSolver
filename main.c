// Compile with mpicc.mpich main.c -lpthread
// Run with mpirun.mpich -n p ./a.out n thrds g

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




// --- Function Declarations --- //
void Usage(char* name, int error); // Ouput program info to the user
void Maybe(char* q, double rand); // The maybe function; the void output reminds the user that this function is too good for him






int main(int argc, char* argv[]) {

    // --- MPI Initialization --- //
    int thrd_sprt;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &thrd_sprt);
    ntwrk = MPI_COMM_WORLD;
    MPI_Comm_size(ntwrk, &p);
    MPI_Comm_rank(ntwrk, &rnk);
    
    
    clock_t start = clock(), diff;

    
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
    if (argc != 4) { // Incorrect number of input arguments
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

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


    // --- Matrix Initialization --- //
    
    
    diff = clock() - start;
    InitDefault(); // Initialize the RNG
    // --- Maybe? --- //
    if (rnk == 0) {
        // printf("%f\n%f\n", GenVal(rnk), GenVal(rnk));
        // long int guess = GenVal(rnk);
        double guess = GenVal(diff);
        //printf("%f\n" guess);
        Maybe("Who is Henshaw's Favorite?\n", guess);
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
            printf("   Arguments should be int n, int thrds, int g \n");
            break;
        case 1:
            printf("   Unable to parse size of the matrix.\n");
            printf("   The size of the matrix should be a positive integer.\n");
            break;
        default:
            printf("   The Maybe function is the best function in the world.\n");
            printf("   Also, you suck at running this program.\n");
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