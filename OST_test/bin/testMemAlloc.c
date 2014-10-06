#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main(int argc, char *argv[])
{
    int myRank, numProc;

    MPI_Init(&argc, &argv);

    size_t size, i;
    size_t sizeMB = atoi(argv[1]);
    size = sizeMB * 1048576;

    char* t;
    t = (char*)malloc(size*sizeof(char));
    if (t == NULL) {
        fprintf(stderr, "Error allocating %dMB...\n", sizeMB);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    else {
        for (i = 0; i < size; i++) {
            t[i] = 'a';
        }
        printf("Allocation for %dMB success\n", sizeMB);
    }


    MPI_Finalize();
    
    return 0;
}
