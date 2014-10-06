#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "mpi.h"

//const int kOST=156; // Hopper
const int kOST=144; // Edison

int main (int argc, char *argv[])
{
    
    int procNum, worldRank, len;
    int ioProc, ioRank, scatterProc, scatterRank;
    int i, j, k;
    int reqSize, procPerNode, dbg_print;
    double start_time, elapsedTime, allTimeMax, allTimeAvg, allTimeMin;

    char* ioData, myData;

    double timeStat[kOST];
    double globalTimeStat[kOST];

    MPI_Status status;
    MPI_File fh;

    MPI_Comm COMMIO;
    MPI_Comm COMMSCATTER;

    // Init
    MPI_Init(&argc, &argv);
       
    if(argc < 3) {
        if (worldRank == 0) {
            printf("Wrong argument number!\n");
            printf("Use %s filepath single_read_size(MB) proc_per_node\n", argv[0]);
        }
        MPI_Finalize();
        exit(-1);
    }
 
    // filename = argv[1];
    reqSize     = atoi(argv[2]);
    procPerNode = atoi(argv[3]);

    // get the number of procs and rank in the comm
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

    int colorIO;
    int colorScatter;

    // color for I/O ranks
    for (i = 0; i < procNum; i++) {
        if (i % procPerNode == 0) {
            colorIO = 0;
        }
        else {
            colorIO = 1;
        }
        colorScatter = (int) (i / procPerNode);
    }

    // split to I/O ranks
    MPI_Comm_split(MPI_COMM_WORLD, colorIO, worldRank, &COMMIO);
    MPI_Comm_size(COMMIO, &ioProc);
    MPI_Comm_rank(COMMIO, &ioRank);

    // split to scatter groups
    MPI_Comm_split(MPI_COMM_WORLD, colorScatter, worldRank, &COMMSCATTER);
    MPI_Comm_size(COMMSCATTER, &scatterProc);
    MPI_Comm_rank(COMMSCATTER, &scatterRank);
    
    // debug print
    printf("World Rank: %d, IO Rank: %d, Scatter Rank:%d\n", worldRank, ioRank, scatterRank);
     
    if (worldRank % procPerNode == 0) {
        ioData  = (char*)malloc(procPerNode * reqSize * sizeof(char));
    }

    /*
    myData  = (char*)malloc(procPerNode * reqSize * sizeof(char));

    MPI_Offset disp;

    disp = 0;
    disp = 2147483648 * worldRank;
    
    for (i = 0; i < 1; i++) {
    
        sprintf(filename, "%s", argv[1]);
        //printf("%s\n", filename);

        MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

        if(fh==NULL){
            if (worldRank == 0) { printf("File not exist\n"); }
            exit(-1);
        }
        
        MPI_File_set_view( fh, disp, MPI_CHAR, contig_type, "native", MPI_INFO_NULL);
        MPI_Barrier(MPI_COMM_WORLD);
        start_time = MPI_Wtime();

        MPI_File_read( fh, myData, reqSize, MPI_CHAR, &status );
        //MPI_File_read_all( fh, myData, reqSize, MPI_CHAR, &status );
     
        elapsedTime = MPI_Wtime() - start_time;
        timeStat[i]  = elapsedTime;
        printf("Rank %d, reading from %s, start %lld, time %.6f\n", worldRank, filename, disp, elapsedTime);
        
        MPI_File_close(&fh);
    }

    //MPI_Reduce(&elapsedTime, &allTimeMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    //MPI_Reduce(&elapsedTime, &allTimeMin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    //MPI_Reduce(&elapsedTime, &allTimeAvg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    int err; 
    err = MPI_Reduce(timeStat, globalTimeStat, kOST*10, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (err != MPI_SUCCESS) {
        printf("Error Reducing, exiting...\n");
        exit(-1);
    }

    free(myData);
    
   */

    MPI_Comm_free(&COMMIO);
    MPI_Comm_free(&COMMSCATTER);
    MPI_Finalize();

    return 0;
}

