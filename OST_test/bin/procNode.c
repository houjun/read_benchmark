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

    char filename[1024];
    char *ioData, *myData;

    double timeStat[kOST];
    double globalTimeStat[kOST];
    double *doublePtr;

    MPI_Offset disp, stripeSize;

    MPI_Status status;
    MPI_File fh;

    MPI_Comm MY_COMM_IO;
    MPI_Comm MY_COMM_SCATTER;

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
    stripeSize  = 2148473648;
    //stripeSize  = 64;

    /* printf("stripeSize: %lld\n", stripeSize); */

    // get the number of procs and rank in the comm
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

    int colorIO;
    int colorScatter;

    // color for I/O ranks
    if (worldRank % procPerNode == 0) {
        colorIO = 0;
    }
    else {
        colorIO = 1;
    }
    colorScatter = (int) (worldRank / procPerNode);

    // debug print
    /* MPI_Barrier(MPI_COMM_WORLD); */
    /* printf("World rank: %d, colorIO: %d, colorScatter: %d\n", worldRank, colorIO, colorScatter); */


    // split to I/O ranks
    MPI_Comm_split(MPI_COMM_WORLD, colorIO, worldRank, &MY_COMM_IO);
    MPI_Comm_size(MY_COMM_IO, &ioProc);
    MPI_Comm_rank(MY_COMM_IO, &ioRank);

    // split to scatter groups
    MPI_Comm_split(MPI_COMM_WORLD, colorScatter, worldRank, &MY_COMM_SCATTER);
    MPI_Comm_size(MY_COMM_SCATTER, &scatterProc);
    MPI_Comm_rank(MY_COMM_SCATTER, &scatterRank);
    
    // debug print
    /* MPI_Barrier(MPI_COMM_WORLD); */
    /* printf("World Rank: %d, IO Rank: %d, Scatter Rank:%d\n", worldRank, ioRank, scatterRank); */

     
    // data allocation
    if (worldRank % procPerNode == 0)
        ioData  = (char*)malloc(procPerNode * reqSize * sizeof(char));

    // local for all proc
    myData  = (char*)malloc(reqSize * sizeof(char));

    disp = stripeSize * worldRank;

    sprintf(filename, "%s", argv[1]);

    // debug print
    /* if (worldRank == 0) { */
    /*     printf("filename: %s\n", filename); */
    /*     printf("Rank,     time,     start offset\n"); */
    /* } */

    if (colorIO == 0) {

        MPI_File_open(MY_COMM_IO, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
        if(fh==NULL){
            if (ioRank == 0) { printf("File not exist\n"); }
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        
        MPI_File_set_view( fh, disp, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);

        // sync time
        MPI_Barrier(MY_COMM_IO);
        start_time = MPI_Wtime();

        MPI_File_read( fh, ioData, reqSize * procPerNode, MPI_CHAR, &status );
        //MPI_File_read_all( fh, ioData, reqSize * procPerNode, MPI_CHAR, &status );
 
        elapsedTime = MPI_Wtime() - start_time;

        /* printf("Rank: %d, time %.6f, start offset: %lld\n", worldRank, elapsedTime, disp); */
        
        MPI_File_close(&fh);

        // debug print
        /* doublePtr = (double*)ioData; */
        /* printf("Rank: %d, data[0]: %f\n", worldRank, *doublePtr); */
    }

    MPI_Reduce(&elapsedTime, &allTimeMax, 1, MPI_DOUBLE, MPI_MAX, 0, MY_COMM_IO);
    MPI_Reduce(&elapsedTime, &allTimeMin, 1, MPI_DOUBLE, MPI_MIN, 0, MY_COMM_IO);
    MPI_Reduce(&elapsedTime, &allTimeAvg, 1, MPI_DOUBLE, MPI_SUM, 0, MY_COMM_IO);
    allTimeAvg /= procPerNode;

    MPI_Barrier(MPI_COMM_WORLD);

    // Scatter data from I/O ranks to COMM_SCATTER
    MPI_Scatter(ioData, reqSize, MPI_CHAR, myData, reqSize, MPI_CHAR, 0, MY_COMM_SCATTER);

    doublePtr = (double*)myData;
    printf("%d: [0]%f [8]%f\n", worldRank, *doublePtr, doublePtr[7]);
    
    //int err; 
    //err = MPI_Reduce(timeStat, globalTimeStat, kOST*10, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    //if (err != MPI_SUCCESS) {
    //    printf("Error Reducing, exiting...\n");
    //    exit(-1);
    //}

    if (worldRank % procPerNode == 0) {
        free(ioData);
    }
    free(myData);

    MPI_Comm_free(&MY_COMM_IO);
    MPI_Comm_free(&MY_COMM_SCATTER);
    MPI_Finalize();

    return 0;
}

