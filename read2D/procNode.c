#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "mpi.h"

//const int kOST=156; // Hopper
const int kOST         =144; // Edison
const int kCorePerNode = 24;    // NERSC

int main (int argc, char *argv[])
{
    int procNum, worldRank, len;
    int ioProc, ioRank, scatterProc, scatterRank;
    int i, j, k;
    int ioProcPerNode, dbg_print;
    int myTotalDoubleNum;

    char filename[1024];
    double *ioData, *myData;

    double start_time, elapsedTime, scatterTime, allTimeMax, allTimeAvg, allTimeMin;
    double timeStat[kOST];
    double globalTimeStat[kOST];
    double mySum;
    double *doublePtr;
    
    int colorIO;
    int colorScatter;
    int scatterGroupSize;

    MPI_Offset disp, myReqSize, stripeSize, ioReadSize;
    MPI_Datatype strided;

    MPI_Status status;
    MPI_File fh;

    MPI_Comm MY_COMM_IO;
    MPI_Comm MY_COMM_SCATTER;

    // Init
    MPI_Init(&argc, &argv);
       
    if(argc < 3) {
        if (worldRank == 0) {
            fprintf(stderr, "Wrong argument number!\n");
            fprintf(stderr, "Usage:\n%s filepath single_read_size(MB) I/O_proc_per_node\n", argv[0]);
        }
        MPI_Finalize();
        exit(-1);
    }
 
    sprintf(filename, "%s", argv[1]);
    myReqSize     = atoi(argv[2]);
    myReqSize    *= 1048576;
    ioProcPerNode = atoi(argv[3]);
    stripeSize    = 2147483648;     //2G

    // adjust to double from char
    myReqSize    /= sizeof(double);
    stripeSize   /= sizeof(double);
    // It turns out MPI seems to have a problem when stride size is greakter than 2148473647 with vector datatype,
    // which is then used in MPI_File_set_view, the first process reads less data than expected.
    // Current solution is to use basic type MPI_DOUBLE instead of MPI_BYTE, so the stride is 1/8 or original.

    scatterGroupSize = kCorePerNode / ioProcPerNode;
    ioReadSize       = myReqSize * scatterGroupSize;

    MPI_Type_vector(scatterGroupSize, myReqSize, stripeSize, MPI_DOUBLE, &strided);
    MPI_Type_commit(&strided);
   
    /* printf("stripeSize: %lld\n", stripeSize); */

    // get the number of procs and rank in the comm
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

    // color for I/O ranks and scatter ranks
    colorIO      = worldRank % scatterGroupSize == 0 ? 0 : MPI_UNDEFINED;
    colorScatter = (int) (worldRank / scatterGroupSize);

    // debug print
    /* MPI_Barrier(MPI_COMM_WORLD); */
    /* printf("World rank: %d, colorIO: %d, colorScatter: %d\n", worldRank, colorIO, colorScatter); */

    // split to I/O ranks
    MPI_Comm_split(MPI_COMM_WORLD, colorIO, worldRank, &MY_COMM_IO);
    if (colorIO == 0) {
        MPI_Comm_size(MY_COMM_IO, &ioProc);
        MPI_Comm_rank(MY_COMM_IO, &ioRank);
    }

    // split to scatter groups
    MPI_Comm_split(MPI_COMM_WORLD, colorScatter, worldRank, &MY_COMM_SCATTER);
    MPI_Comm_size(MY_COMM_SCATTER, &scatterProc);
    MPI_Comm_rank(MY_COMM_SCATTER, &scatterRank);
    
    // debug print
    /* printf("World Rank: %d, IO Rank: %d, Scatter Rank:%d\n", worldRank, ioRank, scatterRank); */

    // local for all proc
    myData  = (double*)malloc(myReqSize * sizeof(double));
    if (myData == NULL) {
        fprintf(stderr, "Error allocating myData\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // actual I/O workers
    if (colorIO == 0) {
        //debug print
        /* printf("[%d] %dMB\n", worldRank, (int)(ioReadSize / 1048576)); */

        // I/O process needs more space
        ioData  = (double*)malloc(ioReadSize * sizeof(double));
        if (ioData == NULL) {
            fprintf(stderr, "Error allocating ioData\n");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        // distribute the data read among I/O ranks
        disp = stripeSize * worldRank;

        MPI_File_open(MY_COMM_IO, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
        if(fh==NULL){
            if (worldRank == 0) { printf("File not exist\n"); }
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        
        MPI_File_set_view( fh, disp, MPI_DOUBLE, strided, "native", MPI_INFO_NULL);


        MPI_Barrier(MY_COMM_IO);
        start_time = MPI_Wtime();

        MPI_File_read( fh, ioData, ioReadSize, MPI_DOUBLE, &status );
 
        elapsedTime = MPI_Wtime() - start_time;
        
        MPI_File_close(&fh);

    }

    MPI_Barrier(MPI_COMM_WORLD);

    start_time = MPI_Wtime();

    // Scatter data from I/O ranks to COMM_SCATTER
    MPI_Scatter(ioData, myReqSize, MPI_DOUBLE, myData, myReqSize, MPI_DOUBLE, 0, MY_COMM_SCATTER);


    scatterTime = MPI_Wtime() - start_time;

    // data correctness verification
    doublePtr = (double*)myData;
    myTotalDoubleNum = myReqSize / sizeof(double);
    mySum = 0.0;
    for (i = 0; i < myTotalDoubleNum; i++) {
        mySum += doublePtr[i];
    }
    printf("[%d] %f\n", worldRank, mySum);

    // Barrier to make sure the timing are not messed up with above prints
    MPI_Barrier(MPI_COMM_WORLD);

    if (colorIO == 0) {
        printf("Rank: %d, time %.6f, start offset: %lld, size: %llu\n", worldRank, elapsedTime, disp, ioReadSize);
    }

    // stat for read time among I/O ranks
    if (colorIO == 0) {
        MPI_Reduce(&elapsedTime, &allTimeMax, 1, MPI_DOUBLE, MPI_MAX, 0, MY_COMM_IO);
        MPI_Reduce(&elapsedTime, &allTimeMin, 1, MPI_DOUBLE, MPI_MIN, 0, MY_COMM_IO);
        MPI_Reduce(&elapsedTime, &allTimeAvg, 1, MPI_DOUBLE, MPI_SUM, 0, MY_COMM_IO);
        allTimeAvg /= ioProcPerNode;

        if (worldRank == 0) {
            printf("[Read time] %f, %f, %f\n", allTimeMin, allTimeAvg, allTimeMax);
            printf("[Scatter time] %f\n", scatterTime);
        }
    }
   
   
    if (colorIO == 0) {
        free(ioData);
        MPI_Comm_free(&MY_COMM_IO);
    }

    free(myData);
    MPI_Comm_free(&MY_COMM_SCATTER);

    MPI_Type_free(&strided);
    MPI_Finalize();

    return 0;
}

