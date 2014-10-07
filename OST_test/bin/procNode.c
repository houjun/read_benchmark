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
    double start_time, elapsedTime, scatterTime, allTimeMax, allTimeAvg, allTimeMin;

    char filename[1024];
    char *ioData, *myData;

    double timeStat[kOST];
    double globalTimeStat[kOST];
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
            printf("Wrong argument number!\n");
            printf("Use %s filepath single_read_size(MB) I/O_proc_per_node\n", argv[0]);
        }
        MPI_Finalize();
        exit(-1);
    }
 
    // filename = argv[1];
    myReqSize     = atoi(argv[2]);
    myReqSize    *= 1048576;
    ioProcPerNode = atoi(argv[3]);
    stripeSize    = 2147483648;
    
    scatterGroupSize = kCorePerNode / ioProcPerNode;
    ioReadSize       = myReqSize * scatterGroupSize;
    //stripeSize  = 64;

    MPI_Type_vector(scatterGroupSize, myReqSize, stripeSize, MPI_CHAR, &strided);
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
    printf("World rank: %d, colorIO: %d, colorScatter: %d\n", worldRank, colorIO, colorScatter);


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

     
    // data allocation
    // I/O process needs more space
    if (colorIO == 0) {
        //debug print
        printf("%d: %dMB\n", worldRank, (int)(ioReadSize / 1048576));
        ioData  = (char*)malloc(ioReadSize * sizeof(char));
        if (ioData == NULL) {
            fprintf(stderr, "Error allocating ioData\n");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }

    // local for all proc
    myData  = (char*)malloc(myReqSize * sizeof(char));
    if (myData == NULL) {
        fprintf(stderr, "Error allocating myData\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    /* printf("%d: %dMB\n", worldRank, (int)(myReqSize / 1048576)); */

    disp = stripeSize * worldRank;

    sprintf(filename, "%s", argv[1]);

    if (colorIO == 0) {

        MPI_File_open(MY_COMM_IO, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
        if(fh==NULL){
            if (ioRank == 0) { printf("File not exist\n"); }
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        
        MPI_File_set_view( fh, disp, MPI_CHAR, strided, "native", MPI_INFO_NULL);


        MPI_Barrier(MY_COMM_IO);
        start_time = MPI_Wtime();

        MPI_File_read( fh, ioData, ioReadSize, MPI_CHAR, &status );
 
        elapsedTime = MPI_Wtime() - start_time;
        
        MPI_File_close(&fh);

        // debug print
        doublePtr = (double*)ioData;
        printf("Rank: %d, data[0]: %f, disp: %lld, size: %lld\n", worldRank, *doublePtr, disp, ioReadSize);
    }

    // debug print
    /* if (worldRank == 0) { */
    /*     doublePtr = (double*)ioData; */
    /*     printf("Data[%lld]%f\n", myReqSize/8, doublePtr[myReqSize/8]); */
    /*     printf("Data[%lld]%f\n", myReqSize/8 - 1, doublePtr[(myReqSize)/8 - 1]); */
    /* } */
    

    start_time = MPI_Wtime();

    // Scatter data from I/O ranks to COMM_SCATTER
    MPI_Scatter(ioData, myReqSize, MPI_CHAR, myData, myReqSize, MPI_CHAR, 0, MY_COMM_SCATTER);

    MPI_Barrier(MPI_COMM_WORLD);

    scatterTime = MPI_Wtime() - start_time;

    // data correctness verification
    doublePtr = (double*)myData;
    int myTotalDoubleNum;
    double mySum;
    myTotalDoubleNum = myReqSize / sizeof(double);
    mySum = 0.0;
    for (i = 0; i < myTotalDoubleNum; i++) {
        mySum += doublePtr[i];
    }
    printf("[%d] %f\n", worldRank, mySum);

    /* if (colorIO == 0) { */
    /*     printf("Rank: %d, time %.6f, start offset: %lld, size: %llu\n", worldRank, elapsedTime, disp, ioReadSize); */
    /* } */

    // stat for read time among I/O ranks
    if (colorIO == 0) {
        MPI_Reduce(&elapsedTime, &allTimeMax, 1, MPI_DOUBLE, MPI_MAX, 0, MY_COMM_IO);
        MPI_Reduce(&elapsedTime, &allTimeMin, 1, MPI_DOUBLE, MPI_MIN, 0, MY_COMM_IO);
        MPI_Reduce(&elapsedTime, &allTimeAvg, 1, MPI_DOUBLE, MPI_SUM, 0, MY_COMM_IO);
        allTimeAvg /= ioProcPerNode;

        if (worldRank == 0) {
            printf("Read time\n%f, %f, %f\n", allTimeMin, allTimeAvg, allTimeMax);
            printf("Scatter time\n%f\n", scatterTime);
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

