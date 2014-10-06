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
    
    int proc_num, my_rank, len;
    int i, j, k;
    int req_size, repeat_time, stride, dbg_print;
    double start_time, elapsed_time;
    double all_time_max, all_time_avg, all_time_min;
    char* read_data;
    char filename[256];

    double timeStat[kOST][10];
    double globalTimeStat[kOST][10];

    MPI_Status status;
    MPI_File fh;
    MPI_Datatype contig_type;
    MPI_Datatype stride_type;
    MPI_Offset OST_proc, start_pos, stripe_size, total_size_proc;

    MPI_Init(&argc, &argv);
    
    // get the number of procs and rank in the comm
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    if(argc < 5) {
        if (my_rank == 0) {
            printf("Wrong argument number!\n");
            printf("Use %s filepath single_read_size(MB) repeat_time\n", argv[0]);
        }
        MPI_Finalize();
        exit(-1);
    }

    req_size    = atoi(argv[2]);
    repeat_time = atoi(argv[3]);
    OST_proc    = atoi(argv[4]);

    req_size   *= (1048576/OST_proc);    // Convert to byte

    /* total_size_proc = atoi(argv[3]); */
    /* stripe_size     = atoi(argv[4]); */

    /* if (argc == 7) { */
    /*     dbg_print = atoi(argv[6]); */
    /* } */

    /* repeat_time     = (int)(total_size_proc / req_size); */
    /* start_pos       = (my_rank/OST_proc) * stripe_size; */
    /* stride          = (int)(kOST * stripe_size / req_size); */

    read_data       = (char*)malloc(req_size*sizeof(char));
     
    MPI_Type_contiguous( req_size, MPI_CHAR, &contig_type);
    MPI_Type_commit(&contig_type);

    /* MPI_Type_vector(repeat_time, 1, stride, contig_type, &stride_type); */
    /* MPI_Type_commit(&stride_type); */

    MPI_Offset disp;

    for (j = 0; j < repeat_time; j++) {
        disp = j * req_size * (my_rank+1);
        
        /* for (i = 0; i < kOST; i++) { */
        for (i = 0; i < 1; i++) {
        
            sprintf(filename, "%s/%d/temp1.bin", argv[1], my_rank);
            /* sprintf(filename, "%s/%d/temp.bin", argv[1], i); */
            /* printf("%s\n", filename); */

            if (OST_proc == 1) {
                MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
            }
            else {
                MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
            }
            if(fh==NULL){
                if (my_rank == 0) { printf("File not exist\n"); }
                exit(-1);
            }
            
            /* MPI_File_set_view( fh, disp, MPI_CHAR, contig_type, "native", MPI_INFO_NULL); */
            MPI_Barrier(MPI_COMM_WORLD);
            start_time = MPI_Wtime();

            /* printf("disp:%lld, size:%d\n", disp, req_size); */
            MPI_File_read( fh, read_data, req_size, MPI_CHAR, &status );
            /* MPI_File_read_all( fh, read_data, req_size, MPI_CHAR, &status ); */
         
            elapsed_time = MPI_Wtime() - start_time;
            timeStat[i][j] = elapsed_time;
            printf("Rank %d, reading from %s, start %d, time %.6f\n", my_rank, filename, disp, elapsed_time);
            
            MPI_File_close(&fh);
        }
    }

    /* MPI_Reduce(&elapsed_time, &all_time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); */
    /* MPI_Reduce(&elapsed_time, &all_time_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD); */
    /* MPI_Reduce(&elapsed_time, &all_time_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); */
    
    int err; 
    err = MPI_Reduce(timeStat, globalTimeStat, kOST*10, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (err != MPI_SUCCESS) {
        printf("Error Reducing, exiting...\n");
        exit(-1);
    }

/*     if (my_rank == 0) { */
/*         for (j = 0; j < kOST; j++) { */
                
/*             printf("%d", j); */
/*             for (i = 0; i < repeat_time; i++) { */
/*                  printf(", %.6f", globalTimeStat[j][i]); */
/*             } */
/*             printf("\n"); */
/*         } */
/*     } */
   
   
    free(read_data);
    MPI_Type_free(&contig_type);
    /* MPI_Type_free(&stride_type); */
    
    MPI_Finalize();

    return 0;
}

