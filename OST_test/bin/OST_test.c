#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "mpi.h"

//const int kOST=156; // Hopper
const int kOST=144; // Edison

int main (int argc, char *argv[])
{
    
    int proc_num, my_rank, len;
    int i, j;
    int req_size, repeat_time, stride, dbg_print;
    double start_time, elapsed_time;
    double all_time_max, all_time_avg, all_time_min;
    char* read_data;

    MPI_Status status;
    MPI_File fh;
    MPI_Datatype contig_type;
    MPI_Datatype stride_type;
    MPI_Offset OST_proc, start_pos, stripe_size, total_size_proc;

    MPI_Init(&argc, &argv);
    
    // get the number of procs and rank in the comm
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    if(argc < 6) {
        if (my_rank == 0) {
            printf("Wrong argument number!\n");
            printf("Use %s filename single_read_size total_size_per_proc stripe_size OST_per_proc\n", argv[0]);
        }
        MPI_Finalize();
        return 0;
    }

    req_size        = atoi(argv[2]);
    total_size_proc = atoi(argv[3]);
    stripe_size     = atoi(argv[4]);
    OST_proc        = atoi(argv[5]);

    if (argc == 7) {
        dbg_print = atoi(argv[6]);
    }

    repeat_time     = (int)(total_size_proc / req_size);
    start_pos       = (my_rank/OST_proc) * stripe_size;
    stride          = (int)(kOST * stripe_size / req_size);


    read_data       = (char*)malloc(total_size_proc*sizeof(char));
     
    MPI_Type_contiguous( req_size, MPI_CHAR, &contig_type);
    MPI_Type_commit(&contig_type);

    MPI_Type_vector(repeat_time, 1, stride, contig_type, &stride_type);
    MPI_Type_commit(&stride_type);

    MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if(fh==NULL){
        if (my_rank == 0) {
    	    printf("File not exist\n");
        }
    	return -1;
    }

    printf("req size:%d, total_size_proc:%d, stripe size:%d, OST_proc:%d, repeat_time:%d, stride:%d\n",req_size, total_size_proc, stripe_size, OST_proc, repeat_time, stride);
    printf("%d: start:%d\n", my_rank, start_pos);
    
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    MPI_File_read_at( fh, start_pos, read_data, 1, stride_type, &status );
 
    elapsed_time = MPI_Wtime() - start_time;
    if(dbg_print==1)
        printf("%d %lf\n", my_rank, elapsed_time);
    
    MPI_Reduce(&elapsed_time, &all_time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&elapsed_time, &all_time_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&elapsed_time, &all_time_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    
    double data_in_mb = (proc_num*(double)req_size*repeat_time)/(1024.0*1024.0);
    if(my_rank == 0) {
        all_time_avg /= proc_num;
        printf("%lf %lf %lf\n", all_time_min, all_time_avg, all_time_max);
    }

  

    MPI_File_close(&fh);


    free(read_data);
    MPI_Type_free(&contig_type);
    MPI_Type_free(&stride_type);
    
    MPI_Finalize();

    return 0;
}

