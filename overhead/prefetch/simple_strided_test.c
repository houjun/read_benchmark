#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "mpi.h"

int main (int argc, char *argv[])
{
    
    int proc_num, my_rank, len;
    int i, j;
    double start_time, elapsed_time, all_time;
    double all_time_max, all_time_avg, all_time_min;
    MPI_Status status;
    MPI_File fh;
    MPI_Datatype contig_type;

    MPI_Init(&argc, &argv);
    
    // get the number of procs and rank in the comm
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    if(argc != 4) {
        printf("Wrong argument number!\n");
        printf("Use %s filename request_size repeat_times\n", argv[0]);
        return 0;
    }

    int req_size       = atoi(argv[2]);
    int repeat_time    = atoi(argv[3]);
    MPI_Offset stride  = proc_num * req_size;
    MPI_Offset tmp_pos = my_rank * req_size;

    char *read_data    = (char*)malloc(req_size);
     
    MPI_Type_contiguous( req_size, MPI_CHAR, &contig_type);
    MPI_Type_commit(&contig_type);


    start_time = MPI_Wtime();
    MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
    if(fh==NULL){
    	printf("File not exist\n");
    	return -1;
    }
       
    for(i = 0; i < repeat_time; i++) {
    
    //    MPI_Barrier(MPI_COMM_WORLD);
        MPI_File_read_at( fh, tmp_pos, read_data, 1, contig_type, &status );
        tmp_pos += stride;
    }
   
    MPI_File_close(&fh);

    elapsed_time = MPI_Wtime() - start_time;

    MPI_Reduce(&elapsed_time, &all_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&elapsed_time, &all_time_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&elapsed_time, &all_time_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    all_time_avg /= proc_num;

    
    double data_in_mb = (proc_num*(double)req_size*repeat_time)/(1024.0*1024.0);
    if(my_rank == 0)
        printf("Total time: %lf Min time: %lf Avg time: %lf Total data: %dM Agg Bandwidth: %lf\n", all_time, all_time_min, all_time_avg, (int)data_in_mb, data_in_mb/all_time);
    

//    printf("%d: %lf\n",my_rank, elapsed_time);
    free(read_data);
    MPI_Type_free(&contig_type);
    
    MPI_Finalize();

    return 0;
}

