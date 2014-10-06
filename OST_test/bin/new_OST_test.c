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
    int i, j;
    int req_size, repeat_time, stride, dbg_print;
    double start_time, elapsed_time;
    double all_time_max, all_time_avg, all_time_min;
    char* read_data;
    char filename[256];

    double timeStat[145][10];

    MPI_Status status;
    MPI_File fh;
    MPI_Datatype contig_type;
    MPI_Datatype stride_type;
    MPI_Offset OST_proc, start_pos, stripe_size, total_size_proc;

    MPI_Init(&argc, &argv);
    
    // get the number of procs and rank in the comm
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    if(argc < 4) {
        if (my_rank == 0) {
            printf("Wrong argument number!\n");
            printf("Use %s filepath single_read_size(MB) repeat_time\n", argv[0]);
        }
        MPI_Finalize();
        exit(-1);
    }

    req_size    = atoi(argv[2]);
    req_size   *= 1048576;    // Convert to byte
    repeat_time = atoi(argv[3]);
    /* total_size_proc = atoi(argv[3]); */
    /* stripe_size     = atoi(argv[4]); */
    /* OST_proc        = atoi(argv[5]); */

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

    
    for (j = 0; j < repeat_time; j++) {
        start_pos = j*req_size;
        
        for (i = 0; i < kOST; i++) {
        
            sprintf(filename, "%s/%d/temp.bin", argv[1], i);
            /* printf("%s\n", filename); */

            MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
            if(fh==NULL){
                if (my_rank == 0) {
                    printf("File not exist\n");
                }
                exit(-1);
            }
            
            start_time = MPI_Wtime();

            MPI_File_read_at( fh, start_pos, read_data, 1, contig_type, &status );
         
            elapsed_time = MPI_Wtime() - start_time;

            timeStat[i][j] = elapsed_time;
            /* printf("%d, %lf\n", i, elapsed_time); */
            
            MPI_File_close(&fh);

            //sleep(1);
        
        }
    }

    /* MPI_Reduce(&elapsed_time, &all_time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); */
    /* MPI_Reduce(&elapsed_time, &all_time_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD); */
    /* MPI_Reduce(&elapsed_time, &all_time_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); */
     

    for (j = 0; j < kOST; j++) {
            
        printf("%d", j);
        for (i = 0; i < repeat_time; i++) {
             printf(", %.6f", timeStat[j][i]);
        }
        printf("\n");
    }
   
    
    /* double data_in_mb = (proc_num*(double)req_size*repeat_time)/(1024.0*1024.0); */
    /* if(my_rank == 0) { */
    /*     all_time_avg /= proc_num; */
    /*     printf("%lf %lf %lf\n", all_time_min, all_time_avg, all_time_max); */
    /* } */

  



    free(read_data);
    MPI_Type_free(&contig_type);
    /* MPI_Type_free(&stride_type); */
    
    MPI_Finalize();

    return 0;
}

