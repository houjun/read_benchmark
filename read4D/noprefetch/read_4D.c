/*
 * read_4D.c
 *
 *  Created on: Sep 25, 2013
 *      Author: houjun
 */

/*
    ROW decompose
	    ____________
    	   /           /|  ...Proc0
    	  /           //|
    	 /___________// |z
myrows I |__________|/  |
    	 |          |  /
    	 |          | /y
    	 |          |/    ...ProcN
    	 ------------
      	       x
    	      T0
*/


#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>

typedef enum {
    ROW,
    COL,
    CUBE,
} decomp_t;

void usage() {
    printf( "Usage:\nread_4D FILENAME Decompose_Type X Y Z B Time_step Wait_time Time_modifier [T_replay_start] [T_replay_end]\n");
    exit(1);
}

int main(int argc, char *argv[])
{
    int proc_num, my_rank;
    int b, x, y, z, time_step, t_replay_start, t_replay_end;
    int i, j, k, t, m, n, tmp, err;
    struct timespec ts;
    double start_time, elapsed_time, all_time, ht_wait_time, ht_modifier;
    double all_time_max, all_time_avg, all_time_min;
    decomp_t decompose_type;

    MPI_Status status;
    MPI_Datatype contig;
    MPI_File fh;

    // MPI init
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // check arguments
    if (argc < 9) {
        printf("%d\n",argc);
        usage();
    }
    // init
    char *fname = argv[1];

    if(strcmp(argv[2], "ROW") == 0){
        decompose_type = ROW;
//        if(my_rank == 0)
//            printf("ROW ");
    }
    else if(strcmp(argv[2], "COL") == 0){
        decompose_type = COL;
//        if(my_rank == 0)
//            printf("COL ");
    }
    else if(strcmp(argv[2], "CUBE") == 0){
        decompose_type = CUBE;
//        if(my_rank == 0)
//            printf("CUBE ");
    }
    else{
        printf("Unsupported decompose type: %s\n", argv[2]);
        return 0;
    } 

    x = atoi(argv[3]);	        	// number of rows of cubic
    y = atoi(argv[4]);
    z = atoi(argv[5]);
    b = atoi(argv[6]);

    time_step    = atoi(argv[7]);		// start time step
    ht_wait_time = atof(argv[8]);
    ht_modifier  = atof(argv[9]);

    if(argc == 12) {
        t_replay_start = atoi(argv[10]); // "replay" start time step
        t_replay_end   = atoi(argv[11]);
    }
    else {
        t_replay_start = 0; // "replay" start time step
        t_replay_end   = time_step;
    }

    if(my_rank == 0)
        printf("x:%d y:%d z:%d b:%d time_step:%d wait_time:%lf multiplier:%lf \n",
               x, y, z, b, time_step, ht_wait_time, ht_modifier);
    /*
    int sleep;
    if(my_rank == 0) {
        sleep=1;
        while(sleep){;}
    }
    */
    
    MPI_Offset start_offset = 0;

    // distribute work to different procs
    int my_read_cnt      = 0;
    int my_one_read_size = 0;
    int my_total_size_xy = 0;
    int my_total_size    = 0;

    switch(decompose_type){
        case ROW:
            my_read_cnt      = 1;
            my_one_read_size = x * b * y / proc_num; 
            my_total_size_xy = my_one_read_size;
            my_total_size    = my_total_size_xy * z;
            break;

        case COL:
            my_read_cnt      = y;
            my_one_read_size = x * b / proc_num; 
            my_total_size_xy = my_one_read_size * y;
            my_total_size    = my_total_size_xy * z;
            break;

        case CUBE:
            ;
            // TODO
    }
   
    MPI_Type_contiguous(my_one_read_size, MPI_BYTE, &contig);
    MPI_Type_commit(&contig);
    
    // allocate buffer for reading data for each time step
    char *buf = (char*)malloc(my_total_size);
    assert(buf != NULL);

    // wait time structure
    ht_wait_time *= ht_modifier;
    ts.tv_sec = (int)ht_wait_time;
    ts.tv_nsec = (ht_wait_time - ts.tv_sec) * 1000000000.0;

    // start together
    MPI_Barrier(MPI_COMM_WORLD);
    

    // Open file
    err = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    assert(err == MPI_SUCCESS);

    start_time = MPI_Wtime();

    if(decompose_type == ROW) {

        for(t = 0; t < time_step; t++) {
    
            // each time step start offset is x*y*z
            start_offset = t * x * y * z * b;
    
            // read a slice of each time step
            for(i = 0; i < z; i++) {
                 MPI_File_read_at(fh, (start_offset + my_rank * my_total_size_xy)
                                     , buf + i*my_total_size_xy, 1, contig, &status);
                 start_offset += (x * b * y);
                 nanosleep(&ts, NULL);
            }
    
        }

    }
    else if (decompose_type == COL) {
       
        for(t = 0; t < time_step; t++) {
                
            // each time step start offset is x*y*z
            start_offset = t * x * y * z * b;
            
            for(j = 0; j < z; j++) {
            
                for(i = 0; i < y; i++) {
                    
                    MPI_File_read_at(fh, (start_offset + my_rank * my_one_read_size + i*x*b + j*x*b*y)
                                      , buf + i*my_one_read_size, 1, contig, &status);
                    nanosleep(&ts, NULL);

                }
                

            }
        }

    }
    else if (decompose_type == CUBE) {

    }
     
    
    elapsed_time = MPI_Wtime() - start_time;
    err = MPI_File_close(&fh);
    assert(err == MPI_SUCCESS);
    

    MPI_Reduce(&elapsed_time, &all_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&elapsed_time, &all_time_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&elapsed_time, &all_time_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    all_time_avg /= proc_num;

    double data_in_mb = ((double)(x)*y*z*b*time_step)/(1024.0*1024.0);
    if(my_rank == 0)
        printf("Total time: %lf Min time: %lf Avg time: %lf Total data: %dM Agg Bandwidth: %lf Wait Time%lf\n"
                , all_time, all_time_min, all_time_avg, (int)data_in_mb, data_in_mb/all_time, ts.tv_sec+ts.tv_nsec/1000000000.0);




    /*
    // check read numbers
    if(my_rank == 1){
    	int cnt = 0;
    	for(i = 0; i < b*x*myrows*z*(t_end-t_start); i++){
    		if(i % (x*b) == 0)
    			printf("\n");
    		if(i % (x*b*myrows) == 0)
    			printf("\n==============%d============\n\n",cnt++);
    		printf(" %3d",buf[i]);
    	}
    	printf("\n");
    }

    */
    free(buf);
    MPI_Type_free(&contig);

    MPI_Finalize();
    return 0;
}


