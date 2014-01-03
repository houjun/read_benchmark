#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_sort.h>

#define SIZE128K 131072
#define SIZE1M 1048576
#define SIZE16M 16777216

#define REPEAT 10
int main(int argc, char *argv[])
{
    int rank, size;
    int i, j;

    MPI_Init (&argc, &argv);  /* starts MPI */
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);  /* get current process id */
    MPI_Comm_size (MPI_COMM_WORLD, &size);    /* get number of processes */
    
    double *double_128K, *double_1M, *double_16M;
    double *double_128K_cp, *double_1M_cp, *double_16M_cp;

    double_128K    = (double*)malloc(SIZE128K);
    double_1M      = (double*)malloc(SIZE1M);
    double_16M     = (double*)malloc(SIZE16M);

    double_128K_cp = (double*)malloc(SIZE128K);
    double_1M_cp   = (double*)malloc(SIZE1M);
    double_16M_cp  = (double*)malloc(SIZE16M);

    int num_128K = SIZE128K / sizeof(double);
    int num_1M   = SIZE1M   / sizeof(double);
    int num_16M  = SIZE16M  / sizeof(double);

    // random number generator using GSL
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    for (i = 0; i < num_128K; i++) {
        double_128K[i] = gsl_rng_uniform (r);
    }
    for (i = 0; i < num_1M; i++) {
        double_1M[i] = gsl_rng_uniform (r);
    }
    for (i = 0; i < num_16M; i++) {
        double_16M[i] = gsl_rng_uniform (r);
    }
    gsl_rng_free (r);

    // make copy of array
    memcpy(double_128K_cp, double_128K, SIZE128K);
    memcpy(double_1M_cp,   double_1M,   SIZE1M  );
    memcpy(double_16M_cp,  double_16M,  SIZE16M );

    //printf("%d %d %d\n", num_128K, num_1M, num_16M);

    double res;
    double start, end, elapsed_time;

    printf("128K\n");
    // find max (linear scan)
    start = MPI_Wtime();
    for(i = 0; i < REPEAT; i++) {
        res = gsl_stats_max(double_128K, 1, num_128K);        
    }
    end = MPI_Wtime();
    elapsed_time = end - start;
    elapsed_time /= REPEAT;
    printf("%.6lf\n", elapsed_time);

    // find mean 
    start = MPI_Wtime();
    for(i = 0; i < REPEAT; i++) {
        res = gsl_stats_mean(double_128K, 1, num_128K);        
    }
    end = MPI_Wtime();
    elapsed_time = end - start;
    elapsed_time /= REPEAT;
    printf("%.6lf\n", elapsed_time);

    // find std 
    start = MPI_Wtime();
    for(i = 0; i < REPEAT; i++) {
        res = gsl_stats_sd(double_128K, 1, num_128K);        
    }
    end = MPI_Wtime();
    elapsed_time = end - start;
    elapsed_time /= REPEAT;
    printf("%.6lf\n", elapsed_time);


    // sort 
    start = MPI_Wtime();
    for(i = 0; i < REPEAT; i++) {
        gsl_sort(double_128K_cp, 1, num_128K);        
    }
    end = MPI_Wtime();
    elapsed_time = end - start;
    elapsed_time /= REPEAT;
    printf("%.6lf\n", elapsed_time);
    memcpy(double_128K_cp, double_128K, SIZE128K);


    // get 100 smallest elements 
    double smallest[100];
    start = MPI_Wtime();
    for(i = 0; i < REPEAT; i++) {
        gsl_sort_smallest(smallest, 100, double_128K_cp, 1, num_128K);        
    }
    end = MPI_Wtime();
    elapsed_time = end - start;
    elapsed_time /= REPEAT;
    printf("%.6lf\n", elapsed_time);
    memcpy(double_128K_cp, double_128K, SIZE128K);

    printf("1M\n");
    // find max (linear scan)
    start = MPI_Wtime();
    for(i = 0; i < REPEAT; i++) {
        res = gsl_stats_max(double_1M, 1, num_1M);
    }
    end = MPI_Wtime();
    elapsed_time = end - start;
    elapsed_time /= REPEAT;
    printf("%.6lf\n", elapsed_time);

    // find mean
    start = MPI_Wtime();
    for(i = 0; i < REPEAT; i++) {
        res = gsl_stats_mean(double_1M, 1, num_1M);
    }
    end = MPI_Wtime();
    elapsed_time = end - start;
    elapsed_time /= REPEAT;
    printf("%.6lf\n", elapsed_time);

    // find std
    start = MPI_Wtime();
    for(i = 0; i < REPEAT; i++) {
        res = gsl_stats_sd(double_1M, 1, num_1M);
    }
    end = MPI_Wtime();
    elapsed_time = end - start;
    elapsed_time /= REPEAT;
    printf("%.6lf\n", elapsed_time);


    // sort
    start = MPI_Wtime();
    for(i = 0; i < REPEAT; i++) {
        gsl_sort(double_1M_cp, 1, num_1M);
    }
    end = MPI_Wtime();
    elapsed_time = end - start;
    elapsed_time /= REPEAT;
    printf("%.6lf\n", elapsed_time);
    memcpy(double_1M_cp, double_1M, SIZE1M);

    // get 100 smallest elements
    start = MPI_Wtime();
    for(i = 0; i < REPEAT; i++) {
        gsl_sort_smallest(smallest, 100, double_1M_cp, 1, num_1M);
    }
    end = MPI_Wtime();
    elapsed_time = end - start;
    elapsed_time /= REPEAT;
    printf("%.6lf\n", elapsed_time);
    memcpy(double_1M_cp, double_1M, SIZE1M);
   
    printf("16M\n");
    // find max (linear scan)
    start = MPI_Wtime();
    for(i = 0; i < REPEAT; i++) {
        res = gsl_stats_max(double_16M, 1, num_16M);
    }
    end = MPI_Wtime();
    elapsed_time = end - start;
    elapsed_time /= REPEAT;
    printf("%.6lf\n", elapsed_time);

    // find mean
    start = MPI_Wtime();
    for(i = 0; i < REPEAT; i++) {
        res = gsl_stats_mean(double_16M, 1, num_16M);
    }
    end = MPI_Wtime();
    elapsed_time = end - start;
    elapsed_time /= REPEAT;
    printf("%.6lf\n", elapsed_time);

    // find std
    start = MPI_Wtime();
    for(i = 0; i < REPEAT; i++) {
        res = gsl_stats_sd(double_16M, 1, num_16M);
    }
    end = MPI_Wtime();
    elapsed_time = end - start;
    elapsed_time /= REPEAT;
    printf("%.6lf\n", elapsed_time);


    // sort
    start = MPI_Wtime();
    for(i = 0; i < REPEAT; i++) {
        gsl_sort(double_16M_cp, 1, num_16M);
    }
    end = MPI_Wtime();
    elapsed_time = end - start;
    elapsed_time /= REPEAT;
    printf("%.6lf\n", elapsed_time);
    memcpy(double_16M_cp, double_16M, SIZE16M);

    // get 100 smallest elements
    start = MPI_Wtime();
    for(i = 0; i < REPEAT; i++) {
        gsl_sort_smallest(smallest, 100, double_16M_cp, 1, num_16M);
    }
    end = MPI_Wtime();
    elapsed_time = end - start;
    elapsed_time /= REPEAT;
    printf("%.6lf\n", elapsed_time);
    memcpy(double_16M_cp, double_16M, SIZE16M);
   


    free(double_128K   );
    free(double_1M     );
    free(double_16M    );
                       
    free(double_128K_cp);
    free(double_1M_cp  );
    free(double_16M_cp );


    MPI_Finalize();
    return 0;
}
