#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
    
#include "mpi.h"
#include "space_filling/space_filling.h"

#define MAXDIM 3

typedef enum {
    XYZ,
    XZY,
    ZYX,
    BLOCK,
    Z,
    HILBERT
} Layout_t;

const char layoutName[][10] = {"XYZ", "XZY", "ZYX", "BLOCK", "Z", "HILBERT"};

typedef struct request {
    Layout_t layout;
    int ndim;
    int dims[MAXDIM];
    int blk_size[MAXDIM];
    int start[MAXDIM];
    int count[MAXDIM];
} Request_t;

Layout_t Get_Layout_Type(char *type)
{
    Layout_t layout_type;

    if(strcmp(type, "XYZ") == 0){
        layout_type = XYZ;
    }
    else if(strcmp(type, "XZY") == 0){
        layout_type = XZY;
    }    
    else if(strcmp(type, "ZYX") == 0){
        layout_type = ZYX;
    }    
    else if(strcmp(type, "BLOCK") == 0){
        layout_type = BLOCK;
    }    
    else if(strcmp(type, "HILBERT") == 0){
        layout_type = HILBERT;
    }    
    else if(strcmp(type, "Z") == 0){
        layout_type = Z;
    }    
    else {
        fprintf(stderr, "Unsupported layout type %s, exiting...\n", type);
        exit(-1);
    }

    return layout_type;
}

void Print_Usage(int argc, char **argv) {
    fprintf(stderr, "Wrong number of arguments: %d\n", argc);
    fprintf(stderr, "Usage:\n%s path/to/file file_layout[XYZ,XZY,YZX,HILBERT,Z,BLOCK] ndim dim_0 ... dim_n block_size_0 ... block_size_n start_0 ... start_n count_0 ... count_n\n",argv[0]);
    exit(-1);
}

void Print_Request(Request_t *req)
{
    int i;
    printf("Layout: %s\n", layoutName[req->layout]);
    printf("Dim: %d\n", req->ndim);
    printf("Dim size: ");
    for (i = 0; i < req->ndim; i++) 
        printf("%d  ", req->dims[i]);
    printf("\nBlock size: ");
    for (i = 0; i < req->ndim; i++) 
        printf("%d  ", req->blk_size[i]);
    printf("\nStart: ");
    for (i = 0; i < req->ndim; i++) 
        printf("%d  ", req->start[i]);
    printf("\nCount: ");
    for (i = 0; i < req->ndim; i++) 
        printf("%d  ", req->count[i]);
    printf("\n");
}

void Parse_Args(int argc, char **argv, Request_t *req)
{
    int i, t;
    req->layout     = Get_Layout_Type(argv[2]);
    req->ndim       = atoi(argv[3]);
    if (argc != 4+req->ndim*4) {
        Print_Usage(argc, argv);
    }
    // set offset after ndim
    t = 4;
    for (i = 0; i < req->ndim; i++) 
        req->dims[i] = atoi(argv[t++]);
    
    for (i = 0; i < req->ndim; i++) 
        req->blk_size[i] = atoi(argv[t++]);
    
    for (i = 0; i < req->ndim; i++) 
        req->start[i] = atoi(argv[t++]);
    
    for (i = 0; i < req->ndim; i++) 
        req->count[i] = atoi(argv[t++]);
}



int main(int argc, char *argv[])
{
    int proc_num, my_rank;
    int err;
    int i, j, k, t;
    double start_time, elapsed_time, all_time;
    char* fname;
    Request_t req;

    // MPI init
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // init with arguments
    fname = argv[1];
    Parse_Args(argc, argv, &req);
    Print_Request(&req);

    // Assuming data are all double values




   /*
    // assume cubic space
    // use fix block size for now
    int cubic_size = block_size * block_size * block_size;
    int nblock_x   = x / block_size; 
    int total_blk  = nblock_x * nblock_x * nblock_x;
    int nblock_xy  = nblock_x * nblock_x;
    
    int* used_blk;
    int* mapping_c2h;
    int* mapping_h2c;

    int  blk_coor_x;
    int  blk_coor_y;
    int  nblock_need;
    int  prev_blk_num;
    
    int* array_of_blocklengths;
    int* array_of_displacements;

    switch(layout_type) {
        case XYZ:
            start_offset = subplane_start_y*x + subplane_start_x + my_rank*x*y + plane_id*x*y; 
            
            if(subplane_start_x == 0 && subplane_end_x == y-1) {
                // read entire plane
                MPI_Type_contiguous(x*y, unit_datatype, &final_datatype);
                MPI_Type_commit(&final_datatype);
            }
            else {
                // read partial plane
                MPI_Type_vector(subplane_end_y-subplane_start_y+1, subplane_end_x-subplane_start_x+1, x, unit_datatype, &final_datatype);
                MPI_Type_commit(&final_datatype);
            }

            // dummy strided_datatype
            MPI_Type_vector(subplane_end_y-subplane_start_y+1, 1, y, unit_datatype, &stride_datatype);
            MPI_Type_commit(&stride_datatype);
            break;

        case XZY:
            start_offset = subplane_start_y*x*z + subplane_start_x + my_rank*x + plane_id*x; 
            
            MPI_Type_vector(subplane_end_y-subplane_start_y+1, subplane_end_x-subplane_start_y+1, x*z, unit_datatype, &final_datatype);
            MPI_Type_commit(&final_datatype);

            // dummy strided_datatype
            MPI_Type_vector(subplane_end_y-subplane_start_y+1, 1, y, unit_datatype, &stride_datatype);
            MPI_Type_commit(&stride_datatype);
            break;

        case ZYX:
            // NOTE: this results in reordered (transposed) result!!!
            start_offset = subplane_start_x*y*z + subplane_start_y*z + my_rank + plane_id; 
            
            // assuming single sub-plane now
            MPI_Type_vector(subplane_end_y-subplane_start_y+1, 1, y, unit_datatype, &stride_datatype);
            MPI_Type_commit(&stride_datatype);

            // tried to used negative stride, but doesn't work
            // also using MPI_Type_Vector over MPI_Type_Vector datatype doesn't work as expected
            // so use MPI_Type_create_hvector over MPI_Type_Vector, which does work.
            MPI_Type_create_hvector(subplane_end_x-subplane_start_x+1, 1, y*z*unit_size, stride_datatype, &final_datatype);
            MPI_Type_commit(&final_datatype);

            break;

        case BLOCK:
            // NOTE: this results in reordered result!!!
            start_offset = subplane_start_x / block_size * cubic_size + subplane_start_y / block_size * nblock_x * cubic_size + (plane_id/block_size+my_rank)*x*y*block_size;

            // assume accessed plane width it the divisible by block size
            // read entire plane of each block
            MPI_Type_vector( (subplane_end_x-subplane_start_x)/block_size + 1, block_size*block_size, cubic_size,
                             unit_datatype, &stride_datatype);

            MPI_Type_commit(&stride_datatype);


            MPI_Type_create_hvector( (subplane_end_y-subplane_start_y)/block_size + 1, 1, y*z*unit_size,
                                     stride_datatype, &final_datatype);
            
            MPI_Type_commit(&final_datatype);

            break;

        case HILBERT:

            // transform from x-y-z to Hilbert curve
            mapping_c2h = (int*)malloc(total_blk*sizeof(int));
            mapping_h2c = (int*)malloc(total_blk*sizeof(int));
            hilbert_init(3, nblock_x, mapping_c2h, mapping_h2c);

            used_blk = (int*)malloc(sizeof(int)*total_blk);
            memset(used_blk, 0, sizeof(int)*total_blk);
            nblock_need = 0;
            for(i = 0; i < nblock_xy; i++) {

                // calculate block coordinate and compare with the requested range
                blk_coor_x = (i % nblock_x) * block_size;
                blk_coor_y = (int)(i / nblock_x) * block_size;

                if( blk_coor_x >= subplane_start_x && blk_coor_x <=subplane_end_x &&
                    blk_coor_y >= subplane_start_y && blk_coor_y <=subplane_end_y   ) {

                    // mark the corresponding block number in Hilbert curve
                    used_blk[mapping_c2h[i+(plane_id/block_size)*nblock_xy]] = 1;
                    nblock_need++;
                }
            }
            

            array_of_blocklengths  = (int*)malloc(nblock_need*sizeof(int));
            array_of_displacements = (int*)malloc(nblock_need*sizeof(int));

            t = 0;
            prev_blk_num = -1;
            for(i = 0; i < total_blk; i++) {
                if(used_blk[i] == 1) {

                    // calculate initial offset
                    if(t == 0) {
                        start_offset = i * cubic_size;
                        prev_blk_num = i;
                        array_of_displacements[0] = 0;
                    }
                    else {
                        array_of_displacements[t] = (i-prev_blk_num) * cubic_size + array_of_displacements[t-1];
                    }
                    // all has the same block length with the assumption that
                    // the requested range contains COMPLETE x-y plane of all needed blocks
                    array_of_blocklengths[t++] = block_size * block_size;
                    prev_blk_num = i;
                }

            }

            if(nblock_need == 1) {
                MPI_Type_contiguous(block_size*block_size, unit_datatype, &final_datatype);
                MPI_Type_commit(&final_datatype);
            }
            else {
                MPI_Type_indexed(nblock_need, array_of_blocklengths, array_of_displacements, unit_datatype, &final_datatype);
                MPI_Type_commit(&final_datatype);
            }

            free(used_blk);
            free(mapping_c2h);
            free(mapping_h2c);

            // dummy commit just to make error of freeing this go away
            MPI_Type_vector( (subplane_end_x-subplane_start_x)/block_size + 1, block_size*block_size, cubic_size,
                             unit_datatype, &stride_datatype);
            MPI_Type_commit(&stride_datatype);
            break;

    }
    

    MPI_Status status;
    MPI_File fh;
    size_t total_size = (subplane_end_x-subplane_start_x+1)*(subplane_end_y-subplane_start_y+1)*unit_size;
    char* buf = (char*)malloc(total_size);

    MPI_Barrier(MPI_COMM_WORLD);

    // Timing
    start_time = MPI_Wtime();
    // ==========================Start of actual Read=========================
    
    err = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    assert(err == MPI_SUCCESS);

    MPI_File_set_view(fh, start_offset*unit_size, unit_datatype, final_datatype, "native", MPI_INFO_NULL);
    MPI_File_read(fh, buf, total_size, MPI_BYTE, &status);

    MPI_File_close(&fh);

    // ============================ End of actual Read ==========================
    elapsed_time = MPI_Wtime() - start_time;
    all_time = timer(elapsed_time, "Read");

    double data_in_mb = (double)(total_size*proc_num)/(1024.0*1024.0);
    if(my_rank == 0)
        printf("Total data: %dM Agg Bandwidth: %lf\n", (int)data_in_mb, data_in_mb/all_time);


    // ================== == Calculate sum for varification ========================
    double my_sum, all_sum;
    uint64_t nelem = (subplane_end_x-subplane_start_x+1)*(subplane_end_y-subplane_start_y+1);
    uint64_t iter;
    
    start_time = MPI_Wtime();
    
    for(iter = 0; iter <nelem; iter++) {
        my_sum += data_full[iter];
    }
    MPI_Reduce(&my_sum, &all_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(my_rank == 0)
        printf("Sum:%lf\n", all_sum);

    // timing of reconstruction
    elapsed_time = MPI_Wtime() - start_time;
    all_time = timer(elapsed_time, "Sum");

        
    MPI_Type_free(&unit_datatype);
    MPI_Type_free(&stride_datatype);
    MPI_Type_free(&final_datatype);

    if(total_size_full != 0)
        free(data_full);
    free(buf);

    */

    MPI_Finalize();
    return 0;
}
