#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
    
#include "mpi.h"
#include "space_filling/space_filling.h"

#define MAXDIM 3
#define KTIMESTEP 4 

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
    int start[MAXDIM];
    int count[MAXDIM];
    int b_dims[MAXDIM];
    int b_sizes[MAXDIM];
    int b_start[MAXDIM];
    int b_off[MAXDIM];
    int b_count[MAXDIM];
} Request_t;

void timer(MPI_Comm Comm, double elapsed_time, char* op, int my_rank, int proc_num)
{
    double all_time_max, all_time_avg, all_time_min;
    MPI_Reduce(&elapsed_time, &all_time_max, 1, MPI_DOUBLE, MPI_MAX, 0, Comm);
    MPI_Reduce(&elapsed_time, &all_time_min, 1, MPI_DOUBLE, MPI_MIN, 0, Comm);
    MPI_Reduce(&elapsed_time, &all_time_avg, 1, MPI_DOUBLE, MPI_SUM, 0, Comm);
    all_time_avg /= proc_num;
    if(my_rank == 0) 
        printf("%s time\n%f, %lf, %lf\n",  op, all_time_min, all_time_avg, all_time_max);
}

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
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    return layout_type;
}

void Print_Usage(int argc, char **argv) {
    fprintf(stderr, "Wrong number of arguments: %d\n", argc);
    fprintf(stderr, "Usage:\n%s path/to/file file_layout[XYZ,XZY,YZX,HILBERT,Z,BLOCK] ndim dim_0 ... dim_n block_size_0 ... block_size_n start_0 ... start_n count_0 ... count_n\n",argv[0]);
    MPI_Abort(MPI_COMM_WORLD, -1);
}

void Print_Request(Request_t *req)
{
    int i;
    printf("Layout: %s\n", layoutName[req->layout]);
    printf("Dim: %d\n", req->ndim);
    printf("Dim size: ");
    for (i = 0; i < req->ndim; i++) 
        printf("%d  ", req->dims[i]);
    printf("\nStart: ");
    for (i = 0; i < req->ndim; i++) 
        printf("%d  ", req->start[i]);
    printf("\nCount: ");
    for (i = 0; i < req->ndim; i++) 
        printf("%d  ", req->count[i]);
    printf("\nBlock size: ");
    for (i = 0; i < req->ndim; i++) 
        printf("%d  ", req->b_sizes[i]);
    printf("\nBlock dims: ");
    for (i = 0; i < req->ndim; i++) 
        printf("%d  ", req->b_dims[i]);
    printf("\nBlock start: ");
    for (i = 0; i < req->ndim; i++) 
        printf("%d  ", req->b_start[i]);
    printf("\nBlock count: ");
    for (i = 0; i < req->ndim; i++) 
        printf("%d  ", req->b_count[i]);
    printf("\n");
}

int Request_To_Block_Based(Request_t *req)
{
    int i;
    int head_count;
    for (i = 0; i < req->ndim; i++) {
        req->b_start[i]  = req->start[i] / req->b_sizes[i];
        req->b_off[i]    = req->start[i] % req->b_sizes[i];
        req->b_count[i]  = req->count[i] / req->b_sizes[i];
        req->b_dims[i]   = req->dims[i]  / req->b_sizes[i];

        if (req->b_off[i] != 0 || req->count[i] % req->b_sizes[i] != 0) {
            fprintf(stderr, "Request_To_Block_Based(): request region contains partial block not supported\n");
            MPI_Abort(MPI_COMM_WORLD, -1);
            // TODO: add support
            /* head_count = req->b_dims[i] - req->b_start[i] % req->b_dims[i]; */
            /* head = req->count[i] */
        }
    }

}

void Parse_Args(int argc, char **argv, Request_t *req, int my_rank)
{
    int i, t;
    if (argc < 4) {
        if (my_rank == 0) 
            Print_Usage(argc, argv);
    }
    req->layout     = Get_Layout_Type(argv[2]);
    req->ndim       = atoi(argv[3]);
    if (argc != 4+req->ndim*4) {
        if (my_rank == 0) 
            Print_Usage(argc, argv);
    }
    // Set offset after ndim
    t = 4;
    for (i = 0; i < req->ndim; i++) 
        req->dims[i] = atoi(argv[t++]);
    
    for (i = 0; i < req->ndim; i++) 
        req->b_sizes[i] = atoi(argv[t++]);
    
    for (i = 0; i < req->ndim; i++) 
        req->start[i] = atoi(argv[t++]);
    
    for (i = 0; i < req->ndim; i++) 
        req->count[i] = atoi(argv[t++]);

    // Convert to block based
    Request_To_Block_Based(req);
}

int Check_Full_Block(int ndim, int *bsizes, int *req_start)
{
    // Return 1 when the request region contains full blocks only
    // Else return -1
    int i;
    for (i = 0; i < ndim; i++) {
         if (req_start[i] % bsizes[i] != 0) {
             return -1;
         }
    }
    return 1;
}

int Get_Req_Block_Count(Request_t *req)
{
    int i, nblk;
    nblk = 1;
    for (i = 0; i < req->ndim; i++) 
        nblk *= req->b_count[i];
    return nblk;
}

int Get_Block_Size(Request_t *req)
{
    int i, size;
    size = 1;
    for (i = 0; i < req->ndim; i++) 
        size *= req->b_sizes[i];
    return size;
}

MPI_Offset Get_Timestep_Size(Request_t *req)
{
    int i;
    MPI_Offset size;
    size = 1;
    for (i = 0; i < req->ndim; i++) 
        size *= req->dims[i];
    return size;
}

MPI_Offset Get_Total_count(Request_t *req)
{
    int i;
    MPI_Offset size;
    size = 1;
    for (i = 0; i < req->ndim; i++) 
        size *= req->count[i];
    return size;
}

MPI_Offset Get_Start_offset(Request_t *req)
{
    int i, mul;
    MPI_Offset start;
    start = 0;
    mul = 1;
    for (i = req->ndim - 1; i >= 0; i--) {
        start += req->start[i] * mul;
        mul   *= req->dims[i];
    }
    return start;
}

int Get_Start_Bid(Request_t *req)
{
    int i, mul;
    int start;
    start = 0;
    mul = 1;
    for (i = req->ndim - 1; i >= 0; i--) {
        start += req->b_start[i] * mul;
        mul   *= req->b_dims[i];
    }
    return start;
}

void Set_Bcoord(Request_t *req, int* bids)
{
    int i, j, k, t, base;
    t = 0;
    if (req->ndim == 2) {
        for (i = 0; i < req->b_count[0]; i++) {
            for (j = 0; j < req->b_count[1]; j++) {
                bids[t++] = req->b_start[0] + i;
                bids[t++] = req->b_start[1] + j;
            }
        }
    }
    else if (req->ndim == 3) {
        // TODO
        ;
    }
 
}

MPI_Offset Get_Reorder_Start_offset(Request_t *req)
{
    MPI_Offset start;


    return start;
}

int main(int argc, char *argv[])
{
    int proc_num, my_rank;
    int err;
    int i, j, k, t;
    int unit_blk_size, total_req_blk_timestep;
    int *bcoords, *blocklengths, *displacements;
    double start_time, elapsed_time, all_time;
    char* fname;
    Request_t req;

    MPI_Status status;
    MPI_File fh;
    MPI_Datatype MY_READTYPE;
    MPI_Offset my_total_count, total_count, total_size, start_offset;

    // MPI init
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // Init with arguments
    fname = argv[1];
    Parse_Args(argc, argv, &req, my_rank);
    if (my_rank == 0) {
        Print_Request(&req);
    }

    // Assuming data are all double values
    total_count = Get_Total_count(&req);
    total_size  = total_count * sizeof(double);
    double* buf = (double*)malloc(total_size);

    unit_blk_size = Get_Block_Size(&req);
    total_req_blk_timestep = Get_Req_Block_Count(&req);

    bcoords      = (int*)malloc(total_req_blk_timestep * req.ndim  * sizeof(int));
    blocklengths = (int*)malloc(total_req_blk_timestep * KTIMESTEP * sizeof(int));
    displacements= (int*)malloc(total_req_blk_timestep * KTIMESTEP * sizeof(int));

    if (Check_Full_Block(req.ndim, req.b_sizes, req.start) != 1) {
        // TODO: add support when the requested region contains partial blocks
        fprintf(stderr, "Currently not support requested region contains partial blocks\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    int decompose_ratio;
    int my_off;
    // Assume ratio is 1 2 4 8 ...
    decompose_ratio = proc_num / KTIMESTEP;

    switch(req.layout) {
    case XYZ:
        unit_blk_size = 1;
        if (req.ndim == 2) {
            // Assumes 2D layout and decompose by timestep
            start_offset  = Get_Start_offset(&req);
            start_offset += (int)(my_rank / decompose_ratio) * Get_Timestep_Size(&req);

            // update actual count
            req.count[0] /= decompose_ratio;
            my_total_count= Get_Total_count(&req);

            for (i = 0; i < decompose_ratio; i++) {
                // Further decompose by row
                if (my_rank % decompose_ratio == i) {
                    my_off = i * req.count[0];
                    break;
                }
            }
            /* printf("[%d] myoff: %d\n", my_rank, my_off); */

            start_offset += my_off * req.dims[1];

            MPI_Type_vector(req.count[0], req.count[1], req.dims[1], MPI_DOUBLE, &MY_READTYPE);
            MPI_Type_commit(&MY_READTYPE);
        }
        else if (req.ndim == 3) {
            // TODO
        }
        else {
            fprintf(stderr, "Unsupported dimension, exiting...\n");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        
        break;

    case BLOCK:
    case Z:
    case HILBERT:
        // NOTE: this results in reordered result!!!
        
        Set_Bcoord(&req, bcoords);
       
        /* printf("[%d]\n", my_rank); */
        for (i = 0; i < total_req_blk_timestep * KTIMESTEP; i++) {
            blocklengths[i]   = unit_blk_size;
            if (req.layout == BLOCK) 
                displacements[i]  = coord_to_idx(req.ndim, req.b_dims, &bcoords[i*req.ndim % (total_req_blk_timestep*req.ndim)]); 
            else if (req.layout == Z) 
                displacements[i]  = coord_to_z(req.ndim, req.b_dims, &bcoords[i*req.ndim % (total_req_blk_timestep*req.ndim)]); 
            else if (req.layout == HILBERT) 
                displacements[i]  = coord_to_hilbert(req.ndim, req.b_dims, &bcoords[i*req.ndim % (total_req_blk_timestep*req.ndim)]); 

            // debug print
            /* printf("%d  ", displacements[i]); */
            displacements[i] *= unit_blk_size;
        }
 
        req.count[0] /= decompose_ratio;

        my_total_count = total_req_blk_timestep / decompose_ratio;
        if (my_rank == proc_num - 1) {
            // Probably need fix
            if (total_req_blk_timestep % decompose_ratio != 0) {
                my_total_count += total_req_blk_timestep - my_total_count * proc_num ;
            }
        }
        
        my_off         = my_rank * my_total_count; 

        start_offset   = (int)(my_rank/decompose_ratio);
        start_offset  *= Get_Timestep_Size(&req);

        /* my_total_count *= unit_blk_size; */

        /* printf("[%d] start_off: %d, my_total_count: %d, my_off: %d, disp[0]: %d, disp[1]: %d\n"                                 \ */
        /*             , my_rank, start_offset, my_total_count, my_off,  displacements[my_off], displacements[my_off+1]); */

        MPI_Type_indexed(my_total_count, blocklengths + my_off, displacements + my_off, MPI_DOUBLE, &MY_READTYPE);
        MPI_Type_commit(&MY_READTYPE);

        break;

    default:
        fprintf(stderr, "Unsupported layout type %s, exiting...\n", req.layout);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    


    MPI_Barrier(MPI_COMM_WORLD);

    // Timing
    start_time = MPI_Wtime();
    // ==========================Start of actual Read=========================
    
    err = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    assert(err == MPI_SUCCESS);

    // start_offset is in bytes
    start_offset *= sizeof(double);
    /* fprintf(stderr, "[%d]: start_offset: %d\n", my_rank, start_offset); */
    MPI_File_set_view(fh, start_offset, MPI_DOUBLE, MY_READTYPE, "native", MPI_INFO_NULL);
    MPI_File_read(fh, buf, my_total_count * unit_blk_size, MPI_DOUBLE, &status);

    MPI_File_close(&fh);

    // ============================ End of actual Read ==========================
    elapsed_time = MPI_Wtime() - start_time;
    timer(MPI_COMM_WORLD, elapsed_time, "Read", my_rank, proc_num);

    double data_in_mb = (double)(total_size*proc_num)/(1024.0*1024.0);
    /* if(my_rank == 0) */
    /*     printf("Total data: %dM Agg Bandwidth: %lf\n", (int)data_in_mb, data_in_mb/all_time); */


    // ================== == Calculate sum for varification ========================
    double my_sum, all_sum;
    uint64_t iter;
    
    /* printf("[%d]\n", my_rank); */
    /* print_data_double(req.ndim, req.count, buf); */

    start_time = MPI_Wtime();
    
    for(iter = 0; iter < my_total_count * unit_blk_size; iter++) 
        my_sum += buf[iter];

    MPI_Reduce(&my_sum, &all_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(my_rank == 0)
        printf("Sum:%lf\n", all_sum);

    // timing of reconstruction
    /* elapsed_time = MPI_Wtime() - start_time; */
    /* timer(MPI_COMM_WORLD, elapsed_time, "Summation", my_rank, proc_num); */


    MPI_Type_free(&MY_READTYPE);

    free(buf);
    free(bcoords);
    free(displacements);
    free(blocklengths);

    MPI_Finalize();
    return 0;

}
