#include <stdio.h>
#include <stdlib.h>
#include "space_filling.h"

// Below code from http://and-what-happened.blogspot.com/2011/08/fast-2d-and-3d-hilbert-curves-and.html
uint Morton_2D_Encode_16bit( uint index2, uint index1 )
{   // pack 2 16-bit indices into a 32-bit Morton code
    index1 &= 0x0000ffff;
    index2 &= 0x0000ffff;
    index1 |= ( index1 << 8 );
    index2 |= ( index2 << 8 );
    index1 &= 0x00ff00ff;
    index2 &= 0x00ff00ff;
    index1 |= ( index1 << 4 );
    index2 |= ( index2 << 4 );
    index1 &= 0x0f0f0f0f;
    index2 &= 0x0f0f0f0f;
    index1 |= ( index1 << 2 );
    index2 |= ( index2 << 2 );
    index1 &= 0x33333333;
    index2 &= 0x33333333;
    index1 |= ( index1 << 1 );
    index2 |= ( index2 << 1 );
    index1 &= 0x55555555;
    index2 &= 0x55555555;
    return( index1 | ( index2 << 1 ) );
}

void Morton_2D_Decode_16bit( const uint morton, uint* index1, uint* index2 )
{   // unpack 2 16-bit indices from a 32-bit Morton code
    uint value1 = morton;
    uint value2 = ( value1 >> 1 );
    value1 &= 0x55555555;
    value2 &= 0x55555555;
    value1 |= ( value1 >> 1 );
    value2 |= ( value2 >> 1 );
    value1 &= 0x33333333;
    value2 &= 0x33333333;
    value1 |= ( value1 >> 2 );
    value2 |= ( value2 >> 2 );
    value1 &= 0x0f0f0f0f;
    value2 &= 0x0f0f0f0f;
    value1 |= ( value1 >> 4 );
    value2 |= ( value2 >> 4 );
    value1 &= 0x00ff00ff;
    value2 &= 0x00ff00ff;
    value1 |= ( value1 >> 8 );
    value2 |= ( value2 >> 8 );
    value1 &= 0x0000ffff;
    value2 &= 0x0000ffff;
    *index1 = value1;
    *index2 = value2;
}

uint Morton_3D_Encode_10bit( uint index3, uint index2, uint index1 )
{   // pack 3 10-bit indices into a 30-bit Morton code
    index1 &= 0x000003ff;
    index2 &= 0x000003ff;
    index3 &= 0x000003ff;
    index1 |= ( index1 << 16 );
    index2 |= ( index2 << 16 );
    index3 |= ( index3 << 16 );
    index1 &= 0x030000ff;
    index2 &= 0x030000ff;
    index3 &= 0x030000ff;
    index1 |= ( index1 << 8 );
    index2 |= ( index2 << 8 );
    index3 |= ( index3 << 8 );
    index1 &= 0x0300f00f;
    index2 &= 0x0300f00f;
    index3 &= 0x0300f00f;
    index1 |= ( index1 << 4 );
    index2 |= ( index2 << 4 );
    index3 |= ( index3 << 4 );
    index1 &= 0x030c30c3;
    index2 &= 0x030c30c3;
    index3 &= 0x030c30c3;
    index1 |= ( index1 << 2 );
    index2 |= ( index2 << 2 );
    index3 |= ( index3 << 2 );
    index1 &= 0x09249249;
    index2 &= 0x09249249;
    index3 &= 0x09249249;
    return( index1 | ( index2 << 1 ) | ( index3 << 2 ) );
}

void Morton_3D_Decode_10bit( const uint morton, uint* index3, uint* index2, uint* index1 )
{   // unpack 3 10-bit indices from a 30-bit Morton code
    uint value1 = morton;
    uint value2 = ( value1 >> 1 );
    uint value3 = ( value1 >> 2 );
    value1 &= 0x09249249;
    value2 &= 0x09249249;
    value3 &= 0x09249249;
    value1 |= ( value1 >> 2 );
    value2 |= ( value2 >> 2 );
    value3 |= ( value3 >> 2 );
    value1 &= 0x030c30c3;
    value2 &= 0x030c30c3;
    value3 &= 0x030c30c3;
    value1 |= ( value1 >> 4 );
    value2 |= ( value2 >> 4 );
    value3 |= ( value3 >> 4 );
    value1 &= 0x0300f00f;
    value2 &= 0x0300f00f;
    value3 &= 0x0300f00f;
    value1 |= ( value1 >> 8 );
    value2 |= ( value2 >> 8 );
    value3 |= ( value3 >> 8 );
    value1 &= 0x030000ff;
    value2 &= 0x030000ff;
    value3 &= 0x030000ff;
    value1 |= ( value1 >> 16 );
    value2 |= ( value2 >> 16 );
    value3 |= ( value3 >> 16 );
    value1 &= 0x000003ff;
    value2 &= 0x000003ff;
    value3 &= 0x000003ff;
    *index1 = value1;
    *index2 = value2;
    *index3 = value3;
}

uint Morton_to_Hilbert2D( const uint morton, const uint bits )
{
    uint hilbert = 0;
    uint remap = 0xb4;
    uint block = ( bits << 1 );
    while( block )
    {
        block -= 2;
        uint mcode = ( ( morton >> block ) & 3 );
        uint hcode = ( ( remap >> ( mcode << 1 ) ) & 3 );
        remap ^= ( 0x82000028 >> ( hcode << 3 ) );
        hilbert = ( ( hilbert << 2 ) + hcode );
    }
    return( hilbert );
}

uint Hilbert_to_Morton2D( const uint hilbert, const uint bits )
{
    uint morton = 0;
    uint remap = 0xb4;
    uint block = ( bits << 1 );
    while( block )
    {
        block -= 2;
        uint hcode = ( ( hilbert >> block ) & 3 );
        uint mcode = ( ( remap >> ( hcode << 1 ) ) & 3 );
        remap ^= ( 0x330000cc >> ( hcode << 3 ) );
        morton = ( ( morton << 2 ) + mcode );
    }
    return( morton );
}

uint Morton_to_Hilbert3D( const uint morton, const uint bits )
{
    uint hilbert = morton;
    if( bits > 1 )
    {
        uint block = ( ( bits * 3 ) - 3 );
        uint hcode = ( ( hilbert >> block ) & 7 );
        uint mcode, shift, signs;
        shift = signs = 0;
        while( block )
        {
            block -= 3;
            hcode <<= 2;
            mcode = ( ( 0x20212021 >> hcode ) & 3 );
            shift = ( ( 0x48 >> ( 7 - shift - mcode ) ) & 3 );
            signs = ( ( signs | ( signs << 3 ) ) >> mcode );
            signs = ( ( signs ^ ( 0x53560300 >> hcode ) ) & 7 );
            mcode = ( ( hilbert >> block ) & 7 );
            hcode = mcode;
            hcode = ( ( ( hcode | ( hcode << 3 ) ) >> shift ) & 7 );
            hcode ^= signs;
            hilbert ^= ( ( mcode ^ hcode ) << block );
        }
    }
    hilbert ^= ( ( hilbert >> 1 ) & 0x92492492 );
    hilbert ^= ( ( hilbert & 0x92492492 ) >> 1 );
    return( hilbert );
}

uint Hilbert_to_Morton3D( const uint hilbert, const uint bits )
{
    uint morton = hilbert;
    morton ^= ( ( morton & 0x92492492 ) >> 1 );
    morton ^= ( ( morton >> 1 ) & 0x92492492 );
    if( bits > 1 )
    {
        uint block = ( ( bits * 3 ) - 3 );
        uint hcode = ( ( morton >> block ) & 7 );
        uint mcode, shift, signs;
        shift = signs = 0;
        while( block )
        {
            block -= 3;
            hcode <<= 2;
            mcode = ( ( 0x20212021 >> hcode ) & 3 );
            shift = ( ( 0x48 >> ( 4 - shift + mcode ) ) & 3 );
            signs = ( ( signs | ( signs << 3 ) ) >> mcode );
            signs = ( ( signs ^ ( 0x53560300 >> hcode ) ) & 7 );
            hcode = ( ( morton >> block ) & 7 );
            mcode = hcode;
            mcode ^= signs;
            mcode = ( ( ( mcode | ( mcode << 3 ) ) >> shift ) & 7 );
            morton ^= ( ( hcode ^ mcode ) << block );
        }
    }
    return( morton );
}

// Above code copied and slightly modified

int is_power_of_2(int n)
{
    return n == 1 || (n & (n-1)) == 0;
}

int get_power_of_2(int n)
{
    int count = 0;
    while (n > 1) {
        n = n >> 1;
        count++;
    }
    return count;
}

int dim_check(int ndim, int *dims)
{
    // currently both z and hilbert convertion is restricted to equal dim size across all dims
    // and dim size must be power of 2
    int i;
    for (i = 0; i < ndim; i++) {
        // check if any is zero
        if (dims[i] == 0) {
            fprintf(stderr, "dim_check() error: dims[%d]=%d\n", i, dims[i]);
            return -1;
        }

        if (is_power_of_2(dims[i]) != 1) {
            fprintf(stderr, "dim_check() error: dims[%d]=%d is not power of 2\n", i, dims[i]);
            return -1;
        }

        // check if all dim size is equal
        if (i > 0) {
            if (dims[i] != dims[i-1]) {
                fprintf(stderr, "dim_check() error: dims[%d]=%d is not equal to dims[%d]=%d\n", i, dims[i], i-1, dims[i-1]);
                return -1;
            }
        }
    }
    
    return 1;
}

int coord_to_z(int ndim, int *dims, int* coord)
{
    if (dim_check(ndim, dims) == -1) {
        fprintf(stderr, "coord_to_z() error\n");
        return -1; 
    }

    if (ndim == 2) {
        // check the coord are within 16bit
        if (coord[0] > 65536 || coord[1] > 65536) {
            fprintf(stderr, "coord2z(2D) error, current coord (%d, %d) exceeds supported 16-bit\n", coord[0], coord[1]);
            return -1;
        }
        else {
            return (int)Morton_2D_Encode_16bit((uint)coord[0], (uint)coord[1]);
        }

    }
    else if (ndim == 3) {
        // check the coord are within 16bit
        if (coord[0] > 1024 || coord[1] > 1024) {
            fprintf(stderr, "coord2z(3D) error, current coord (%d, %d, %d) exceeds supported 10-bit\n", coord[0], coord[1], coord[2]);
            return -1;
        }
        else {
            return (int)Morton_3D_Encode_10bit((uint)coord[0], (uint)coord[1], (uint)coord[2]);
        }
    }
    else {
        fprintf(stderr,"coord2z() error, current dim :%d is not supported\n", ndim);
        return -1;
    }

    return -1;
}

int coord_to_hilbert(int ndim, int *dims, int* coord)
{
    if (dim_check(ndim, dims) == -1) {
        fprintf(stderr, "coord_to_hilbert() error\n");
        return -1; 
    }

    // get z order index first
    int zid;
    zid = coord_to_z(ndim, dims, coord);

    int size = get_power_of_2(dims[0]);

    if (ndim == 2) {
        return (int)Morton_to_Hilbert2D((uint)zid, (uint)size);
    }
    else if (ndim == 3) {
        return (int)Morton_to_Hilbert2D((uint)zid, (uint)size);
    }

    return -1;
}

int coord_to_idx(int ndim, int *dims, int* coord)
{
    // linearize the coord to pos
    int i, pos, mul;

    pos = 0;
    mul = 1;
    for (i = ndim - 1; i >= 0; i--) {
        pos += coord[i] * mul;
        mul *= dims[i];
    }   

    return pos; 
}

int idx_to_coord(int ndim, int *dims, int idx, int* coord)
{
    int i, mod;
    int *dimSize;

    dimSize = (int*)malloc(ndim*sizeof(int));
    dimSize[0] = dims[0];

    // total number of elements
    for (i = 1; i < ndim; i++) {
        dimSize[i] = dimSize[i-1] * dims[i];
    }


    mod = idx;
    for (i = 0; i < ndim-1; i++) {
        coord[i] = mod / dimSize[ndim-i-2];
        mod      = mod % dimSize[ndim-i-2];
    }
    coord[ndim-1] = mod % dimSize[ndim-i-1];

    free(dimSize);
}

void test_idx_to_coord(int ndim, int *dims)
{
    int i;
    int coord[3];
    int tmp, size;
    size = 1;
    // total number of elements
    for (i = 0; i < ndim; i++) {
        size *= dims[i];
    }

    for (i = 0; i < size; i++) {
        idx_to_coord(ndim, dims, i, coord);
        tmp = coord_to_idx(ndim, dims, coord);

        if (ndim == 3) {
            printf("%d => (%d, %d, %d) => %d\n", i, coord[0], coord[1], coord[2], tmp);
        }
        else if (ndim == 2) {
            printf("%d => (%d, %d) => %d\n", i, coord[0], coord[1], tmp);
        }
        if (i != tmp) {
            fprintf(stderr, "ERROR!!\n");
        }
    }
}

void print_data_double(int ndim, int *dims, double *data)
{
    int i, j;
    int *dimSize;

    dimSize = (int*)malloc(ndim*sizeof(int));
    dimSize[0] = dims[0];

    // total number of elements
    for (i = 1; i < ndim; i++) {
        dimSize[i] = dimSize[i-1] * dims[i];
    }

    printf("\n");
    for (i = 0; i < dimSize[ndim-1]; i++) {
        printf("%4.1f  ", data[i]);   

        for (j = 0; j < ndim; j++) {
            if ( (i+1) % dimSize[j] == 0) {
                printf("\n");
            }
        }
    }

    free(dimSize);
}


void print_data_int(int ndim, int *dims, int *data)
{
    int i, j;
    int *dimSize;

    dimSize = (int*)malloc(ndim*sizeof(int));
    dimSize[0] = dims[0];

    // total number of elements
    for (i = 1; i < ndim; i++) {
        dimSize[i] = dimSize[i-1] * dims[i];
    }

    printf("\n");
    for (i = 0; i < dimSize[ndim-1]; i++) {
        printf("%4d  ", data[i]);   

        for (j = 0; j < ndim; j++) {
            if ( (i+1) % dimSize[j] == 0) {
                printf("\n");
            }
        }
    }

    free(dimSize);
}

int test(int ndim, int *dims)
{
    int i, j, k;
    int count, size, pos;
    int coord[3];

    int *data_ori, *data_z, *data_h;


    size = 1;
    // total number of elements
    for (i = 0; i < ndim; i++) {
        size *= dims[i];
    }

    data_ori = (int*)malloc(size * sizeof(int));
    data_z   = (int*)malloc(size * sizeof(int));
    data_h   = (int*)malloc(size * sizeof(int));

    // init
    count = 0;
    for (i = 0; i < size; i++)
        data_ori[i] = count++;

    // print original data
    print_data_int(ndim ,dims, data_ori);


    // transform to z-ordered data
    for (i = 0; i < size; i++) {
        idx_to_coord(ndim, dims, i, coord);

        pos = coord_to_z(ndim, dims, coord);
        data_z[pos] = data_ori[i];
    }

    printf("Z-ordered:\n");
    print_data_int(ndim, dims, data_z);
    
    
    // transform to hilbert-ordered data
    for (i = 0; i < size; i++) {
        idx_to_coord(ndim, dims, i, coord);

        pos = coord_to_hilbert(ndim, dims, coord);
        data_h[pos] = data_ori[i];
    }
    
    printf("Hilbert-ordered:\n");
    print_data_int(ndim, dims, data_h);
}
/*
int main(int argc, char* argv[])
{
    int ndim, i;
    int dims[3];

    ndim = atoi(argv[1]);
    for (i = 0; i < ndim; i++) {
        dims[i] = atoi(argv[2+i]);
    }

    test_idx_to_coord(ndim, dims);
    
    test(ndim, dims);


    return 0;
}
*/
