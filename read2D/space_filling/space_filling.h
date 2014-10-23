#ifndef _space_filling_h_
#define _space_filling_h_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <strings.h>

// for Hilbert
int coord_to_hilbert(int ndim, int *dims, int *coord);

// for Morton (z-order)
int coord_to_z(int ndim, int *dims, int *coord);

#endif
