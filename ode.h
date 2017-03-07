#ifndef INC_ODE_h // #include guard
#define INC_ODE_h // ensure this file is only included once

#include <stdlib.h>

// flattened sparse representation of derivative matrix f s.th. u' = f * u
struct derivative_sparse
{
    int *j; // column index
    double *c; // element (coefficient of u_j in f_i * u)
};

typedef void *(*derivative_function)(double t, double *u_n, double *out)

#endif // #include guard
