#ifndef INC_ODE_h // #include guard
#define INC_ODE_h // ensure this file is only included once

#ifdef __cplusplus // if this header was included in a C++ context
extern "C" { // use C linkage. Forbids symbol mangling (and thus overloading)
#endif // This makes it possible to cross-link C and C++ objects.

#include <stdlib.h>

typedef void (*derivative_function)(double t_n, double *u_n, double *out);

void euler(double *u, double *u_init, int dim, int numsteps, 
           double h, double t, derivative_function get_f);

#ifdef __cplusplus
} // closing brace for extern "C"
#endif

#endif // #include guard
