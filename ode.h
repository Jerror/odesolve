#ifndef INC_ODE_h // #include guard
#define INC_ODE_h // ensure this file is only included once

#include <stdlib.h>

typedef void (*derivative_function)(double t_n, double *u_n, double *out);

void euler(double h, int numsteps, int dim, 
           double t0, double *u, derivative_function get_f);

#endif // #include guard
