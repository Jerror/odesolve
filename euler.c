#include "ode.h"

void euler(double h, int numsteps, int dim, 
           double t0, double *u, derivative_function get_f)
{
    double f_buffer[dim]; // variable length array on the stack; requires C99
    for (int n = 0; n < numsteps; ++n;)
    {
        get_f(t0 + h * n, &u[D * n], f_buffer);
        for (int i = 0; i < dim; ++i) {
            u[dim * (n + 1) + i] += h * f_buffer[i];
        }
    }
}
