#include "ode.h"

void euler(double *u, double *u_init, int dim, int numsteps, 
           double h, double t, derivative_function get_f)
{
    double f_prev[dim]; // variable length array buffer on stack - requires C99
    double *u_prev = u_init;
    for (int n = 0; n < numsteps; ++n)
    {
        get_f(t, u_prev, f_prev);
        for (int i = 0; i < dim; ++i) {
            u[i] = u_prev[i] + h * f_prev[i];
        }
        t += h; u_prev = u; u += dim; // advance time and pointers
    }
}
