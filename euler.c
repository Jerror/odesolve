/** @brief Euler method implementation file
 * @author Jeremiah O'Neil
 * @copyright GNU Public License. */
#include "euler.h"

/** @brief Euler method
 * @param u [out] The memory to write the solution to
 * @param u_init The initial state array of the system.
 * @param dim The dimension of the system.
 * @param numsteps The number of iterations to run.
 * @param h The system parameter iteration step size.
 * @param t The initial value of the system parameter.
 * @param get_f A callback function get_f(t, *u_t, *f) which writes the
 * derivative of u at system parameter t and state u_t to array f.*/
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
