/** @brief Euler method interface header
 * @author Jeremiah O'Neil
 * @copyright GNU Public License. */
#ifndef INC_EULER_h // #include guard
#define INC_EULER_h // ensure this file is only included once

#include <stdlib.h>

typedef void (*derivative_function)(double t_n, double *u_n, double *out);

void euler(double *u, double *u_init, int dim, int numsteps, 
           double h, double t, derivative_function get_f);

#endif // #include guard
