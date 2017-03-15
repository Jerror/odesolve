#include "../euler.h"
#include "../adaptive_step_rk.h"
#include <stdlib.h>
#include <stdio.h>

int main()
{
    void get_f_sho(double t_n, double *u_n, double *f)
    {
        f[0] = u_n[1];
        f[1] = -u_n[0];
    }
    
    double u0[] = {0, 1};
    int maxsteps = 1000;
    double tstart = 0;
    double tend = 50;
    
    double *u = malloc(2 * maxsteps * sizeof(double));
    double h = (tend - tstart) / maxsteps;
    euler(u, u0, 2, maxsteps, h, tstart, get_f_sho);
    for (int n = 0; n < maxsteps; ++n){
        printf("%.4e: (%.4e, %.4e)\n", tstart + h*(n+1), u[2*n], u[2*n+1]);
    }
    
    results_rkab *res = rk45(u0, 2, maxsteps, 1e-3, tstart, tend, get_f_sho);
    printf("%d, %d\n", res->numsteps, res->numfailures);
    for (int n = 0; n < res->numsteps; ++n){
        printf("%.4e: (%.4e, %.4e)\n", res->t[n], res->u[2*n], res->u[2*n+1]);
    }
}
