#include "../adaptive_step_rk.h"
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
    results_rkab *res = rk45(u0, 2, 1e-3, maxsteps, tstart, tend, get_f_sho);
}
