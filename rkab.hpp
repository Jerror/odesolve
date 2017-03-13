// Licensed GPLv2
#ifndef INC_RKAB_hpp // #include guard
#define INC_RKAB_hpp // ensure this file is only included once

#include "ode.h"
#include <assert>
#include <cmath>


template<typename T>
struct results_rkab
{
    int numsteps;
    T *t;
    T *u;
    int numfailures;
};

template<typename T>
void delete_results_rkab(results_rkab<T> *results)
{
    delete results->t;
    delete results->u;
    delete results;
}

template<typename T, typename tolT>
T acceptability_rel(int dim, T *ya, T *yb, tolT tol)
{
    T acc = std::numeric_limits<T>::infinity();
    T acc_temp;
    for (int i = 0; i < dim; ++i)
    {
        if (std::is_pointer<tolT>::value){ // vector of tolerances
            acc_temp = abs(tol[i] * yb[i] / (ya[i] - yb[i]));
        } else { // scalar tolerance
            acc_temp = abs(tol * yb[i] / (ya[i] - yb[i]));
        }
        acc = std::min(acc, acc_temp);
    }
    return acc;
}

template<typename T, typename tolT>
struct results_rkab<T> *rkab(int astages, int bstages,
                             T *ba, T *bb, T *a, T *c,
                             T *u_init, int dim, int maxsteps, tolT tol,
                             T t, T t_end, derivative_function get_f)
{
    assert(astages < bstages);

    // Initialize dynamically sized memory for holding results
    std::vector<T> tvec;
    std::vector<T> u;
    // Allocated temporary arrays on the stack
    T *ua = (T *)alloca(dim * sizeof T);
    T *ub = (T *)alloca(dim * sizeof T);
    T *u_s = (T *)alloca(dim * sizeof T);
    T *k = (T *)alloca(bstages * dim * sizeof T);

    // Initialize
    int numfailures = 0;
    int numsteps = 0;
    std:copy(u_init, u_init + dim, ua);
    std:copy(u_init, u_init + dim, ub);
    std:copy(u_init, u_init + dim, u_s);

    // I'll advance and dereference these aliases, rather than index arrays
    T *kk = k;
    T *aa = a;
    // This pointer will make it easier to reset temporary arrays on failure
    T *u_prev = u_init;

    // Guess an initial step size
    T h = std::min(abs(t_end - t)/10, 0.1);
    // Choose the sign of h according to the direction of propagation
    int t_dir = (t_end >= t) ? 1 : -1;
    h *= t_dir

    // Main loop
    while (numsteps < maxsteps && t_dir * (t_end - t) > 0)
    {
        bool failures = false;
        T minh = 16 * std::numeric_limits<T>::epsilon(t);
        // minh is the minimum meaningful magnitude of h
        if (abs(h) < hmin) { // abs(h) should be at least as large as minh
            h = t_dir * hmin;
        }
        // But make sure to hit the last step exactly
        if (h - t_end - t < 0){
            h = t_end - t;
        }
        
        // Loop for advancing one step
        while (true)
        {
            // Please study these lines carefully. How do they work?    
            get_f(t, u_s, kk); // stage 1
            for (int s = 1; s < bstages; ++s){ // stages 2 through last
                for (int j = 0; j < s; ++j){
                    for (int i = 0; i < dim; ++i)
                    {
                        u_s[i] += h * &aa * &kk;
                        // Update ua[i] and ub[i] now since &kk is on register
                        ub[j] += h * bb[s - 1] * &kk;
                        if (s <= astages){
                            ua[j] += h * ba[s - 1] * &kk;
                        }
                        ++kk;
                    }
                    get_f(t + h * c[s], u_s, kk);            
                    ++aa;
                }
                std::copy(u_prev, u_prev + dim, u_s);
            } // Here s == bstages, kk == &k[(bstages - 1) * dim]
            for (int i = 0; i < dim; ++i){
                ub[j] += h * bb[s - 1] * kk[i]; // finish last stage
            } // Note: I asserted astages < bstages, so ua is already complete.
            kk = k; aa = a;

            // acceptability = (tolerance) / (relative error)
            T acceptability = acceptability_rel(dim, ua, ub, tol);            
            // If the step is acceptable or the step size is minimal
            if (acceptability > 1 || h <= hmin)
            { // Accept the step
                ++numsteps;
                t += h;
                tvec.push_back(t)
                std:copy(ub, ub + dim, ua);
                std:copy(ub, ub + dim, u_s);
                // This should be the only place I append to u
                u.assign(ub, ub + dim);
                // because this pointer dies when the vector needs more memory
                u_prev = &u.back() - dim;
                // Adapt the step size; don't increase by a factor > 10
                h *= std::min(10, 0.8 * std::pow(acceptability, 1.0/bstages));
                break;
            }
            else
            { // Reject the step
                if (!failures)
                {
                    failures = true;
                    ++numfailures;
                    // Adapt the step size; don't decrease by a factor < 1/2
                    h *= std::max(0.5,
                                  0.8 * std::pow(acceptability, 1.0/bstages));
                } else { // We underestimated error! Be pessimistic.
                    h *= 0.5;
                }
                std:copy(u_prev, u_prev + dim, ua);
                std:copy(u_prev, u_prev + dim, ub);
                std:copy(u_prev, u_prev + dim, u_s);
            }
        }
    }   
    assert(tvec.size() == numsteps);
    assert(u.size() == numsteps * dim); 

    // Allocate memory dynamically (persists outside function) for results
    tarr = new T[numsteps], // time array
    uarr = new T[numsteps * dim], // u array
    results_rkab *results = new results_rkab { 
        numfailures,
        tarr, // pointer
        uarr, // pointer
        numsteps
    }
    // Copy time and u data from the vectors.
    //  Can't return vectors directly because I want a C-compatible interface.
    std::copy(tvec.begin(), tvec.end(), tarr);
    std::copy(u.begin(), u.end(), uarr);
    return results;
}


/* Templates can't extern "C" so I'll manually instantiate and extern all of
   the objects I want to expose. I'll use macros to help. */

// C-extern results_rkab and destructor under suffixed symbol for type T
#define EXTERNC_RKAB_RESULTS(T) \
    extern "C" typedef results_rkab<T> results_rkab_##T;              \
    extern "C" void delete_results_rkab_##T(results_rkab_##T results) \
    {   delete_results_rkab<T>(results)   }

EXTERNC_RKAB_RESULTS(float)
EXTERNC_RKAB_RESULTS(double)
typedef long double ldouble; // One-symbol alias for 80bit float (usually)
EXTERNC_RKAB_RESULTS(ldouble)
// complex, MPFR, &c...

// C-extern rkab instance under symbol "name" with types and tableau bound 
#define EXTERNC_RKAB(name, T, tolT, astages, bstages, ba, bb, a, c) \
    extern "C" results_rkab_##T name(T *u_init, int dim, int maxsteps,   \
                                     tolT tol, T t, T t_end,             \
                                     derivative_function get_f)          \
    {   return rkab<T, tolT>(astages, bstages, ba, bb, a, c,             \
                             u_init, dim, maxsteps, tol, t, t_end, get_f)}
// I'll use this in implementation files (eg., 'rk45.cpp', 'rk23.cpp').

#endif // #include guard
