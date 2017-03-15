#ifndef INC_RKAB_hpp // #include guard
#define INC_RKAB_hpp // ensure this file is included at most once per unit

#include <assert.h> // for assertions
#include <alloca.h> // for allocating memory on the stack
#include <algorithm> // for copy, min, max
#include <boost/math/special_functions/ulp.hpp> // for ulp (MATLAB's eps)
#include <limits> // for numeric_limits
#include <type_traits> // for is_pointer
#include <vector> // for vectors (dynamic size arrays)
#include <cmath> // for abs, pow
#include <iostream> // for debugging
using namespace std;


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

template<typename T>
T acceptability_rel(int dim, T *ya, T *yb, T tol) // scalar tolerance
{
    T acc = numeric_limits<T>::infinity();
    for (int i = 0; i < dim; ++i) {
        acc = min(acc, abs(tol * yb[i] / (ya[i] - yb[i])));
    }
    return acc;
}
// Overload for an array of tolerances for each component of u
template<typename T>
T acceptability_rel(int dim, T *ya, T *yb, T *tol)
{
    T acc = numeric_limits<T>::infinity();
    for (int i = 0; i < dim; ++i) {
        acc = min(acc, abs(tol[i] * yb[i] / (ya[i] - yb[i])));
    }
    return acc;
}
    
template<typename T, typename tolT>
struct results_rkab<T> *rkab(int astages, int bstages,
                             const T *ba, const T *bb, const T *a, const T *c,
                             T *u_init, int dim, int maxsteps, tolT tol,
                             T t, T t_end, void (*get_f)(T, T*, T*))
{
    assert(astages < bstages);

    // Initialize dynamically sized memory for holding results
    vector<T> tvec;
    vector<T> u;
    // Allocate temporary arrays on the stack
    T *ua = (T *)alloca(dim * sizeof(T));
    T *ub = (T *)alloca(dim * sizeof(T));
    T *u_s = (T *)alloca(dim * sizeof(T));
    T *u_prev = (T *)alloca(dim * sizeof(T)); // for easy re-init on failure
    T *k = (T *)alloca(bstages * dim * sizeof(T));

    // Initialize
    int numfailures = 0;
    int numsteps = 0;
    copy(u_init, u_init + dim, ua);
    copy(u_init, u_init + dim, ub);
    copy(u_init, u_init + dim, u_s);
    copy(u_init, u_init + dim, u_prev);

    // I'll advance and dereference these aliases, rather than index arrays
    T *kk = k;
    const T *aa = a;

    // Guess an initial step size
    T h = min(abs(t_end - t)/10, (T)0.1);
    // Choose the sign of h according to the direction of propagation
    int t_dir = (t_end >= t) ? 1 : -1;
    h *= t_dir;

    // Main loop
    while (numsteps < maxsteps && t_dir * (t_end - t) > 0)
    {
        bool failures = false;
        T hmin = 16 * boost::math::ulp(t);
        // hmin is the minimum meaningful magnitude of h
        if (abs(h) < hmin) { // abs(h) should be >= hmin
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
                        u_s[i] += h * *aa * *kk;
                        // Update ua[i] and ub[i] now since *kk is on register
                        ub[i] += h * bb[s - 1] * *kk;
                        if (s <= astages){
                            ua[i] += h * ba[s - 1] * *kk;
                        }
                        ++kk;
                    }
                    get_f(t + h * c[s], u_s, kk);            
                    ++aa;
                }
                copy(u_prev, u_prev + dim, u_s);
            }
            for (int i = 0; i < dim; ++i){ // finish last stage
                ub[i] += h * bb[bstages - 1] * kk[i]; 
            } // Note: I asserted astages < bstages, so ua is already complete.
            kk = k; aa = a;

            // acceptability = (tolerance) / (relative error)
            T acceptability = acceptability_rel(dim, ua, ub, tol);            
            // If the step is acceptable or the step size is minimal
            if (acceptability > 1 || h <= hmin)
            { // Accept the step
                ++numsteps;
                t += h;
                tvec.push_back(t);
                copy(ub, ub + dim, ua);
                copy(ub, ub + dim, u_s);
                copy(ub, ub + dim, u_prev);
                u.insert(u.end(), ub, ub + dim);
                // Adapt the step size; don't increase by a factor > 10
                h *= min((T)10, (T)(0.8 * pow(acceptability, 1.0/bstages)));
                break;
            }
            else
            { // Reject the step
                if (!failures)
                {
                    failures = true;
                    ++numfailures;
                    // Adapt the step size; don't decrease by a factor < 1/2
                    h *= max((T)0.5, (T)(0.8*pow(acceptability, 1.0/bstages)));
                } else { // We underestimated error! Be pessimistic.
                    h *= 0.5;
                }
                copy(u_prev, u_prev + dim, ua);
                copy(u_prev, u_prev + dim, ub);
                copy(u_prev, u_prev + dim, u_s);
            }
        }
    }   
    assert(tvec.size() == (size_t)numsteps);
    assert(u.size() == (size_t)numsteps * dim); 

    // Allocate memory dynamically (persists outside function) for results
    T *tarr = new T[numsteps]; // time array
    T *uarr = new T[numsteps * dim]; // u array
    results_rkab<T> *results = new results_rkab<T> { 
        numfailures,
        tarr, // pointer
        uarr, // pointer
        numsteps
    };
    // Copy time and u data from the vectors.
    //  Can't return vectors directly because I want a C-compatible interface.
    copy(tvec.begin(), tvec.end(), tarr);
    copy(u.begin(), u.end(), uarr);
    return results;
}

#endif // #include guard
