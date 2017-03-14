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
    T acc = std::numeric_limits<T>::infinity();
    for (int i = 0; i < dim; ++i) {
        acc = std::min(acc, std::abs(tol * yb[i] / (ya[i] - yb[i])));
    }
    return acc;
}
// Overload for an array of tolerances for each component of u
template<typename T>
T acceptability_rel(int dim, T *ya, T *yb, T *tol)
{
    T acc = std::numeric_limits<T>::infinity();
    for (int i = 0; i < dim; ++i) {
        acc = std::min(acc, std::abs(tol[i] * yb[i] / (ya[i] - yb[i])));
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
    std::vector<T> tvec;
    std::vector<T> u;
    // Allocate temporary arrays on the stack
    T *ua = (T *)alloca(dim * sizeof(T));
    T *ub = (T *)alloca(dim * sizeof(T));
    T *u_s = (T *)alloca(dim * sizeof(T));
    T *k = (T *)alloca(bstages * dim * sizeof(T));

    // Initialize
    int numfailures = 0;
    int numsteps = 0;
    std::copy(u_init, u_init + dim, ua);
    std::copy(u_init, u_init + dim, ub);
    std::copy(u_init, u_init + dim, u_s);

    // I'll advance and dereference these aliases, rather than index arrays
    T *kk = k;
    const T *aa = a;
    // This pointer will make it easier to reset temporary arrays on failure
    T *u_prev = u_init;

    // Guess an initial step size
    T h = std::min(std::abs(t_end - t)/10, (T)0.1);
    // Choose the sign of h according to the direction of propagation
    int t_dir = (t_end >= t) ? 1 : -1;
    h *= t_dir;

    // Main loop
    while (numsteps < maxsteps && t_dir * (t_end - t) > 0)
    {
        bool failures = false;
        T hmin = 16 * boost::math::ulp(t);
        // hmin is the minimum meaningful magnitude of h
        if (std::abs(h) < hmin) { // std::abs(h) should be >= hmin
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
                std::copy(u_prev, u_prev + dim, u_s);
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
                std::copy(ub, ub + dim, ua);
                std::copy(ub, ub + dim, u_s);
                // This should be the only place I append to u
                u.assign(ub, ub + dim);
                // because this pointer dies when the vector needs more memory
                u_prev = &u.back() - dim;
                // Adapt the step size; don't increase by a factor > 10
                h *= std::min((T)10,
                              (T)(0.8 * std::pow(acceptability, 1.0/bstages)));
                break;
            }
            else
            { // Reject the step
                if (!failures)
                {
                    failures = true;
                    ++numfailures;
                    // Adapt the step size; don't decrease by a factor < 1/2
                    h *= std::max((T)0.5,
                             (T)(0.8 * std::pow(acceptability, 1.0/bstages)));
                } else { // We underestimated error! Be pessimistic.
                    h *= 0.5;
                }
                std::copy(u_prev, u_prev + dim, ua);
                std::copy(u_prev, u_prev + dim, ub);
                std::copy(u_prev, u_prev + dim, u_s);
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
    std::copy(tvec.begin(), tvec.end(), tarr);
    std::copy(u.begin(), u.end(), uarr);
    return results;
}

#endif // #include guard
