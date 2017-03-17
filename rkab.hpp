/** @file
 * @brief Templates for adaptive step size Runge-Kutta solver.
 * @details Provides a template rkab() for adaptive Runge-Kutta functions of
 * arbitrary Butcher tableau and floating-point compatible data type, which
 * accept tolerance as either a number or an array.
 * Provides templates for the return type of rkab() and its API, and
 * templates for auxilliary functions used by rkab().
 * @author Jeremiah O'Neil
 * @copyright GNU Public License. */
#ifndef INC_RKAB_hpp // #include guard
#define INC_RKAB_hpp // ensure this file is included at most once per unit

#include <assert.h> // for assertions
#include <alloca.h> // for allocating memory on the stack
#include <algorithm> // for copy, min, max, fill
#include <boost/math/special_functions/ulp.hpp> // for ulp (MATLAB's eps)
#include <limits> // for numeric_limits
#include <type_traits> // for is_pointer
#include <vector> // for vectors (dynamic size arrays)
#include <cmath> // for abs, pow
using namespace std;


/** @brief Structure template for the return type of adaptive step size 
 * Runge-Kutta methods.
 * @remark The destructor is defined separately; see delete_results_rkab().
 * @tparam T Floating-point compatible data type. */
template<typename T>
struct results_rkab
{
    /// The number of accepted steps
    int numsteps;
    /// Pointer to the array of parameters corresponding to the solution
    T *t;
    /// Pointer to the solution trajectory (array of state arrays)
    T *u;
    /// The number of steps where a failure occured
    int numfailures;
};

/** @brief Function template for deleting results_rkab instances.
 * @details The destructor of results_rkab, separated from the structure
 * for C programs to release the memory via callback. */
template<typename T>
void delete_results_rkab(results_rkab<T> *results)
{
    delete [] results->t;
    delete [] results->u;
    delete results;
}

/** @brief Function template for calculating relative acceptability for scalar
 * tolerance.
 * @details I define acceptability as the minimum over all elements of the
 * ratio of tolerance to error.
 * @param dim The dimension of the system.
 * @param ua A proposed state of the system.
 * @param ub A proposed state of the system computed to a higher order.
 * @param tol The relative tolerance for the local error of the system.
 * @tparam T Floating-point compatible data type. */
template<typename T>
T acceptability_rel(int dim, T *ua, T *ub, T tol)
{
    T acc = numeric_limits<T>::infinity();
    for (int i = 0; i < dim; ++i) {
        acc = min(acc, abs(tol * ub[i] / (ua[i] - ub[i])));
    }
    return acc;
}

/** @brief Function template for calculating relative acceptability for an array
 * of tolerances.
 * @details Overload of the scalar method for "tol" a pointer type, assumedly 
 * to an array of tolerances, one for each component of u. I define 
 * acceptability as the minimum over all elements of the ratio of tolerance
 * to error.
 * @param dim The dimension of the system.
 * @param ua A proposed state of the system.
 * @param ub A proposed state of the system computed to a higher order.
 * @param tol A pointer to an array of relative tolerances for the local
 * error of the system.
 * @tparam T Floating-point compatible data type. */
template<typename T>
T acceptability_rel(int dim, T *ua, T *ub, T *tol)
{
    T acc = numeric_limits<T>::infinity();
    for (int i = 0; i < dim; ++i) {
        acc = min(acc, abs(tol[i] * ub[i] / (ua[i] - ub[i])));
    }
    return acc;
}
    
/** @brief Function template for adaptive step size Runge-Kutta methods.
 * @details Solves a given system over parameterized domain to given relative
 * tolerance in local error. Dynamically allocates memory for the solution and
 * returns a pointer to a results_rkab instance containing that solution.
 * @param astages,bstages,ba,bb,a,c Modified extended Butcher tableau to be
 * bound on instantiation, defining a particular method. In relation to a
 * standard extended Butcher tableau, 'a' is expected to be transposed and
 * flattened with the  zero half removed and 'ba' and 'c' to have the trailing
 * and leading zeroes respectively removed. 'astages' and 'bstages' are the
 * number of stages of the low- and high-order methods respectively.
 * @param u_init The initial state array of the system.
 * @param dim The dimension of the system.
 * @param maxsteps The maximum number of iterations to run.
 * @param tol The relative tolerance or a pointer to an array of relative
 * tolerances for the local error of the system at each step.
 * @param t The initial value of the system parameter.
 * @param t_end The target value of the system parameter.
 * @param get_f A callback function get_f(t, *u_t, *f) which writes the
 * derivative of u at system parameter t and state u_t to array f.
 * @tparam T Floating-point compatible data type.
 * @tparam Ttol Tolerance type: should be (scalar) T or (array) T*. */
template<typename T, typename tolT>
struct results_rkab<T> *rkab(const int astages, const int bstages,
                             const T *ba, const T *bb, const T *a, const T *c,
                             T *u_init, int dim, int maxsteps, tolT tol,
                             T t, T t_end, void (*get_f)(T, T*, T*))
{
    assert(astages < bstages);
    // Set some resonable acceptance patterns
    const T acc_scale = pow(0.9, bstages);
    const T max_adapt = 10;
    const T min_adapt = 0.5;

    // Initialize dynamically sized memory for holding results
    vector<T> tvec;
    vector<T> u;
    // Allocate temporary arrays. On heap for safety; don't forget to delete!
    T *ua = new T[dim];
    T *ub = new T[dim];
    T *u_k = new T[dim];
    T *u_prev = new T[dim]; // for easy re-init on failure
    T *f = new T[bstages * dim];
    T *haf = new T[dim * (bstages - 1)];

    // Initialize
    int numfailures = 0;
    int numsteps = 0;
    copy(u_init, u_init + dim, ua);
    copy(u_init, u_init + dim, ub);
    copy(u_init, u_init + dim, u_k);
    copy(u_init, u_init + dim, u_prev);

    // 'a' is a flattened sparse rep of a triangular matrix; I'll iterate
    //  through its elements by advancing this pointer.
    const T *a_k = a;

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
        if (t_dir * (t_end - t - h) < 0){
            h = t_end - t;
        }
        
        // Loop for advancing one step
        while (true)
        {
            // Please study these lines carefully. How do they work?    
            // (My tableau has leading zeroes dropped with 'a' transposed)
            fill(haf, haf + (bstages - 1) * dim, 0); // reset temp 'haf' array
            a_k = a; // reset to row k = 0
            get_f(t, u_k, f); // stage 1
            for (int k = 0; k < bstages - 1; ++k) // stages 2 through last
            { 
                for (int i = 0; i < dim; ++i){
                    T hf_ki = h * f[k * dim + i];
                    // Update ua[i] and ub[i] now since hf_ki computed
                    ub[i] += bb[k] * hf_ki;
                    if (k < astages){ // Likely always true and optimized out
                        ua[i] += ba[k] * hf_ki;
                    }
                    for (int j = k; j < bstages - 1; ++j){
                        haf[i*(bstages - 1) + j] += hf_ki * a_k[j - k];
                    } 
                    u_k[i] = u_prev[i] + haf[i*(bstages - 1) + k];    
                }
                get_f(t + h * c[k], u_k, &f[(k + 1) * dim]);            
                a_k += bstages - 1 - k; // advance to next row
            }
            for (int i = 0; i < dim; ++i){ // finish last stage
                ub[i] += h * bb[bstages - 1] * f[(bstages - 1) * dim + i]; 
            } // Note: I asserted astages < bstages, so ua is already complete.

            // acceptability = (tolerance) / (relative error)
            T acceptability = acc_scale * acceptability_rel(dim, ua, ub, tol);
            // If the step is acceptable or the step size is minimal
            if (acceptability > 1 || abs(h) <= hmin)
            { // Accept the step
                ++numsteps;
                t += h;
                tvec.push_back(t);
                copy(ub, ub + dim, ua);
                copy(ub, ub + dim, u_prev);
                copy(ub, ub + dim, u_k);
                u.insert(u.end(), ub, ub + dim);
                // Adapt step size; don't increase by a factor > max_adapt
                h *= min(max_adapt, (T)(pow(acceptability, 1.0/bstages)));
                break;
            }
            else
            { // Reject the step
                if (!failures)
                {
                    failures = true;
                    ++numfailures;
                    // Adapt step size; don't decrease by a factor < min_adapt
                    h *= max(min_adapt, (T)(pow(acceptability, 1.0/bstages)));
                } else { // We underestimated error! Be pessimistic.
                    h *= min_adapt;
                }
                copy(u_prev, u_prev + dim, ua);
                copy(u_prev, u_prev + dim, ub);
                copy(u_prev, u_prev + dim, u_k);
            }
        }
    }   
    assert(tvec.size() == (size_t)numsteps);
    assert(u.size() == (size_t)numsteps * dim); 

    // Allocate memory dynamically (persists outside function) for results
    T *tarr = new T[numsteps]; // time array
    T *uarr = new T[numsteps * dim]; // u array
    results_rkab<T> *results = new results_rkab<T> { 
        numsteps,
        tarr, // pointer
        uarr, // pointer
        numfailures
    };
    // Copy time and u data from the vectors.
    //  Can't return vectors directly because I want a C-compatible interface.
    copy(tvec.begin(), tvec.end(), tarr);
    copy(u.begin(), u.end(), uarr);
    
    // delete temporary arrays
    delete [] ua;
    delete [] ub;
    delete [] u_k;
    delete [] u_prev;
    delete [] f;
    delete [] haf;
    // vectors are RAII; they'll be deleted automatically.

    return results;
}

#endif // #include guard
