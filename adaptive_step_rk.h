/** @file
 * @brief Adaptive step size Runge-Kutta interface
 * @details Templates can't extern "C" so I have to manually instantiate and 
 * extern all of the objects I want to expose. I'll use macros to help. This 
 * header is meant to be included at library compile-time by C++ source to 
 * produce the C-compatible interface, and included at program compile-time by 
 * C to define that interface.
 * @author Jeremiah O'Neil
 * @copyright GNU Public License. */
#ifndef INC_ADAPTIVE_STEP_RK_h // #include guard
#define INC_ADAPTIVE_STEP_RK_h // ensure this is included at most once per unit

/// Map macro functions of the form F(T, Tid) to target types and type IDs
#define MAP_TARGETS_TO(F) \
    F(float, _f)          \
    F(double, _d)         \
    F(long double, _g)    \
    F(double, )
//  complex, MPFR, etc...
// The naming convention matches the Numpy C-compatible dtype character codes.
// The last target simply provides a suffix-free alias for double-precision
//  types and methods (double is the "default").
/* For posterity: a C-incompatible type could be used through a C++ API via:
     ...
      F(C_incompatible_type, _suffix)
  #ifndef __cplusplus
  typedef struct C_incompatible_type C_incompatible_type;
  #endif */


// Macros for implementation

/** @brief Instantiate delete_results_rkab under suffixed symbol for type T
 * @details I'll use this in results_rkab.cpp through MAP_TARGETS_TO().*/
#define INST_DELETE_RESULTS_RKAB(T, Tid) \
    void delete_results_rkab##Tid(results_rkab<T> *results) \
    {   delete_results_rkab<T>(results);   }

/** @brief Instantiate rkab under suffixed symbol with types and tableau bound
 * @details I'll use this in implementation files (eg., 'rk45.cpp', 'rk23.cpp')
 * through MAP_TARGETS_TO(). */
#define INST_RKAB(sfx, T, tolT, astages, bstages, ba, bb, a, c) \
    results_rkab<T> *rk##sfx(T *u_init, int dim, int maxsteps, tolT tol,  \
                             T t, T t_end, void (*get_f)(T, T*, T*))      \
    {   return rkab<T, tolT>(astages, bstages, ba, bb, a, c,              \
                             u_init, dim, maxsteps, tol, t, t_end, get_f);}

// Expose C-extern interfaces of instantiated functions
// @cond EXPOSE

#ifndef __cplusplus 
    // We're included in a C context to define the library interface.
    #define TYPEDEF_RESULTS_RKAB(T, Tid) \
        typedef struct results_rkab##Tid {             \
            int numsteps; T *t; T *u; int numfailures; \
        } results_rkab##Tid;
    MAP_TARGETS_TO(TYPEDEF_RESULTS_RKAB)
    // Very sorry about this. C'est la C.
    #define RESULTS_RKAB(T, Tid) results_rkab##Tid
#else
// We're included in a C++ context for library compilation.
#define RESULTS_RKAB(T, Tid) results_rkab<T> // use the structure template
extern "C" { // use C linkage. Forbids symbol mangling (and thus overloading)
#endif

#define EXPOSE_DELETE_RESULTS_RKAB(T, Tid) \
    void delete_results_rkab##Tid(RESULTS_RKAB(T, Tid)*);
// for results_rkab.cpp
MAP_TARGETS_TO(EXPOSE_DELETE_RESULTS_RKAB)

#define EXPOSE_RKAB(AB, T, Tid) \
    RESULTS_RKAB(T, Tid) *rk##AB##Tid                                      \
                                (T *u_init, int dim, int maxsteps, T tol,  \
                                 T t, T t_end, void (*get_f)(T, T*, T*));  \
    RESULTS_RKAB(T, Tid) *rk##AB##_arrtol##Tid                             \
                                (T *u_init, int dim, int maxsteps, T *tol, \
                                 T t, T t_end, void (*get_f)(T, T*, T*));
//  for rk45.cpp
#define EXPOSE_RK45(T, Tid) EXPOSE_RKAB(45, T, Tid)
MAP_TARGETS_TO(EXPOSE_RK45)
//  for rk23.cpp
#define EXPOSE_RK23(T, Tid) EXPOSE_RKAB(23, T, Tid)
MAP_TARGETS_TO(EXPOSE_RK23)
//  for rk45.cpp
#define EXPOSE_RK12(T, Tid) EXPOSE_RKAB(12, T, Tid)
MAP_TARGETS_TO(EXPOSE_RK12)

#ifdef __cplusplus
} // closing brace for extern "C"
#endif

// @endcond
#endif // #include guard
