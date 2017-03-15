/* Templates can't extern "C" so I have to manually instantiate and extern all
   of the objects I want to expose. I'll use macros to help. */
#ifndef INC_ADAPTIVE_STEP_RK_h // #include guard
#define INC_ADAPTIVE_STEP_RK_h // ensure this is included at most once per unit

// Map macro functions of the form F(T, Tid) to target types and type IDs
#define MAP_TARGETS_TO(F) \
    F(double, )           \
    F(float, _f)          \
    F(long double, _ld)
// complex, MPFR, etc...

//// Macros for implementation

// To instantiate delete_rkab_results under suffixed symbol for type T
#define INST_DELETE_RKAB_RESULTS(T, Tid) \
    void delete_results_rkab##Tid(results_rkab<T> *results) \
    {   delete_results_rkab<T>(results);   }
// I'll use this in rkab_results.cpp through MAP_TARGETS_TO().

// To instantiate rkab under suffixed symbol with types and tableau bound
#define INST_RKAB(sfx, T, tolT, astages, bstages, ba, bb, a, c) \
    results_rkab<T> *rk##sfx(T *u_init, int dim, int maxsteps, tolT tol,  \
                             T t, T t_end, void (*get_f)(T, T*, T*))      \
    {   return rkab<T, tolT>(astages, bstages, ba, bb, a, c,              \
                             u_init, dim, maxsteps, tol, t, t_end, get_f);}
// I'll use this in implementation files (eg., 'rk45.cpp', 'rk23.cpp') through
//  MAP_TARGETS_TO(). astages and bstages must be literal integers!

//// Expose C-extern interfaces of instantiated functions

#ifdef __cplusplus // if this header was included in a C++ context
extern "C" { // use C linkage. Forbids symbol mangling (and thus overloading)
#define RESULTS_RKAB(T, Tid) results_rkab<T> // templated for safe compiling
#else // included in C context to provide the C interface to the library
#define TYPEDEF_RESULTS_RKAB(T, Tid) \
    typedef struct results_rkab##Tid {             \
        int numsteps; T *t; T *u; int numfailures; \
    } results_rkab##Tid;
MAP_TARGETS_TO(TYPEDEF_RESULTS_RKAB) // very sorry about this.
#define RESULTS_RKAB(T, Tid) results_rkab##Tid // C'est la C
#endif // Thus this header tells C++ how to compile and C how to use.

// rkab_results.cpp
#define EXPOSE_DELETE_RKAB_RESULTS(T, Tid) \
    void delete_results_rkab##Tid(RESULTS_RKAB(T, Tid)*);

MAP_TARGETS_TO(EXPOSE_DELETE_RKAB_RESULTS)

//  rk45.cpp
#define EXPOSE_RK45(T, Tid) \
    RESULTS_RKAB(T, Tid) *rk45##Tid                                        \
                                (T *u_init, int dim, int maxsteps, T tol,  \
                                 T t, T t_end, void (*get_f)(T, T*, T*));  \
    RESULTS_RKAB(T, Tid) *rk45_arrtol##Tid                                 \
                                (T *u_init, int dim, int maxsteps, T *tol, \
                                 T t, T t_end, void (*get_f)(T, T*, T*));

MAP_TARGETS_TO(EXPOSE_RK45)

#ifdef __cplusplus
} // closing brace for extern "C"
#endif

#endif // #include guard
