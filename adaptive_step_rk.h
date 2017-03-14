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
    results_rkab<T> *rk##astages##bstages##sfx                            \
                        (T *u_init, int dim, int maxsteps, tolT tol,      \
                         T t, T t_end, void (*get_f)(T, T*, T*))          \
    {   return rkab<T, tolT>(astages, bstages, ba, bb, a, c,              \
                             u_init, dim, maxsteps, tol, t, t_end, get_f);}
// I'll use this in implementation files (eg., 'rk45.cpp', 'rk23.cpp') through
//  MAP_TARGETS_TO(). astages and bstages must be literal integers!

//// Expose C-extern interfaces of instantiated functions

// rkab_results.cpp
#define EXPOSE_DELETE_RKAB_RESULTS(T, Tid) \
    extern "C" void delete_results_rkab##Tid(results_rkab<T>*);

MAP_TARGETS_TO(EXPOSE_DELETE_RKAB_RESULTS)

//  rk45.cpp
#define EXPOSE_RK45(T, Tid) \
    extern "C" results_rkab<T> *rk45##Tid                             \
                                   (T *u_init, int dim, int maxsteps, \
                                    T tol, T t, T t_end,              \
                                    void (*get_f)(T, T*, T*));        \
    extern "C" results_rkab<T> *rk45_arrtol##Tid                      \
                                   (T *u_init, int dim, int maxsteps, \
                                    T *tol, T t, T t_end,             \
                                    void (*get_f)(T, T*, T*));

MAP_TARGETS_TO(EXPOSE_RK45)

#endif // #include guard
