// Bogacki-Shampine method
#include "rkab.hpp" // templates
#include "adaptive_step_rk.h" // interface

// Vs Butcher tableau: trailing zeroes removed
// Additionally, A is transposed and flattened: convenient for C/C++ iteration.
#define A {1/2.L,   0,  2/9.L,\
                 3/4.L, 1/3.L,\
                        4/9.L}
#define C {1/2.L, 3/4.L, 1}
#define BA {2/9.L, 1/3.L, 4/9.L}
#define BB {7/24.L, 1/4.L, 1/3.L, 1/8.L}

// Instantiate as defined in adaptive_step_rk.h

#define INST_RK23(T, Tid) \
    static const T *ba##Tid = (T[])BA, *bb##Tid = (T[])BB, \
                   *a##Tid = (T[])A, *c##Tid = (T[])C;     \
    INST_RKAB(23##Tid, T, T, 3, 4,                         \
              ba##Tid, bb##Tid, a##Tid, c##Tid)            \
    INST_RKAB(23_arrtol##Tid, T, T *, 3, 4,                \
              ba##Tid, bb##Tid, a##Tid, c##Tid)

MAP_TARGETS_TO(INST_RK23)
