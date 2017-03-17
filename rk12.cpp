/// Heun-Euler method: adaptive step Runge-Kutta implementation
/// I'm providing a C interface; see adaptive_step_rk.h for details.
#include "rkab.hpp" // templates
#include "adaptive_step_rk.h" // interface

// Modified Butcher tableau: trailing zeroes removed
// Additionally, A is transposed and flattened: convenient for C/C++ iteration.
#define A {1}
#define C {1}
#define BA {1}
#define BB {1/2.L, 1/2.L}

// Instantiate as defined in adaptive_step_rk.h

#define INST_RK12(T, Tid) \
    static const T *ba##Tid = (T[])BA, *bb##Tid = (T[])BB, \
                   *a##Tid = (T[])A, *c##Tid = (T[])C;     \
    INST_RKAB(12##Tid, T, T, 1, 2,                         \
              ba##Tid, bb##Tid, a##Tid, c##Tid)            \
    INST_RKAB(12_arrtol##Tid, T, T *, 1, 2,                \
              ba##Tid, bb##Tid, a##Tid, c##Tid)

MAP_TARGETS_TO(INST_RK12)
