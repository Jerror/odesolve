/** @file
 * @brief Bogacki-Shampine method: adaptive step Runge-Kutta implementation
 * @details Provides a C interface; see adaptive_step_rk.h for details.
 * @author Jeremiah O'Neil
 * @copyright GNU Public License. */
#include "rkab.hpp" // templates
#include "adaptive_step_rk.h" // interface

/// Runge-Kutta matrix nonzero part, transposed and flattened.
#define A \
    {1/2.L,     0, 2/9.L,\
            3/4.L, 1/3.L,\
                   4/9.L}
/// Weights, with leading zero removed.
#define C \
    {1/2.L, 3/4.L, 1}
/// Nodes of low-order method.
#define BA \
    {2/9.L, 1/3.L, 4/9.L}
/// Nodes of high-order method.
#define BB \
    {7/24.L, 1/4.L, 1/3.L, 1/8.L}

/** Instantiate as defined in adaptive_step_rk.h, binding the modified Butcher
 * tableau to an rkab function instance.
 * Should be used through MAP_TARGETS_TO(). */
#define INST_RK23(T, Tid) \
    static const T *ba##Tid = (T[])BA, *bb##Tid = (T[])BB, \
                   *a##Tid = (T[])A, *c##Tid = (T[])C;     \
    INST_RKAB(23##Tid, T, T, 3, 3, 4,                      \
              ba##Tid, bb##Tid, a##Tid, c##Tid)            \
    INST_RKAB(23_arrtol##Tid, T, T *, 3, 3, 4,             \
              ba##Tid, bb##Tid, a##Tid, c##Tid)

MAP_TARGETS_TO(INST_RK23)
