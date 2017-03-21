/** @file 
 * @brief Heun-Euler method: adaptive step Runge-Kutta implementation
 * @details Provides a C interface; see adaptive_step_rk.h for details.
 * @author Jeremiah O'Neil
 * @copyright GNU Public License. */
#include "rkab.hpp" // templates
#include "adaptive_step_rk.h" // interface

/// Runge-Kutta matrix nonzero part, transposed and flattened.
#define A \
    {1}
/// Weights, with leading zero removed.
#define C \
    {1}
/// Nodes of low-order method.
#define BA \
    {1}
/// Nodes of high-order method.
#define BB \
    {1/2.L, 1/2.L}

/** Instantiate as defined in adaptive_step_rk.h, binding the modified Butcher
 * tableau to an rkab function instance.
 * Should be used through MAP_TARGETS_TO(). */
#define INST_RK12(T, Tid) \
    static const T *ba##Tid = (T[])BA, *bb##Tid = (T[])BB, \
                   *a##Tid = (T[])A, *c##Tid = (T[])C;     \
    INST_RKAB(12##Tid, T, T, 2, 1, 2,                      \
              ba##Tid, bb##Tid, a##Tid, c##Tid)            \
    INST_RKAB(12_arrtol##Tid, T, T *, 2, 1, 2,             \
              ba##Tid, bb##Tid, a##Tid, c##Tid)

MAP_TARGETS_TO(INST_RK12)
