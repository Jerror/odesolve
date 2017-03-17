/** @file
 * @brief Runge-Kutta-Fehlberg method: adaptive step Runge-Kutta implementation
 * @details Provides a C interface; see adaptive_step_rk.h for details.
 * @author Jeremiah O'Neil
 * @copyright GNU Public License. */
#include "rkab.hpp" // templates
#include "adaptive_step_rk.h" // interface

/// Runge-Kutta matrix nonzero part, transposed and flattened.
#define A \
    {1/4.L, 3/32.L,  1932/2197.L,  439/216.L,      -8/27.L, \
            9/32.L, -7200/2197.L,         -8,            2, \
                     7296/2197.L, 3680/513.L, -3544/2565.L, \
                                 -845/4104.L,  1859/4104.L, \
                                                  -11/40.L}
/// Weights, with leading zero removed.
#define C \
    {1/4.L, 3/8.L, 12/13.L, 1, 1/2.L}
/// Nodes of low-order method.
#define BA \
    {25/216.L, 0, 1408/2565.L, 2197/4104.L, -1/5.L}
/// Nodes of high-order method.
#define BB \
    {16/135.L, 0, 6656/12825.L, 28561/56430.L, -9/50.L, 2/55.L}

/** Instantiate as defined in adaptive_step_rk.h, binding the modified Butcher
 * tableau to an rkab function instance.
 * Should be used through MAP_TARGETS_TO(). */
#define INST_RK45(T, Tid) \
    static const T *ba##Tid = (T[])BA, *bb##Tid = (T[])BB, \
                   *a##Tid = (T[])A, *c##Tid = (T[])C;     \
    INST_RKAB(45##Tid, T, T, 5, 6,                         \
              ba##Tid, bb##Tid, a##Tid, c##Tid)            \
    INST_RKAB(45_arrtol##Tid, T, T *, 5, 6,                \
              ba##Tid, bb##Tid, a##Tid, c##Tid)

MAP_TARGETS_TO(INST_RK45)
