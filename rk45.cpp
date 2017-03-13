#include "rkab.hpp"

#define C {1/4L, 3/8L, 12/13L, 1, 1/2L}
#define A {1/4L,\
           3/32L, 9/32L,\
           1932/2197L, -7200/2197L, 7296/2197L,\
           439/216L, -8, 3680/513L, -845/4104L,\
           -8/27L, 2, -3544/2565L, 1859/4104L, -11/40L}
#define BA {25/216L, 0, 1408/2565L, 2197/4104L, -1/5L}
#define BB {16/135L, 0, 6656/12825L, 28561/56430L, -9/50L, 2/55L}

#define EXTERNC_RK45(T, Tid) \
    const T *ba##Tid = (T[])BA, *bb##Tid = (T[])BB,  \
            *a##Tid = (T[])A, *c##Tid = (T[])C;      \
    EXTERNC_RKAB(rk45##Tid, T, T, 4, 5,              \
                 ba##Tid, bb##Tid, a##Tid, c##Tid)   \
    EXTERNC_RKAB(rk45_arrtol##Tid, T, T *, 4, 5,     \
                 ba##Tid, bb##Tid, a##Tid, c##Tid)

EXTERNC_RK45(double, )
//EXTERNC_RK45(float, _f)
//EXTERNC_RK45(long double, _ld)
