# libode {#mainpage}
(Dox hosted at https://jerror.github.io/odesolve/)

A mixed C/C++ library of custom ODE solvers. Set up for C linkage and static compilation, as suitable for use in Python through ctypes, but written mostly as templated C++, for safe and simple generalization of the methods.

Presently, the library provides a dead simple fixed-step euler method (source euler.c) and a variety of adaptive step size Runge-Kutta methods for scalar and vector relatve local tolerance (instantiated from a template in rkab.hpp). Adaptive step size Runge-Kutta methods of any order for any floating-point compatible data type can be trivially instantiated given the Butcher tableau; see rk45.cpp for an example. The C interface of the euler method is provided by inclusion of euler.h, and those of the Runge-Kutta methods are provided by inclusion of adaptive\_step\_rk.h.

The memory required for the solution of the euler method is known at runtime, so the euler method takes the preallocated memory as an argument, writes it, and returns nothing; but the memory required for the adaptive Runge-Kutta methods is not known until the method completes, so those methods dynamically allocate the memory required and return a pointer to a class encapsulating the results, namely, results\_rkab. The class is templated to instantiate for any requested floating-point compatible data type in rkab.hpp, and in compilation the file results\_rkab.cpp instantiates the class according to the definitions in adaptive\_step\_rk.h, which then exports the definitions with C linkage under suffixed names. The same is done for the adaptive step size Runge-Kutta methods for different types, via the implementation files rk12.cpp, rk23.cpp and rk45.cpp. Instances of results\_rkab can be freed by calling delete\_results\_rkab, suffixed for the appropriate type.

Presently, the Runge-Kutta methods and results class are exported for data types float, double and long double under symbols suffixed by \_f, \_d and \_g respectively; no suffix aliases \_d (double). The euler method is only provided for double type data. Runge-Kutta methods accepting array tolerance (as opposed to scalar) are exported with the suffix \_arrtol in addition to (preceding) the suffix denoting the data type.

The symbols currently exported are:
- euler
- rk12
- rk23
- rk45
- rk12\_f
- rk23\_f
- rk45\_f
- rk12\_d
- rk23\_d
- rk45\_d
- rk12\_f
- rk23\_f
- rk45\_f
- rk12\_arrtol
- rk23\_arrtol
- rk45\_arrtol
- rk12\_arrtol\_f
- rk23\_arrtol\_f
- rk45\_arrtol\_f
- rk12\_arrtol\_d
- rk23\_arrtol\_d
- rk45\_arrtol\_d
- rk12\_arrtol\_f
- rk23\_arrtol\_f
- rk45\_arrtol\_f
- results\_rkab
- results\_rkab_f
- results\_rkab_d
- results\_rkab_g
- delete_results\_rkab
- delete_results\_rkab_f
- delete_results\_rkab_d
- delete_results\_rkab_g

