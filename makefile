SHELL=zsh

INC_DIR = .
CFLAGS = -Wall -Winline -I$(INC_DIR) -O3

all: libode.so

# I'll compile statically for anaconda python ctypes import, lest anaconda gcc4 cause complications.
%.o: %.c $(INC_DIR)/ode.h
	gcc -c $(CFLAGS) -std=c99 -Wl,static -fPIC $< -o $@

%.cpp.o: %.cpp $(INC_DIR)/rkab.hpp
	g++ -c $(CFLAGS) -std=c++11 -Wl,static -fPIC $< -o $@

libode.so: $(INC_DIR)/ode.h rk45.cpp.o euler.o
	g++ -std=c++11 -shared -Wl,-soname,libode.so -o libode.so euler.o rk45.cpp.o -lc -lm
