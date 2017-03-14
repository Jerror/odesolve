SHELL=zsh

INC_DIR = .
CFLAGS = -Wall -Winline -I$(INC_DIR) -O3 -g3

EXEC = valgrind -v 
# time gdb 

all: libode.so

# I'll compile statically for anaconda python ctypes import, lest anaconda gcc4 cause complications.
euler.o: euler.c $(INC_DIR)/euler.h
	gcc -c $(CFLAGS) -std=c99 -Wl,static -fPIC $< -o $@

%.cpp.o: %.cpp $(INC_DIR)/rkab.hpp $(INC_DIR)/adaptive_step_rk.h
	g++ -static-libstdc++ -c $(CFLAGS) -std=c++11 -Wl,static -fPIC $< -o $@

libode.so: rkab_results.cpp.o rk45.cpp.o euler.o
	g++ -static-libstdc++ -std=c++11 -shared -Wl,-soname,libode.so -o libode.so euler.o rk45.cpp.o -lc -lm

test/test: test/main.c libode.so
	- cp libode.so test
	cd test && \
	gcc -o test main.c -lode

.PHONY: run_test

run_test: test/test
	cd test && $(EXEC) ./test
