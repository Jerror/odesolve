SHELL=zsh

PAROPTS=-fopenmp
OMP_NUM_THREADS=12
OMP_SCHEDULE="auto"

INC_DIR = .
CFLAGS = -Wall -Winline -I$(INC_DIR) -O3

all: libode.so

# I need to link statically for anaconda python ctypes import, else anaconda gcc4 will fruitlessly search for old GOMP headers.
libode.so: $(INC_DIR)/ode.h rk45.c euler.c
	gcc -c $(CFLAGS) -Wl,static -DOMP_NUM_THREADS=$(OMP_NUM_THREADS) $(PAROPTS) -fPIC -o euler.o  euler.c && \
	gcc -shared -Wl,-soname,libode.so -o libode.so euler.o rk45.o -lc
