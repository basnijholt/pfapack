CC = gcc
CFLAGS = -O3 -I ../
FORT = gfortran
FFLAGS = -O3 -fimplicit-none -I ../../fortran/ -Wall
LFLAGS = -llapack -lblas -lgfortran -lm

all: test

test: test.out
	./test.out

test.out: check.o ../libcpfapack.a ../../fortran/libpfapack.a
	$(CC) $(CFLAGS) test.c -o test.out check.o ../libcpfapack.a \
	    ../../fortran/libpfapack.a $(LFLAGS)

clean:
	rm -f *.o
	rm -f *.out

.SUFFIXES: .o .f90

.f90.o:
	$(FORT) $(FFLAGS) -c $< -o $@
