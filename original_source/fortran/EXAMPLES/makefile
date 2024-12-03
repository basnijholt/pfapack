PROGS = example1.out example1_ws.out example2.out example3.out		\
example4.out example4_upper.out example5.out example6.out example7.out

FORT = gfortran
FFLAGS = -O3 -I.. -fimplicit-none
LFLAGS = ../libpfapack.a -llapack -lblas print_mat.o

all: print_mat.mod $(PROGS) ../libpfapack.a

.SUFFIXES: .f .f90 .mod .out

.f.out:
	$(FORT) $(FFLAGS) $< -o $@ $(LFLAGS)
.f90.out:
	$(FORT) $(FFLAGS) $< -o $@ $(LFLAGS)
.f90.mod:
	$(FORT) $(FFLAGS) $< -c

clean:
	rm -f *.out
	rm -f *.o
	rm -f *.mod