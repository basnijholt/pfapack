OBJECTS = error.o matrix_tools.o matrix_tools_s.o matrix_tools_d.o	\
	  matrix_tools_c.o matrix_tools_z.o check_decomp.o		\
	  check_decomp_s.o check_decomp_d.o check_decomp_c.o		\
	  check_decomp_z.o check_pfaffian.o check_pfaffian_s.o		\
	  check_pfaffian_d.o check_pfaffian_c.o check_pfaffian_z.o	\
	  test_ilaenv.o

FORT = gfortran
FFLAGS = -O3 -fimplicit-none -I..
LFLAGS = ../libpfapack.a -llapack -lblas

test: test_pfaff.out test_skew.out
	./test_skew.out
	./test_pfaff.out

test_pfaff.out: test_pfaff.f90 $(OBJECTS) $(MODULES) ../libpfapack.a
	$(FORT) $(FFLAGS) test_pfaff.f90 -o test_pfaff.out $(OBJECTS) $(LFLAGS)

test_skew.out: test_skew.f90 $(OBJECTS) $(MODULES) ../libpfapack.a
	$(FORT) $(FFLAGS) test_skew.f90 -o test_skew.out $(OBJECTS) $(LFLAGS)

clean:
	rm -f *.out
	rm -f *.o
	rm -f *.mod

.SUFFIXES: .o .f .f90

.f90.o:
	$(FORT) $(FFLAGS) -c $< -o $@

.f.o:
	$(FORT) $(FFLAGS) -c $< -o $@
