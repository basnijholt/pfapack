CC = g++
ANSIC = gcc -ansi
C99 = gcc -std=c99
CFLAGS = -O3 -I ../
LFLAGS = ../../fortran/libpfapack.a -llapack -lblas -lgfortran

all: example1_f77.out example3_f77.out example1_c.out \
	example1_c99.out example1_c++.out example2_c.out \
	example3_c.out example4_c.out example5_c.out \
	example6_c.out

example1_f77.out: example1_f77.cc
	$(CC) $(CFLAGS) example1_f77.cc -o example1_f77.out $(LFLAGS)

example3_f77.out: example3_f77.cc
	$(CC) $(CFLAGS) example3_f77.cc -o example3_f77.out $(LFLAGS)

example1_c.out: example1_c.c
	$(ANSIC) $(CFLAGS) example1_c.c -o example1_c.out ../libcpfapack.a $(LFLAGS)

example1_c99.out: example1_c99.c
	$(C99) $(CFLAGS) example1_c99.c -o example1_c99.out ../libcpfapack.a $(LFLAGS)

example1_c++.out: example1_c++.cc
	$(CC) $(CFLAGS) example1_c++.cc -o example1_c++.out ../libcpfapack.a $(LFLAGS)

example2_c.out: example2_c.c
	$(ANSIC) $(CFLAGS) example2_c.c -o example2_c.out ../libcpfapack.a $(LFLAGS)

example3_c.out: example3_c.c
	$(ANSIC) $(CFLAGS) example3_c.c -o example3_c.out ../libcpfapack.a $(LFLAGS)

example4_c.out: example4_c.c
	$(ANSIC) $(CFLAGS) example4_c.c -o example4_c.out ../libcpfapack.a $(LFLAGS)

example5_c.out: example5_c.c
	$(ANSIC) $(CFLAGS) example5_c.c -o example5_c.out ../libcpfapack.a $(LFLAGS)

example6_c.out: example6_c.c
	$(ANSIC) $(CFLAGS) example6_c.c -o example6_c.out ../libcpfapack.a $(LFLAGS)

clean:
	rm -f *.out

