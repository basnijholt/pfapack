OBJECTS= skpfa.o skpf10.o skbpfa.o skbpf10.o sktrf.o sktrd.o skbtrd.o

CC = gcc
CFLAGS = -O3

all: $(OBJECTS)
	$(AR) rvu libcpfapack.a $(OBJECTS)

clean:
	rm -f *.o
	rm -f libcpfapack.a

.SUFFIXES: .o .c

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

