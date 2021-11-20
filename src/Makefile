# Makefile for primary src
.SUFFIXES:
.SUFFIXES: .o .f90 .c

# Compiler (ifort, gfortran)
FC = ifort
CC = icc

# Libraries
COMPILER = $(shell $(FC) --version | head -n1 | cut -d' ' -f1)

# ifort and gfortran take different compiler flags
ifeq ($(COMPILER),ifort)
   # Intel
   COMMONFLAGS = -qopenmp
   PRODFLAGS = -O3
endif

ifeq ($(COMPILER),GNU)
   # gfortran
   COMMONFLAGS = -ffree-line-length-0
   PRODFLAGS = -O2
endif

COMPFLAGS = $(COMMONFLAGS) $(PRODFLAGS)

# Libraries 
LIBS =  ${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

INC = -I/home/markherndon/vortex-stability/GSL_INTERFACE/ -I${MKLROOT}/include/intel64/ilp64 -I"${MKLROOT}/include"
# Object list
OBJECTS = mod_global.o                 \
	  mod_numerical_routines.o     \
	  mod_file_io.o                \
	  main.o

.c.o:; $(CC) $(CFLAGS) -c -o $@ $<
.f.o:; $(FC) $(COMPFLAGS) -c -o $@ $<
.f90.o:; $(FC) $(COMPFLAGS) -c -o $@ $< $(INC) $(LIBS)

src_objs: $(OBJECTS)

clean:
	rm -rf *.o *.mod $(EXEC_NAME)

## Object dependencies
mod_global.o: mod_global.f90
mod_numerical_routines.o: mod_numerical_routines.f90
mod_file_io.o: mod_file_io.f90 mod_global.o
main.o: main.f90 mod_file_io.o mod_global.o mod_numerical_routines.o
# Makefile:1 ends here
