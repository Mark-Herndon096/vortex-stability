# Commond source directories
SRC_DIR = src
GSL_DIR = GSL_INTERFACE

# Compiler (ifort, gfortran)
FC = ifort
CC = icc

# Libraries
COMPILER = $(shell $(FC) --version | head -n1 | cut -d' ' -f1)

# These flags needed for GNU GSL library -- path dependent on your system
INCLUDE   = -I/custom_builds/GSL/include
LDFLAGS   = -L/custom_builds/GSL/lib -lgsl -lgslcblas -lm

# Libraries 
LIBS = ${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

INC =  -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include

LIBRARIES = $(INCLUDE) $(INC) $(LDFLAGS) $(LIBS)
# ifort and gfortran take different compiler flags
ifeq ($(COMPILER),ifort)
   # Intel
   COMMONFLAGS = -r8 -traceback -qopenmp
   PRODFLAGS = -O3
endif


# Set flags for debug or release
ifeq ($(MAKECMDGOALS),debug)
   COMPFLAGS = ${COMMONFLAGS} ${DEBUGFLAGS}
else
   COMPFLAGS = ${COMMONFLAGS} ${PRODFLAGS}
endif
# Subsection for executable variables and objects


# Executable name
EXEC_NAME = vortex_solver.exe

# Object list
OBJECTS = $(GSL_DIR)/special_function_wrapper.o   \
	  $(GSL_DIR)/special_function_interface.o \
	  $(SRC_DIR)/mod_global.o                 \
	  $(SRC_DIR)/mod_numerical_routines.o     \
	  $(SRC_DIR)/mod_file_io.o                \
	  $(SRC_DIR)/main.o
solver:
	$(MAKE) -C $(GSL_DIR) gsl_objs
	$(MAKE) -C $(SRC_DIR) src_objs
	$(FC) -o $(EXEC_NAME) $(COMPFLAGS) $(OBJECTS) $(LIBRARIES)

clean:
	$(MAKE) -C $(GSL_DIR) clean
	$(MAKE) -C $(SRC_DIR) clean
	rm $(EXEC_NAME)
# Makefile:1 ends here
