
# LIB
CHOLMOD_LIB= -L/usr/lib/ -lcholmod -lcblas 

# INCLUDE
CHOLMOD_INC= -I/usr/include/suitesparse/

# COLLECT
INC= $(CHOLMOD_INC)
LIB= $(CHOLMOD_LIB)

# LAPACK AND BLAS
USE_LAPACK=-DNO_LAPACK

# Compiler flags
#FLAGS= -O3 -std=c++0x -m64
FLAGS= -O3 -m64 ${USE_LAPACK}

# Compilers
CC=gcc
CXX=g++ -std=c++0x
