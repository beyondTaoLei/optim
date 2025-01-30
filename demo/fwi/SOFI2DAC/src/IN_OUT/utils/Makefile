# Makefile for SOFI2D

#--------------------------------------------------------
# edit here:


# LINUX with OpenMPI / IntelMPI and GCC Compiler
# LINUX with OpenMPI / IntelMPI and GCC Compiler
#CC = mpicc
#FC = mpif90
#CFLAGS = -O3

#GCC
#CC = mpicc
#FC = mpif90
#FFLAGS =  -O3 -g
#CFLAGS =  -O3 -g

#INTEL
#CC = mpiicc
#FC = mpiifort
#FFLAGS =  -nofor_main -O3 -g 
#CFLAGS =  -O3 -g

# after this line, no further editing should be necessary
# --------------------------------------------------------
include ../../Makefile
# variables for librairies

CC0 = $(CC)
CFLAGS0 = $(CFLAGS)
LFLAGS0 = $(LFLAGS)

INC_PATH = -I ../include/
LIB_PATH = -L../lib -L$(SACHOME)/lib -lIN_OUT -lsac -lsacio -lm -lfftw3

CXX_SOURCES= $(wildcard *.c)
CXX_OBJECTS = $(CXX_SOURCES:.c=)

% : %.c
	$(CC0) $(CFLAGS0) -o $@ $< $(INC_PATH) $(LIB_PATH)

donothing : # for giving the information when typing make
	@echo please type -make bin-
	
run : $(CXX_OBJECTS)
	@echo
	@echo $^
	mv $^ ../bin
	@echo
	
all : clean lib

clean :
	cd ../bin; rm -f $(CXX_OBJECTS)


