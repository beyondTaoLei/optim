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

AR= ar
ARFUNCT= crs

CC0 = $(CC)
CFLAGS0 = $(CFLAGS)
LFLAGS0 = $(LFLAGS)

INC_PATH = -I ../include/
LIB_PATH = ../lib/

COM = ../COMMON/

# printf '$(COM)%s \\\n' $(dir -1 *.c)
CXX_SOURCES_COMMON = \
		$(COM)envelope.c \
		$(COM)freqdomain_filt_vector.c \
		$(COM)interpolation.c \
		$(COM)mat_operations.c \
		$(COM)project.c \
		$(COM)util.c \
		$(COM)vector_operations.c \
		$(COM)vector_operations2.c \
		$(COM)windows.c \
		$(COM)xcorr.c \
		$(COM)fft_use.c \
		$(COM)json_parser.c \
		$(COM)string_opts.c \
		$(COM)update_shot_num.c

CXX_SOURCES= $(CXX_SOURCES_COMMON)
CXX_OBJECTS = $(CXX_SOURCES:.c=.o)	    

%.o : %.c
	$(CC0) $(CFLAGS0) $(INC_PATH) -c -o $@ $<
	
donothing : # for giving the information when typing make
	@echo please type -make lib- for compiling the library

lib : $(CXX_OBJECTS)
	@echo "***LIB***"
	$(AR) $(ARFUNCT) $(LIB_PATH)libComm.a $(CXX_OBJECTS)

all : clean lib

clean :
	find . -name "*.o" -exec rm {} \;
	find . -name "*.c%" -exec rm {} \;
	find $(LIB_PATH) -name "libComm.a" -exec rm {} \;


