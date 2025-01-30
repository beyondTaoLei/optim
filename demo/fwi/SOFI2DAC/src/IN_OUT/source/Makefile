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
CXX_SOURCES_COMMON = $(wildcard $(COM)*.c)
		
CXX_SOURCES_kernel = \
		coord_grid.c \
		insert_src.c \
		read_coordinates.c \
		seis_glob2loc.c \
		sources.c \
		source_subfix.c \
		splitpos.c \
		src_term_shot.c \
		src_term_time.c \
		su_operators.c \
		sac_operators.c\
		wavelet.c \
		write_snap_boundary.c
		

CXX_WAVELET = write_wavelet.c
		
		
CXX_SOURCES= $(CXX_SOURCES_COMMON) $(CXX_SOURCES_kernel)
CXX_OBJECTS = $(CXX_SOURCES:.c=.o)	

CXX_WAVELET_SRC = $(CXX_SOURCES_COMMON) $(CXX_WAVELET)
CXX_WAVELET_OBJ = $(CXX_WAVELET_SRC:.c=.o)

%.o : %.c
	$(CC0) $(CFLAGS0) $(INC_PATH) -c -o $@ $<
	
donothing : # for giving the information when typing make
	@echo please type -make lib- for compiling the library

lib : $(CXX_OBJECTS)
	@echo "***LIB***"
	cd $(COM); $(MAKE) lib; cd -
	$(AR) $(ARFUNCT) $(LIB_PATH)libSource.a $(CXX_OBJECTS)
	
write_wavelet: $(CXX_WAVELET_OBJ)
	$(CC) $(CXX_WAVELET_OBJ) -o ../bin/write_wavelet $(LIB_PATH)libSource.a

all : clean lib

clean :
	find . -name "*.o" -exec rm {} \;
	find . -name "*.c%" -exec rm {} \;
	find $(COM) -name "*.o" -exec rm {} \;
	find $(LIB_PATH) -name "libComm.a" -exec rm {} \;
	find $(LIB_PATH) -name "libSource.a" -exec rm {} \;


