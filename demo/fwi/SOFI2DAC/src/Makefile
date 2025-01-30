# Makefile for SOFI2D

#--------------------------------------------------------
# edit here:

EXEC= ../bin

# Compiler (LAM: CC=hcc, CRAY T3E: CC=cc)

# ON Linux cluster running LAM
#CC=hcc
#LFLAGS=-lm -lmpi 
#CFLAGS=-Wall -O4

# ON Linux cluster running OpenMPI and ON MAC
CC=mpicc
#LFLAGS=-lm  
CFLAGS=-Wall -O3 

# On CRAY T3E
# CC=cc

# On SCHARnet system
#CC=mpicc
#LFLAGS=-lm  

# On HLRN system
#CC=mpcc
#LFLAGS=-lm  

# ALTIX
#CC=icc
#CFLAGS=-mp -O3 -ip0
#LFLAGS=-lmpi -lm 

IN_OUT= ./IN_OUT
IN_OUT_INC = $(IN_OUT)/include
IN_OUT_LIB_PATH = -L $(IN_OUT)/lib
IN_OUT_LIBS = -lIN_OUT

INC_PATH = -I $(IN_OUT_INC)
LIBS_PATH = $(IN_OUT_LIB_PATH)
LIBS = $(IN_OUT_LIBS) -lm
# after this line, no further editing should be necessary
# --------------------------------------------------------

SNAPMERGE_SCR = \
	snapmerge.c \
	json_parser.c\
	merge.c \
	read_par_json.c \
	readdsk.c \
	writedsk.c \
	util.c 

SNAPMERGE2_SCR = \
	snapmerge_boundary.c \
	json_parser.c\
	read_par_json.c \
	readdsk.c \
	writedsk.c \
	util.c 

SOFI2DAC_COMM= \
	av_rho.c \
	av_tau.c \
	calc_grid_for_inv.c \
	catseis.c \
	checkfd.c \
	comm_ini.c \
	cpml_update.c \
	exchange_par.c \
	exchange_p.c \
	exchange_v.c \
	extract_current_sources.c \
	holbergcoeff_all.c \
	holbergcoeff.c \
	info.c \
	initproc.c \
	json_parser.c \
	manage_grad_term.c \
	manage_seis_term.c \
	matcopy_acoustic.c \
	matcopy_viscac.c \
	matrix_operations.c \
	mergemod.c \
	note.c \
	outseis_glob.c \
	PML_pro.c \
	prepare_update_p.c \
	psource.c \
	readdsk.c \
	readmod_acoustic.c \
	readmod_viscac.c \
	read_par_json.c \
	receiver.c \
	receiver_boundary.c \
	saveseis_glob.c \
	seismo_ssg.c \
	snap_ssg.c \
	snap_ssg_freq.c \
	splitrec.c \
	subgrid_bounds.c \
	surface_acoustic_PML.c \
	time_count.c \
	update_p_interior.c \
	update_p_PML.c \
	update_p_visc_PML.c \
	update_v_interior.c \
	update_v_PML.c \
	util.c \
	writedsk.c \
	writemod.c \
	write_par.c \
	zero_acoustic.c \
	zero_fdveps_ac.c \
	zero_fdveps_viscac.c
	
SOFI2DAC_SRC= \
	$(SOFI2DAC_COMM) \
	sofi2Dac.c 
	

SOFI2DAC_BACK_SRC= \
	$(SOFI2DAC_COMM) \
	source_boundary.c \
	sofi2Dac_back.c
	

SOFI2DAC_GRAD_SRC= \
	$(SOFI2DAC_COMM) \
	source_boundary.c \
	sofi2Dac_grad.c
	
	
SNAPMERGE_OBJ = $(SNAPMERGE_SCR:%.c=%.o)
SNAPMERGE2_OBJ = $(SNAPMERGE2_SCR:%.c=%.o)

SOFI2DAC_OBJ = $(SOFI2DAC_SRC:%.c=%.o)
SOFI2DAC_BACK_OBJ = $(SOFI2DAC_BACK_SRC:%.c=%.o)
SOFI2DAC_GRAD_OBJ = $(SOFI2DAC_GRAD_SRC:%.c=%.o)

.c.o:
	$(CC) $(CFLAGS) -c $< $(INC_PATH)
	
lib:
	cd $(IN_OUT)  ; $(MAKE) lib
	
snapmerge:	$(SNAPMERGE_OBJ)
	$(CC) $(LFLAGS) $(SNAPMERGE_OBJ) -o ../bin/snapmerge $(LIBS_PATH) $(LIBS)

snapmerge2:	$(SNAPMERGE2_OBJ)
	$(CC) $(CFLAGS) $(SNAPMERGE2_OBJ) -o ../bin/snapmerge2 $(LIBS_PATH) $(LIBS)

sofi2Dac:		lib $(SOFI2DAC_OBJ)
	$(CC) $(SOFI2DAC_OBJ) -o $(EXEC)/sofi2Dac $(LIBS_PATH) $(LIBS)
	
sofi2Dac_back:		$(SOFI2DAC_BACK_OBJ)
	$(CC) $(SOFI2DAC_BACK_OBJ) -o $(EXEC)/sofi2Dac_back $(LIBS_PATH) $(LIBS)
	
sofi2Dac_grad:		$(SOFI2DAC_GRAD_OBJ)
	$(CC) $(SOFI2DAC_GRAD_OBJ) -o $(EXEC)/sofi2Dac_grad $(LIBS_PATH) $(LIBS)


clean:
	find . -name "*.o" -exec rm {} \; 
	find . -name "*.c%" -exec rm {} \;
	find . -name "*.bck" -exec rm {} \;

all: clean sofi2Dac snapmerge
