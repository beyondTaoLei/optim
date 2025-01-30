#MASTER Makefile for the OPTIMIZATION
#include ../Makefile.inc
#include ../Makefile
# variables for librairies

AR= ar
ARFUNCT= crs

LIB_PATH = ./lib/
COM = ./COMMON/
KER_source = ./source/
KER_grad = ./grad/
KER_adj = ./adj_src/
KER_process = ./process/
KER_par = ./par/
KER_ssg = ./ssg/
KER_util = ./utils/

all : compile

lib : compile

compile : 
	cd $(KER_source); 	$(MAKE) lib
	cd $(KER_grad); 	$(MAKE) lib
	cd $(KER_adj); 		$(MAKE) lib
	cd $(KER_process); 	$(MAKE) lib
	cd $(KER_par);   	$(MAKE) lib
	cd $(KER_ssg);   	$(MAKE) lib
	$(AR) $(ARFUNCT) $(LIB_PATH)libIN_OUT.a $(COM)*.o $(KER_source)*.o $(KER_grad)*.o $(KER_adj)*.o $(KER_process)*.o $(KER_par)*.o $(KER_ssg)*.o
	# cd $(KER_util);  	$(MAKE) run
clean :	
	cd $(KER_source); 	$(MAKE) clean
	cd $(KER_grad); 	$(MAKE) clean
	cd $(KER_adj); 		$(MAKE) clean
	cd $(KER_process); 	$(MAKE) clean
	cd $(KER_par); 		$(MAKE) clean
	cd $(KER_ssg); 		$(MAKE) clean
	cd $(KER_util); 	$(MAKE) clean
	rm $(LIB_PATH)*.a

