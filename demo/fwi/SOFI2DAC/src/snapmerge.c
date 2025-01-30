/*------------------------------------------------------------------------
 * Copyright (C) 2011 For the list of authors, see file AUTHORS.
 *
 * This file is part of SOFI2D.
 * 
 * SOFI2D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * SOFI2D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with SOFI2D. See file COPYING and/or 
  * <http://www.gnu.org/licenses/gpl-2.0.html>.
--------------------------------------------------------------------------*/
/*------------------------------------------------------------------------
 *   loop over snapshotfiles which have to be merged.                                   
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "globvar.h"      /* definition of global variables  */

int main(int argc, char **argv){

int  nsnap;
//char *fileinp="", modestr[10], infostr[70];
char *fileinp="";

/* ============================================== */
/* Open parameter-file to check if auto mode or not*/
fileinp = argv[1];

printf(" ***********************************************************\n");
printf(" This is program SNAPMERGE. \n");
printf(" Merge of snapshot files from the parallel  \n 2-D Viscoelastic Finite Difference Modeling      \n");
printf("                                                            \n");
printf(" written by  T. Bohlen                          \n");
printf(" Geophysical Institute, Department of Physics,         \n");
printf(" Institute of Technology, Karlsruhe, Germany         \n");
printf(" http://www.gpi.kit.edu \n");
printf(" ***********************************************************\n");
printf("\n");
printf(" Syntax example if excecuted from ./par directory: ../bin/snapmerge in_and_out/sofi2D.json \n");
printf(" Input file for the snapmerge process from command line : %s \n",fileinp);

//FP = fopen(fileinp,"r");
if ((FP=fopen(fileinp,"r"))==NULL) declare_error(" Opening input file failed.");
else printf(" Opening input file was successful.\n\n");
fclose(FP);
//fscanf(FP, "%s %s = %i", infostr, modestr, &RUNMODE);

/* =================================================== */
if (RUNMODE == 0)
	/* read standard input file */
	if (strstr(fileinp,".json"))
		//read json formated input file
		read_par_json(stdout, fileinp);
/*	else
		/read "old" input file *.inp, might not work in future
		read_par(stdout, fileinp);
else
	 auto mode: read input files 
	read_par_auto(stdout, fileinp);*/
/* =================================================== */


NXG=NX;
NYG=NY;	
NX = NXG/NPROCX;
NY = NYG/NPROCY;

nsnap=1+iround((TSNAP2-TSNAP1)/TSNAPINC);

FP=stdout;

switch(SNAP){
case 1 : /*particle velocity*/
   merge(nsnap,1);
   merge(nsnap,2);
   break;
case 2 : /*pressure */
   merge(nsnap,6);
   break;
case 4 : /*particle velocity*/
   merge(nsnap,1);
   merge(nsnap,2);
   merge(nsnap,6);
case 3 :
   merge(nsnap,4);
   break;
default :
   warning(" snapmerge: cannot identify content of snapshot !");
   break;

}	
return 0;	

}
