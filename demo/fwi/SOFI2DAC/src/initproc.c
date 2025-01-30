/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2016  For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS.
 * 
 * IFOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * IFOS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IFOS. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
-----------------------------------------------------------------------------------------*/

/* This is function initproc.
   Dividing the 3-D FD grid into domains and assigning the
   PEs to these domains,
*/

#include "fd.h"

void initproc(void)	{

	extern int NX, NY, IENDX, IENDY, POS[3], INDEX[5];
	extern int NP, NPROC, NPROCX, NPROCY, MYID;
	extern FILE *FP;

	if (NPROC != NP)
		declare_error("Number of processors specified in the parameter file \n and at command line (NP) differ !");


	/*C-- determine the length of the subarray on this processor*/
	IENDX = NX/NPROCX;
	IENDY = NY/NPROCY;

	/* POS(1) indicates x POSition of the processor in the 
		     logical 3D processor array*/
	if ((NX%NPROCX)>0)
		declare_error(" NX%NPROX must be zero  !");
	if ((NY%NPROCY)>0)
		declare_error(" NY%NPROY must be zero  !");


	if (MYID==0){
		fprintf(FP,"\n **Message from initprocs (printed by PE %d):\n",MYID);
		fprintf(FP," Size of subarrays in gridpoints:\n");
		fprintf(FP," IENDX= %d\n",IENDX);
		fprintf(FP," IENDY (vertical) = %d\n",IENDY);
	}



	/*---------------   index is indicating neighbouring processes	--------------------*/
	INDEX[1]=MYID-1;  		 /* left	*/
	INDEX[2]=MYID+1;  		 /* right	*/
	INDEX[3]=MYID-NPROCX;  		 /* upper	*/
	INDEX[4]=MYID+NPROCX;  		 /* lower	*/
	
	/*---------------   POS indicates the processor location in the 3D logical processor array	---------*/
	POS[1] = MYID % NPROCX;			/*  x coordinate */
	POS[2] = (MYID/NPROCX); 	/*  y coordinate */

	if (POS[1] == 0)        INDEX[1]=INDEX[1] + NPROCX;        	  
	if (POS[1] == NPROCX-1) INDEX[2]=INDEX[2] - NPROCX;          	 
	if (POS[2] == 0)        INDEX[3]=(NPROCX*NPROCY)+MYID-NPROCX; 	 
	if (POS[2] == NPROCY-1) INDEX[4]=MYID+NPROCX-(NPROCX*NPROCY);	 

	fprintf(FP,"\n");
	fprintf(FP," **Message from initprocs (written by PE %d):\n",MYID);
	fprintf(FP," Processor locations in the 2D logical processor array\n");
	fprintf(FP," MYID \t POS(1):left,right \t POS(2): top, bottom\n");
	
	fprintf(FP," %d \t\t %d: %d,%d \t\t %d: %d,%d \n",
	    MYID,POS[1],INDEX[1],INDEX[2], POS[2], INDEX[3],INDEX[4]);
}
