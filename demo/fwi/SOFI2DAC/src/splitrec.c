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
/* ----------------------------------------------------------------------
 * Computation of local receiver coordinates
 * (within each subgrid)
 *	 
 * ----------------------------------------------------------------------*/
 
 #include "fd.h"
int **splitrec(int **recpos,int *ntr_loc, int ntr, int *recswitch)
{

	extern int IENDX, IENDY, POS[4];

	int a,b,i=0,j,k;
	//int found=0;
	int ** recpos_dummy, **recpos_local=NULL;
	recpos_dummy = imatrix(1,3,1,ntr);

	for (j=1;j<=ntr;j++)
	{
		recswitch[j] = 0;
		a=(recpos[1][j]-1)/IENDX;
		b=(recpos[2][j]-1)/IENDY;

		if ((POS[1]==a)&&(POS[2]==b))
		{
			//found = 1;
			recswitch[j] = 1;
			i++;	/* Anzahl der Empf?nger i jedes Prozesses ermitteln */
			recpos_dummy[1][i] = ((recpos[1][j]-1)%IENDX)+1;
			recpos_dummy[2][i] = ((recpos[2][j]-1)%IENDY)+1;
			recpos_dummy[3][i] = j;
		}
	}
   
	if (i>0) recpos_local = imatrix(1,3,1,i);
	for (k=1;k<=i;k++){
		recpos_local[1][k] = recpos_dummy[1][k];
		recpos_local[2][k] = recpos_dummy[2][k];
		recpos_local[3][k] = recpos_dummy[3][k];
	}
	free_imatrix(recpos_dummy,1,3,1,ntr);

/*	
	fprintf(FP,"\n **Message from split_rec:\n");
	fprintf(FP," Splitting of receivers from global to local grids finished.\n");
	fprintf(FP," MYID \t \t no. of receivers\n");
	fprintf(FP," %d \t\t %d\n",MYID,i);

	fprintf(FP,"\n **Message from split_rec:\n");
	fprintf(FP," Table of local receiver positions (in gridpoints):\n");
	fprintf(FP," MYID \t x \t y\n");

	for (j=1;j<=i;j++)
		fprintf(FP," %d \t %d \t %d\n",
		    MYID,recpos_local[1][j],recpos_local[2][j]);
*/
   	*ntr_loc=i;
	return recpos_local;

}
