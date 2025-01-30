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

/*------------------------------------------------------------------------
 *   stress free surface condition
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void surface_acoustic_PML(int ndepth, float ** sp){


	int i,j, m;
	int fdoh;
	extern int NX;
	extern int FDORDER;
	
	fdoh = FDORDER/2;

	j=ndepth;     /* The free surface is located exactly in y=1/2*dh !! */
	for (i=1;i<=NX;i++){
        sp[j][i] = 0;
		for (m=1; m<=fdoh; m++) {
			sp[j-m][i] = -sp[j+m][i];
		}
	}
}
