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
 *   Write program name etc to stdout                          
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void info(FILE *fp){

	fprintf(fp," ***********************************************************\n");
	fprintf(fp," This is program IFOS2D. Version 2.0.3                      \n");
	fprintf(fp," Parallel 2-D elastic Full Waveform Inversion code.         \n");
	fprintf(fp,"                                                            \n");
	fprintf(fp," ***********************************************************\n");
	fprintf(fp,"\n");

}
