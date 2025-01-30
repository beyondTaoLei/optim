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
 *   zero wavefield
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void zero_acoustic(int ny1, int ny2, int nx1, int nx2, float **vx, float **vy, float **sp){

    M_init(nx1, nx2, ny1, ny2, vx);
    M_init(nx1, nx2, ny1, ny2, vy);
    M_init(nx1, nx2, ny1, ny2, sp);

}
