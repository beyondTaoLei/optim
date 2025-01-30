/*---------------------------------------------------------------------------------
 * Copyright (C) 2015 For the list of authors, see file AUTHORS.
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
---------------------------------------------------------------------------------*/

/* $Id: cpml_update.c 819 2015-04-17 11:07:06Z tmetz $ */

/*--------------------------------------------------------------------------------
 * CPML update functions (ABS=1) for the stress and particle velocitiy components used
 * in the update functions of the CPML-boundary area.
 *    
 * -------------------------------------------------------------------------------*/


#include "fd.h"


/* CPML Functions for update_s ---------------------------------------------------*/

void cpml_update_p_x ( int i, int j, float  * vxx, float ** psi_vxx,
                       float * K_x, float * a_x, float * b_x)
{

	psi_vxx[j][i] = b_x[i] * psi_vxx[j][i] + a_x[i] * ( *vxx );
	*vxx = ( *vxx ) / K_x[i] + psi_vxx[j][i];
}

void cpml_update_p_y ( int i, int j, float * vyy, float ** psi_vyy,
                       float * K_y, float * a_y, float * b_y)
{

	psi_vyy[j][i] = b_y[j] * psi_vyy[j][i] + a_y[j] * ( *vyy );
	*vyy = ( *vyy ) / K_y[j] + psi_vyy[j][i];
        
}


/* CPML Functions for update_v ---------------------------------------------------*/

void cpml_update_v_x ( int i, int j,float  * sp_x, float ** psi_sxx_x,
                       float * K_x_half, float * a_x_half, float * b_x_half)
{

	psi_sxx_x[j][i] = b_x_half[i] * psi_sxx_x[j][i] + a_x_half[i] * ( *sp_x );
	*sp_x = ( *sp_x ) / K_x_half[i] + psi_sxx_x[j][i];

}

void cpml_update_v_y ( int i, int j,float * sp_y, float ** psi_syy_y,
                       float * K_y_half, float * a_y_half, float * b_y_half)
{

	psi_syy_y[j][i] = b_y_half[j] * psi_syy_y[j][i] + a_y_half[j] * ( *sp_y );
	*sp_y = ( *sp_y ) / K_y_half[j] + psi_syy_y[j][i];

}

