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
 *   store amplitudes (particle velocities or pressure or curl and div) 
 *    at receiver positions in arrays
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void seismo_ssg(int lsamp, int ntr, int **recpos, float **sectionvx,
                float **sectionvy, float **sectionp, float **sectiondiv,
                float **vx, float **vy, float **sp, float **pi, float *hc){
    
    extern int NDT, SEISMO, FDORDER;
    extern float DH;
    extern FILE *FP;
    int i,j, itr, ins, nxrec, nyrec, m, fdoh;
    float dh24, dhi, vxx, vyy, vxy, vyx;
    
    
//     dh24=1.0/(DH*24.0);
    fdoh = FDORDER/2;
    dhi = 1.0/DH;
    
    ins=lsamp/NDT;
//     ins=lsamp;
    for (itr=1;itr<=ntr;itr++){
        nxrec=recpos[1][itr];
        nyrec=recpos[2][itr];
        switch (SEISMO){
            case 1 :
                sectionvx[itr][ins]=vx[nyrec][nxrec];
                sectionvy[itr][ins]=vy[nyrec][nxrec];
                break;
                
            case 2 :
                sectionp[itr][ins]=-sp[nyrec][nxrec];
                break;
                
            case 3 :
                i=nxrec; j=nyrec;
                
                vxx = 0;
                vyy = 0;
                vyx = 0;
                vxy = 0;
                for (m=1; m<=fdoh; m++) {
                    vxx += hc[m]*(vx[j][i+m-1] -vx[j][i-m]  );
                    vyy += hc[m]*(vy[j+m-1][i] -vy[j-m][i]  );
                    vyx += hc[m]*(vy[j][i+m]   -vy[j][i-m+1]);
                    vxy += hc[m]*(vx[j+m][i]   -vx[j-m+1][i]);
                }
                vxx *= dhi;
                vyy *= dhi;
                vyx *= dhi;
                vxy *= dhi;
                
                sectiondiv[itr][ins]=(vxx+vyy)*sqrt(pi[j][i]);
//                 sectioncurl[itr][ins]=0;
                break;
                
            case 4 :
                i=nxrec; j=nyrec;
                
                vxx = 0;
                vyy = 0;
                vyx = 0;
                vxy = 0;
                for (m=1; m<=fdoh; m++) {
                    vxx += hc[m]*(vx[j][i+m-1] -vx[j][i-m]  );
                    vyy += hc[m]*(vy[j+m-1][i] -vy[j-m][i]  );
                    vyx += hc[m]*(vy[j][i+m]   -vy[j][i-m+1]);
                    vxy += hc[m]*(vx[j+m][i]   -vx[j-m+1][i]);
                }
                vxx *= dhi;
                vyy *= dhi;
                vyx *= dhi;
                vxy *= dhi;
                
                sectiondiv[itr][ins]=(vxx+vyy)*sqrt(pi[j][i]);
//                 sectioncurl[itr][ins]=0;
                sectionvx[itr][ins]=vx[nyrec][nxrec];
                sectionvy[itr][ins]=vy[nyrec][nxrec];
                sectionp[itr][ins]=-sp[nyrec][nxrec];
                break;
                
            case 5 :
                sectionvx[itr][ins]=vx[nyrec][nxrec];
                sectionvy[itr][ins]=vy[nyrec][nxrec];
                sectionp[itr][ins]=-sp[nyrec][nxrec];
                
                break;
                
        }

	}
}

void seismo_ssg2(int lsamp, int ntr, int **recpos, float **sectionvx,
                 float **sectionvy, float **sectionp,float **vx, 
                 float **vy, float **sp){
    
    int itr, ins, nxrec, nyrec;
    
    ins=lsamp;
    for (itr=1;itr<=ntr;itr++){
        nxrec=recpos[1][itr];
        nyrec=recpos[2][itr];
        sectionvx[itr][ins]=vx[nyrec][nxrec];
        sectionvy[itr][ins]=vy[nyrec][nxrec];
        sectionp[itr][ins]=sp[nyrec][nxrec];
    }
}