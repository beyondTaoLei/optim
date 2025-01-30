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
 *   Read acoustic model properties (vp,density) from files
 *
 *  ----------------------------------------------------------------------*/


/* This file contains function readmod, which has the purpose
 to read data from model-files for viscoelastic simulation */

#include "fd.h"

void readmod_acoustic(float  **  rho, float **  pi){
    
    extern int NX, NY, NXG, NYG,  POS[3], MYID;//, PARAMETERIZATION;
    extern char  MFILE[STRING_SIZE];
    extern FILE *FP;
    
    
    /* local variables */
    float rhov, piv, vp;
    int i, j, ii, jj;
    FILE *fp_vp, *fp_rho = NULL;
    char filename[STRING_SIZE];
    
    
    fprintf(FP,"\n...reading model information from model-files...\n");

    /* read density and seismic velocities */
    /* ----------------------------------- */
// 	   if(PARAMETERIZATION==1){
    fprintf(FP,"\t Vp:\n\t %s.vp\n\n",MFILE);
    sprintf(filename,"%s.vp",MFILE);
    fp_vp=fopen(filename,"r");
    if (fp_vp==NULL) declare_error(" Could not open model file for Vp ! ");
    
    fprintf(FP,"\t Density:\n\t %s.rho\n\n",MFILE);
    sprintf(filename,"%s.rho",MFILE);
    fp_rho=fopen(filename,"r");
    if (fp_rho==NULL) declare_error(" Could not open model file for densities ! ");
//        }
	   
    
    /* loop over global grid */
    for (i=1;i<=NXG;i++){
        for (j=1;j<=NYG;j++){
            
            if(feof(fp_vp) && feof(fp_rho)){
                declare_error("Model file VP or RHO is to small. Check dimensions NX*NY of file.");
            }
            
            fread(&vp, sizeof(float), 1, fp_vp);
            fread(&rhov, sizeof(float), 1, fp_rho);
            
            if ((isnan(vp)) && (MYID==0)) {
            declare_error(" Found NaN-Values in Vp-Model !");}

            if ((isnan(rhov)) && (MYID==0)) {
            declare_error(" Found NaN-Values in Rho-Model !");}


            piv=vp*vp*rhov;
            
            /* only the PE which belongs to the current global gridpoint
             is saving model parameters in his local arrays */
            if ((POS[1]==((i-1)/NX)) &&
                (POS[2]==((j-1)/NY))){
                ii=i-POS[1]*NX;
                jj=j-POS[2]*NY;
                
                rho[jj][ii]=rhov;    
                pi[jj][ii]=piv;
                
                
            }
        }
    }
    
    
    fread(&vp, sizeof(float), 1, fp_vp);
    if(!feof(fp_vp)){
        declare_error("Model file VP is to big. Check dimensions NX*NY of file.");
    }
    fclose(fp_vp);
    
    fread(&rho, sizeof(float), 1, fp_rho);
    if(!feof(fp_rho)){
        declare_error("Model file RHO is to big. Check dimensions NX*NY of file.");
    }
    fclose(fp_rho);
    
//     /* Write model to file */
//     sprintf(filename,"%s.out.lam",MFILE);
//     write_matrix_disk(pi, filename);
//     
//     sprintf(filename,"%s.out.rho",MFILE);
//     write_matrix_disk(rho, filename);
    
}
