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
 *   merge model files written by the different processes to 
 *   a single file                                 
 *   - If VERBOSE==1, then this function will be silent
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void mergemod(char modfile[STRING_SIZE], int format){



	extern int  NYG, MYID, NPROCX, NPROCY;
	extern int NX, NY, NPROC, IDX, IDY;//,VERBOSE;
	extern FILE *FP;


	char file[STRING_SIZE];
	FILE *fp[NPROCY][NPROCX], *fpout;
	int i, j, ip, jp;
	float a;

//     if(VERBOSE==1){
//         fprintf(FP,"\n **Message from mergemod (printed by PE %d):\n",MYID);
//         fprintf(FP," PE %d starts merge of %d model files \n",MYID,NPROC);
//         fprintf(FP,"\n writing merged model file to  %s \n",modfile);
//     }
	fpout=fopen(modfile,"w");

	
	
	
// 	if(VERBOSE==1) fprintf(FP," Opening model files: %s.??? ",modfile);
	for (ip=0;ip<=NPROCX-1; ip++)
  	for (jp=0;jp<=NPROCY-1; jp++){
      		sprintf(file,"%s.%i.%i",modfile,ip,jp);
      		fp[jp][ip]=fopen(file,"r");
      		if (fp[jp][ip]==NULL) {
                /* If it fails, try again one time */
                fprintf(FP,"\n\n Cannot merge model files: %s.???, try again... \n",modfile);
                fp[jp][ip]=fopen(file,"r");
                if (fp[jp][ip]==NULL) declare_error("merge: can't read model file !");
            }
    }

// 	if(VERBOSE==1) fprintf(FP," ... finished. \n");



// 	if(VERBOSE==1) fprintf(FP," Copying...");

  	 for (ip=0;ip<=NPROCX-1; ip++){
      		for (i=1;i<=NX;i+=IDX){
			for (jp=0;jp<=NPROCY-1; jp++){
	    			for (j=1;j<=NY;j+=IDY){
	      			a=readdsk(fp[jp][ip],format);
	      			writedsk(fpout,a,format);
	       			}
	   		}
	 	}
      	}
// 	if(VERBOSE==1) fprintf(FP," ... finished. \n");

	for (ip=0;ip<=NPROCX-1; ip++)
   	for (jp=0;jp<=NPROCY-1; jp++){
      		fclose(fp[jp][ip]);
      	}
    fclose(fpout);
//     if(VERBOSE==1){
//         fprintf(FP," Use \n");
//         fprintf(FP," ximage n1=%d < %s  label1=Y label2=X title=%s \n",
//                 ((NYG-1)/IDY)+1,modfile,modfile);
//         fprintf(FP," to visualize model. \n");
//     }
    
    

}


