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
 *  compute receiver positions or read them from external file
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include <stdbool.h>


int **receiver(FILE *fp, int *ntr){

    /* declaration of extern variables */
    extern  char REC_FILE[STRING_SIZE];
    extern float DH, REFREC[4], REC_ARRAY_DEPTH, REC_ARRAY_DIST;
    extern int READREC, DRX, REC_ARRAY, NXG, FW;
    /*extern int NGEOPH;*/
    extern float NGEOPH;
    extern int MYID;
    extern float XREC1, YREC1, XREC2, YREC2;

    int **recpos1, **recpos=NULL;
    int   itr=1, itr1=0, itr2=0, recflag=0, n, i, j;
    float nxrec=0, nyrec=0;
    float nxrec1, nxrec2, nyrec1, nyrec2;
    float xrec, yrec;
    bool testbuff1, testbuff2, testbuff3;
    char bufferstring[10], buffer[STRING_SIZE];

    FILE *fpr;


    if (MYID==0){
        fprintf(fp,"-------------------- RECEIVER POSITIONS ---------------------\n");
        fprintf(fp," Message from function receiver (written by PE %d):\n\n",MYID);
        if (READREC){ /* read receiver positions from file */
                fprintf(fp," Reading receiver positions from file: '%s'\n",REC_FILE);
                fpr=fopen(REC_FILE,"r");
                if (fpr==NULL) declare_error(" Receiver file could not be opened !");
                *ntr=0;

                /* counts the number of receivers in the receiver file */
                while(fgets(buffer, STRING_SIZE, fpr)){
                    testbuff1=strchr(buffer,'#');
                    testbuff2=strchr(buffer,'%');
                    testbuff3=sscanf(buffer,"%s",bufferstring)==1;

                    /*the following output is for debugging*/
                    /*testbuff4=(testbuff1==1 || testbuff2==1);
                    fprintf(fp," buffer: _%s_with testbuff1=_%i_ testbuff2=_%i_testbuff3=_%i_ testbuff4=_%i_\n",buffer,testbuff1, testbuff2, testbuff3,testbuff4);*/
                    /* checks if the line contains a '%' or '#' character which indicates a
                    comment line, and if the reading of a string was successful, 
                    which is not the case for an empty line*/
                    if (((testbuff1==1 || testbuff2==1)==0) && testbuff3==1) ++(*ntr);
                }

                rewind(fpr);

                recpos1=imatrix(1,3,1,*ntr);
                for (itr=1;itr<=*ntr;itr++){
                    fscanf(fpr,"%f%f\n",&xrec, &yrec);
                    recpos1[1][itr]=iround((xrec+REFREC[1])/DH);
                    recpos1[2][itr]=iround((yrec+REFREC[2])/DH);
                    recpos1[3][itr]=itr;
                }
                fclose(fpr);
                fprintf(fp," Number of receiver positions found: %i\n",*ntr);

                /* check if more than one receiver is located
                                at the same gridpoint */
                for (itr=1;itr<=(*ntr-1);itr++)
                    for (itr1=itr+1;itr1<=*ntr;itr1++)
                        if ((recpos1[1][itr]==recpos1[1][itr1])
                            && (recpos1[2][itr]==recpos1[2][itr1]))
                                recpos1[1][itr1]=-(++recflag);

                recpos=imatrix(1,3,1,*ntr-recflag);
                for (itr=1;itr<=*ntr;itr++)
                    if (recpos1[1][itr]>0){
                        recpos[1][++itr2]=recpos1[1][itr];
                        recpos[2][itr2]=recpos1[2][itr];
                        recpos[3][itr2]=recpos1[3][itr];
                    }

                *ntr=itr2;
                if (recflag>0){
                    fprintf(fp,"\n\n");
                    fprintf(fp," Warning:\n");
                    fprintf(fp," Several receivers located at the same gridpoint !\n");
                    fprintf(fp," Number of receivers reduced to %i\n", *ntr);
                    fprintf(fp,"\n\n");
                }

                free_imatrix(recpos1,1,3,1,*ntr);

        }else if (REC_ARRAY>0){
            fprintf(fp," Generating receiver planes as specified in input file.\n");


            *ntr=(1+(NXG-2*FW)/DRX)*REC_ARRAY;
            recpos=imatrix(1,3,1,*ntr);
            itr=0;
            for (n=0;n<=REC_ARRAY-1;n++){
                    j=iround((REC_ARRAY_DEPTH+REC_ARRAY_DIST*(float)n)/DH);
                    for (i=FW;i<=NXG-FW;i+=DRX){
                            itr++;
                            recpos[1][itr]=i;
                            recpos[2][itr]=j;
                            recpos[3][itr]=itr;
                    }
            }
        }else{         /* straight horizontal or vertical line of receivers */
            fprintf(fp," Reading receiver positions from input file. \n");
            nxrec1=XREC1/DH;   /* (nxrec1,nyrec1) and (nxrec2,nyrec2) are */
            nyrec1=YREC1/DH;   /* the positions of the first and last receiver*/
            nxrec2=XREC2/DH;	 /* in gridpoints */
            nyrec2=YREC2/DH;


            /* only 1 receiver */
            if (((iround(nxrec2)-iround(nxrec1))==0) && ((iround(nyrec2)-iround(nyrec1))==0)) {
                    fprintf(fp," A single receiver position read from input file. \n");
                    *ntr = 1;
                    recpos = imatrix(1,3,1,*ntr);
                    recpos[1][1] = iround(nxrec1);
                    recpos[2][1] = iround(nyrec1);
                    recpos[3][1] = 1;
            } else if (((iround(nyrec2)-iround(nyrec1))==0)) {
                    fprintf(fp," A horizontal receiver line (in x-direction) is specified in the input file. \n");
                    *ntr=iround((nxrec2-nxrec1)/NGEOPH)+1;
                    recpos=imatrix(1,3,1,*ntr);
                    nxrec = nxrec1;
                    for (n=1;n<=*ntr;n++){
                            nyrec=nyrec1+((nyrec2-nyrec1)/(nxrec2-nxrec1)*(nxrec-nxrec1));
                            itr=iround((nxrec-nxrec1)/NGEOPH)+1;
                            recpos[1][itr] = iround(nxrec);
                            recpos[2][itr] = iround(nyrec);
                            recpos[3][itr] = itr;
                            nxrec += NGEOPH;
                    }

            }
            else if (((iround(nxrec2)-iround(nxrec1))==0)) {
                    fprintf(fp," A vertical receiver line (in y-direction) is specified in the input file. \n");
                    *ntr=iround((nyrec2-nyrec1)/NGEOPH)+1;
                    recpos=imatrix(1,3,1,*ntr);
                    nyrec = nyrec1;
                    for (n=1;n<=*ntr;n++){
                            nxrec=nxrec1+((nxrec2-nxrec1)/(nyrec2-nyrec1)*(nyrec-nyrec1));
                            itr=iround((nyrec-nyrec1)/NGEOPH)+1;
                            recpos[1][itr] = iround(nxrec);
                            recpos[2][itr] = iround(nyrec);
                            recpos[3][itr] = itr;
                            nyrec += NGEOPH;
                    }
            }
            else {
                    /* arbitrary geophone-line */
                    fprintf(fp," No horizontal or vertical receiver line is specified in the input file. \n");
                    fprintf(fp," In order to define an arbitrary receiver line, please make use of an external receiver file (READREC=1). \n");
                    declare_error(" Error in specifying receiver coordinates in the input file !");
            }
        } /* end of if receivers specified in input file */
        fprintf(fp," Number of receiver positions found: %i\n",*ntr);


        /*fprintf(fp,"\n Position of %i receiver(s) (x,y,z) in gridpoints:\n",*ntr);*/
        for (itr=1;itr<=*ntr;itr++) {
                /* check for 0's */
                for (n=1;n<=2;n++)
                        if (recpos[n][itr]==0)
                                recpos[n][itr] = 1;
                /*outputs receiver locations as a list, can be quite long if a large number of receivers is used */
                /*fprintf(fp,"   Receiver %i: (%i, %i)\n",itr,recpos[1][itr],recpos[2][itr]);*/
        }
        fprintf(fp,"-------------------------------------------------------------\n\n");
    } /* End of if(MYID==0) */


    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(ntr,1,MPI_INT,0,MPI_COMM_WORLD);
    if (MYID!=0) recpos=imatrix(1,3,1,*ntr);
    MPI_Bcast(&recpos[1][1],(*ntr)*3,MPI_INT,0,MPI_COMM_WORLD);

    if (MYID==0){
        fprintf(fp,"\n **Message from function receiver (written by PE %d):\n",MYID);
        fprintf(fp," Number of receiver positions found: %i\n",*ntr);

        if (*ntr>50) fprintf(fp," The following table is quite large (%i lines) and will, thus, be truncated to the first 50 entries! \n\n",*ntr);

        fprintf(fp," Receiver positions in the global model-system:\n");
        fprintf(fp," x (gridpoints) y (gridpoints) \t x (in m) \t y (in m) \n");
        fprintf(fp," -------------  -------------- \t ---------\t -------- \n");
        if (*ntr>50) { for (itr=1;itr<=50;itr++)
                fprintf(fp," %i\t\t %i \t\t %6.2f \t %6.2f \n",recpos[1][itr],recpos[2][itr],recpos[1][itr]*DH,recpos[2][itr]*DH);
        }
        else for (itr=1;itr<=*ntr;itr++)
                fprintf(fp," %i\t\t %i \t\t %6.2f \t %6.2f \n",recpos[1][itr],recpos[2][itr],recpos[1][itr]*DH,recpos[2][itr]*DH);

    }

    fprintf(fp,"\n\n");

    return recpos;
}
