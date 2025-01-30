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
 * globvar.h - global variables of viscoelastic 2D FD program
 * generally, for the names of the global variables uppercase letters are used
 *
 * ----------------------------------------------------------------------*/


float XS, YS, DH=0.0, TIME=0.0, DT=0.0, TS=0.0, DAMPING=0.0, PLANE_WAVE_DEPTH=0.0, PLANE_WAVE_ANGLE=0.0, SRCPOSXYZ[3]={0.0, 0.0, 0.0};
float TSNAP1=0.0, TSNAP2=0.0, TSNAPINC=0.0, *FL=NULL, TAU=0.0;
float XREC1=0.0, XREC2=0.0, YREC1=0.0, YREC2=0.0;
float REC_ARRAY_DEPTH=0.0, REC_ARRAY_DIST=0.0;
float REFREC[4]={0.0, 0.0, 0.0, 0.0};
float NGEOPH; /*in auto mode NGEOPH will be calculated from model dimensions, type integer is incorrect*/
int   RUNMODE, WRITE_MODELFILES=0, ABS_TYPE;
int   OUTNTIMESTEPINFO=1; /*every OUTNTIMESTEPINFO th timestep, information on the time step will be given to screen/file */
int   SEISMO=0, NDT=1, NSRC=1, SEIS_FORMAT=0, FREE_SURF=0, READMOD=0, READREC=0, SRCREC=0, RSG=0, FW=0;
int   NX=1, NY=1, NT=0, SOURCE_TYPE=0, SOURCE_SHAPE=0, SNAP=0, SNAP_FORMAT=0, LOG=0, REC_ARRAY=0;
int   L=0, BOUNDARY=0, DC=0, DRX=0, NXG=0, NYG=0, IDX=1, IDY=1, CHECKPTREAD=0, CHECKPTWRITE=0, FDORDER=0, FDORDER_TIME=0, MAXRELERROR=0;
int   RUN_MULTIPLE_SHOTS=0; /* Added for multiple shots */
char  SNAP_FILE[STRING_SIZE]="", SOURCE_FILE[STRING_SIZE]="", SIGNAL_FILE[STRING_SIZE]="";
char  MFILE[STRING_SIZE]="", REC_FILE[STRING_SIZE]="", CHECKPTFILE[STRING_SIZE]="";
char  SEIS_FILE[STRING_SIZE]="", LOG_FILE[STRING_SIZE]="";
FILE *FP;


/* MPI-variables */
int   NP=0, NPSP=0, NPROC=1, NPROCX=1, NPROCY=1, MYID, IENDX=0, IENDY=0;
int   POS[3], INDEX[5];     
const int TAG1=1,TAG2=2, TAG3=3, TAG4=4, TAG5=5,TAG6=6; 

/* spatial adaptive Code variables*/
int   check_id, cfgt_id, cfgt, jumpid;
float DH1;

/* PML-Parameters */
float NPOWER, K_MAX_CPML, FPML, VPPML;

/*******************************************************/
/* modified by Tao, in order to insert different source 
    * for FWI, RTM or some other cases */
/*******************************************************/
int NSHOT1, NSHOT2, NSHOT_INC, FLAG_SU_BIN, SRC_SNAP_TYPE;
int ADD_REPLACE_VX, ADD_REPLACE_VY, ADD_REPLACE_SP, ADD_REPLACE_SNAP;
int SAVE_LAST_SNAP, READ_LAST_SNAP, NX1, NX2, NY1, NY2, NXI, NYI, SEISMO2;
float MEMORY_DEPTH, SHOT_SPACING_X, SHOT_SPACING_Y;
char SRC_POS_FILE[STRING_SIZE]="", SNAP_COORDS_FILE[STRING_SIZE]="";
char SRC_F_VX[STRING_SIZE]="", SRC_F_VY[STRING_SIZE]="";
char SRC_F_SP[STRING_SIZE]="", SRC_F_SNAP[STRING_SIZE]="";
char FREQ_FILE[STRING_SIZE]="";

/* For the gradient */
int FORWARD_ONLY, INV_VP, INV_VS, INV_RHO, EPRECOND, EPRECOND_PER_SHOT;
int OMEGA1, OMEGA2, OMEGA_INC;
float EPSILON, TSHIFT;