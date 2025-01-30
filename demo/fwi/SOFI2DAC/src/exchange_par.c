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
 *   Exchange FD-Parameters between PEs                         
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void exchange_par(void){

	/* declaration of extern variables */
	extern int   NX, NY, FDORDER, MAXRELERROR, SOURCE_TYPE, SOURCE_SHAPE, SNAP, SNAP_FORMAT, L, FDORDER_TIME;
	extern float DH, TIME, DT, TS, *FL, TAU, DAMPING, FPML, NPOWER, K_MAX_CPML, VPPML, PLANE_WAVE_DEPTH, PLANE_WAVE_ANGLE, SRCPOSXYZ[3];
	extern float XREC1, XREC2, YREC1, YREC2;
	extern float REC_ARRAY_DEPTH, REC_ARRAY_DIST;
	/*extern int NGEOPH;*/
	extern float NGEOPH;
	extern int SEISMO, NDT, SEIS_FORMAT, FREE_SURF, READMOD, READREC, SRCREC;
	extern int BOUNDARY, REC_ARRAY, DRX, LOG, RSG, WRITE_MODELFILES;
	extern float TSNAP1, TSNAP2, TSNAPINC, REFREC[4];
	extern char  MFILE[STRING_SIZE], SIGNAL_FILE[STRING_SIZE], LOG_FILE[STRING_SIZE];
	extern char SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE], REC_FILE[STRING_SIZE];
	extern char SEIS_FILE[STRING_SIZE], CHECKPTFILE[STRING_SIZE];
	extern int NPROC, NPROCX, NPROCY, MYID, IDX, IDY, CHECKPTREAD, CHECKPTWRITE, RUN_MULTIPLE_SHOTS, ABS_TYPE, FW; 
        
        /*******************************************************/
        /* modified by Tao, in order to insert different source 
            * for FWI, RTM or some other cases */
        /*******************************************************/
        extern int NSHOT1, NSHOT2, NSHOT_INC, FLAG_SU_BIN, SRC_SNAP_TYPE;
        extern int ADD_REPLACE_VX, ADD_REPLACE_VY, ADD_REPLACE_SP, ADD_REPLACE_SNAP;
        extern int SAVE_LAST_SNAP, READ_LAST_SNAP, SEISMO2;
        extern float MEMORY_DEPTH, SHOT_SPACING_X, SHOT_SPACING_Y;
        extern char SNAP_COORDS_FILE[STRING_SIZE];
        extern char SRC_F_VX[STRING_SIZE], SRC_F_VY[STRING_SIZE];
        extern char SRC_F_SP[STRING_SIZE], SRC_F_SNAP[STRING_SIZE];
        extern char FREQ_FILE[STRING_SIZE];
        
        /* For the gradient */
        extern int FORWARD_ONLY, INV_VP, INV_VS, INV_RHO, EPRECOND, EPRECOND_PER_SHOT;
        extern int OMEGA1, OMEGA2, OMEGA_INC;
        extern float EPSILON, TSHIFT;

	/* definition of local variables */
	int idum[NPAR];
	float fdum[NPAR];


	if (MYID == 0){ 

		fdum[1]  = DH;
		fdum[2]  = TIME;
		fdum[3]  = DT;
		fdum[4]  = TS;
		fdum[5]  = 0.0;
		fdum[6]  = 0.0;

		fdum[8]  = TAU;
		
		fdum[10]  = TSNAP1;
		fdum[11]  = TSNAP2;
		fdum[12]  = TSNAPINC;
		fdum[13]  = REFREC[1];
		fdum[14]  = REFREC[2];
		fdum[15]  = PLANE_WAVE_ANGLE;

		fdum[16]  = XREC1;
		fdum[17]  = YREC1;

		fdum[19]  = XREC2;
		fdum[20]  = YREC2;

		fdum[22]  = DAMPING;
		fdum[23]  = REC_ARRAY_DEPTH;
		fdum[24]  = REC_ARRAY_DIST;
		fdum[25]  = PLANE_WAVE_DEPTH;

		fdum[26]  = NGEOPH;

		fdum[27]  = SRCPOSXYZ[0];
		fdum[28]  = SRCPOSXYZ[1];
		fdum[29]  = SRCPOSXYZ[2];
		
		fdum[30]  = FPML;
		fdum[31]  = VPPML;
		fdum[32]  = NPOWER;
		fdum[33]  = K_MAX_CPML;
                
                fdum[34]  = MEMORY_DEPTH;
                fdum[35]  = SHOT_SPACING_X;
                fdum[36]  = SHOT_SPACING_Y;
                
                fdum[37] = EPSILON;
                fdum[38] = TSHIFT;


		/*************************************/

		idum[1]  = NPROCX;
		idum[2]  = NPROCY;
		idum[3]  = LOG;

		idum[4]  = NPROCX*NPROCY;
		idum[5]  = NX;
		idum[6]  = NY;
		idum[7] = FW;
		idum[8]  = SOURCE_SHAPE;
		idum[9]  = SOURCE_TYPE;
		idum[10]  = READMOD;
		idum[11]  = L;
		idum[12]  = FREE_SURF;
		idum[13]  = SNAP;
		idum[14]  = DRX;

		idum[16]  = BOUNDARY;
		idum[17]  = REC_ARRAY;
		idum[18]  = SRCREC;
		idum[19]  = IDX;
		idum[20]  = IDY;
		idum[21]  = 0;
		idum[22]  = 0;
		idum[23]  = SNAP_FORMAT;
		idum[24]  = SEISMO;
		idum[25]  = READREC;
		idum[26]  = RSG;
		idum[27]  = NDT;
		idum[28]  = SEIS_FORMAT;
		idum[29]  = CHECKPTREAD;
		idum[30]  = CHECKPTWRITE;

		idum[31]  = FDORDER;
		idum[32]  = MAXRELERROR;
		idum[33]  = RUN_MULTIPLE_SHOTS;
		idum[34]  = WRITE_MODELFILES;
		
		idum[35] = ABS_TYPE;
        
                idum[36] = FDORDER_TIME;
                
                idum[37] = NSHOT1;
                idum[38] = NSHOT2;
                idum[39] = NSHOT_INC;
                idum[40] = SRC_SNAP_TYPE;
                idum[41] = ADD_REPLACE_VX;
                idum[42] = ADD_REPLACE_VY;
                idum[43] = ADD_REPLACE_SP;
                idum[44] = ADD_REPLACE_SNAP;
                idum[45] = FLAG_SU_BIN;
                idum[46] = SAVE_LAST_SNAP;
                idum[47] = READ_LAST_SNAP;
                
                idum[48] = FORWARD_ONLY;
                idum[49] = INV_VP;
                idum[50] = INV_VS;
                idum[51] = INV_RHO;
                idum[52] = EPRECOND;
                idum[53] = EPRECOND_PER_SHOT;
                idum[54] = OMEGA1;
                idum[55] = OMEGA2;
                idum[56] = OMEGA_INC;
                idum[57] = SEISMO2;

	} /** if (MYID == 0) **/

	if (MYID != 0) FL=vector(1,L);
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Bcast(&idum,NPAR,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&fdum,NPAR,MPI_FLOAT,0,MPI_COMM_WORLD);

// 	MPI_Bcast(&SOURCE_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&MFILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SNAP_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&REC_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SEIS_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&LOG_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
// 	MPI_Bcast(&SIGNAL_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&CHECKPTFILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
        
        MPI_Bcast(&SNAP_COORDS_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
        MPI_Bcast(&SRC_F_SP,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
        MPI_Bcast(&SRC_F_VX,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
        MPI_Bcast(&SRC_F_VY,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
        MPI_Bcast(&SRC_F_SNAP,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
        MPI_Bcast(&FREQ_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
        
	MPI_Barrier(MPI_COMM_WORLD);

	DH=fdum[1];
	TIME=fdum[2];
	DT=fdum[3];
	TS=fdum[4];

	TAU=fdum[8];
	
	TSNAP1=fdum[10];
	TSNAP2=fdum[11];
	TSNAPINC=fdum[12];
	REFREC[1]=fdum[13];
	REFREC[2]=fdum[14];
	PLANE_WAVE_ANGLE=fdum[15];
	XREC1=fdum[16];
	YREC1=fdum[17];

	XREC2=fdum[19];
	YREC2=fdum[20];

	DAMPING=fdum[22];
	REC_ARRAY_DEPTH=fdum[23];
	REC_ARRAY_DIST=fdum[24];
	PLANE_WAVE_DEPTH=fdum[25];

	NGEOPH = fdum[26];

	SRCPOSXYZ[0] = fdum[27];
	SRCPOSXYZ[1] = fdum[28];
	SRCPOSXYZ[2] = fdum[29];

	FPML = fdum[30];
	VPPML = fdum[31];
	NPOWER = fdum[32];
	K_MAX_CPML = fdum[33];
        
        MEMORY_DEPTH = fdum[34];
        SHOT_SPACING_X = fdum[35];
        SHOT_SPACING_Y = fdum[36];
        
        EPSILON = fdum[37];
        TSHIFT = fdum[38];

	/********************************************/

	NPROCX = idum[1];
	NPROCY = idum[2];
	LOG=idum[3];
	NPROC  = idum[4];
	NX = idum[5];
	NY = idum[6];
	FW = idum [7];
	SOURCE_SHAPE = idum[8];
	SOURCE_TYPE = idum[9];
	READMOD = idum[10];
	L = idum[11];
	FREE_SURF = idum[12];
	SNAP = idum[13];
	DRX = idum[14];

	BOUNDARY = idum[16];
	REC_ARRAY = idum[17];
	SRCREC = idum[18];
	IDX = idum[19];
	IDY = idum[20];


	SNAP_FORMAT = idum[23];
	SEISMO = idum[24];
	READREC = idum[25];
	RSG = idum[26];
	NDT = idum[27];
	SEIS_FORMAT = idum[28];
	CHECKPTREAD = idum[29];
	CHECKPTWRITE = idum[30];

	FDORDER = idum[31];
	MAXRELERROR = idum[32];
	RUN_MULTIPLE_SHOTS = idum[33];
	WRITE_MODELFILES = idum[34];
	
	ABS_TYPE = idum[35];
    
        FDORDER_TIME = idum[36];
    
        /*******************************************************/
        /* modified by Tao, in order to insert different source 
            * for FWI, RTM or some other cases */
        /*******************************************************/
        
        NSHOT1 = idum[37];
        NSHOT2 = idum[38];
        NSHOT_INC = idum[39];
        SRC_SNAP_TYPE = idum[40];
        ADD_REPLACE_VX = idum[41];
        ADD_REPLACE_VY = idum[42];
        ADD_REPLACE_SP = idum[43];
        ADD_REPLACE_SNAP = idum[44];
        FLAG_SU_BIN = idum[45];
        SAVE_LAST_SNAP = idum[46];
        READ_LAST_SNAP = idum[47];
        
        FORWARD_ONLY = idum[48];
        INV_VP = idum[49];
        INV_VS = idum[50];
        INV_RHO = idum[51];
        EPRECOND = idum[52];
        EPRECOND_PER_SHOT = idum[53];
        
        OMEGA1 = idum[54];
        OMEGA2 = idum[55];
        OMEGA_INC = idum[56];
        SEISMO2 = idum[57];
        
	MPI_Bcast(&FL[1],L,MPI_FLOAT,0,MPI_COMM_WORLD);

}
