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
 *   program SOFI2D, reading input-parameters from input-file or stdin
 *   that are formatted according to the json standard
 *
 *  ----------------------------------------------------------------------*/

#include <unistd.h>
#include "fd.h"

char ** varname_list,** value_list;

void read_par_json(FILE *fp, char *fileinp){

	/* declaration of extern variables */
	extern int   NX, NY, FDORDER, FDORDER_TIME, MAXRELERROR, SOURCE_TYPE, SOURCE_SHAPE, SNAP, SNAP_FORMAT, L;
	extern int SEISMO, NDT, SEIS_FORMAT, FREE_SURF, READMOD, READREC, SRCREC;
	extern int BOUNDARY, REC_ARRAY, DRX, LOG,  WRITE_MODELFILES; //RSG
	extern int  NPROCX, NPROCY, MYID, IDX, IDY, CHECKPTREAD, CHECKPTWRITE, RUN_MULTIPLE_SHOTS, ABS_TYPE, FW;
	extern int OUTNTIMESTEPINFO;
	/*extern int NGEOPH;*/
	extern float NGEOPH;
	extern float DH, TIME, DT, TS, *FL, TAU, DAMPING, PLANE_WAVE_DEPTH, PLANE_WAVE_ANGLE, FPML, VPPML,NPOWER, K_MAX_CPML;
	extern float XREC1, XREC2, YREC1, YREC2;
	extern float REC_ARRAY_DEPTH, REC_ARRAY_DIST;
	extern float TSNAP1, TSNAP2, TSNAPINC, REFREC[4];
	extern char  MFILE[STRING_SIZE], SIGNAL_FILE[STRING_SIZE], LOG_FILE[STRING_SIZE];
	extern char SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE], REC_FILE[STRING_SIZE];
	extern char SEIS_FILE[STRING_SIZE], CHECKPTFILE[STRING_SIZE];
        
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

	int number_readobjects=0,fserr=0;
	char errormessage[STRING_SIZE];

	char ** varname_list, ** value_list;
        int i;

	if (MYID == 0){

		//allocate first object in list
		varname_list = malloc(STRING_SIZE*sizeof(char*));
		value_list = malloc(STRING_SIZE*sizeof(char*));

		//read in objects from file
		number_readobjects=read_objects_from_intputfile(fp, fileinp, varname_list, value_list);
		fprintf(fp,"\nFrom input file %s, %i objects have been read in. \n",fileinp, number_readobjects);

		//print objects to screen
		fprintf(fp, "\n===========================================================");
		fprintf(fp, "\n=   List of Parameters read by the built in Json Parser   =");
		print_objectlist_screen(fp, number_readobjects, varname_list, value_list);

		//extract variables form object list

		/*=================================
		section general grid and discretization parameters
		  =================================*/
		if (get_int_from_objectlist("NPROCX",number_readobjects,&NPROCX,varname_list, value_list))
			declare_error("Variable NPROCX could not be retrieved from the json input file!");
		if (get_int_from_objectlist("NPROCY",number_readobjects,&NPROCY,varname_list, value_list))
			declare_error("Variable NPROCY could not be retrieved from the json input file!");
	/*	if (get_int_from_objectlist("RSG",number_readobjects,&RSG,varname_list, value_list))
			declare_error("Variable RSG could not be retrieved from the json input file!");*/
		if (get_int_from_objectlist("FDORDER",number_readobjects,&FDORDER,varname_list, value_list))
			declare_error("Variable FDORDER could not be retrieved from the json input file!");
        if (get_int_from_objectlist("FDORDER_TIME",number_readobjects,&FDORDER_TIME,varname_list, value_list)) {
            FDORDER_TIME=2;
        } else {
            if(FDORDER_TIME!=2 && FDORDER_TIME!=4) {
                declare_error("Only FDORDER_TIME 2 or 4 are supported!");
            }
        }
		if (get_int_from_objectlist("MAXRELERROR",number_readobjects,&MAXRELERROR,varname_list, value_list))
			declare_error("Variable MAXRELERROR could not be retrieved from the json input file!");
		if (get_int_from_objectlist("NX",number_readobjects,&NX,varname_list, value_list))
			declare_error("Variable NX could not be retrieved from the json input file!");
		if (get_int_from_objectlist("NY",number_readobjects,&NY,varname_list, value_list))
			declare_error("Variable NY could not be retrieved from the json input file!");
		if (get_float_from_objectlist("DH",number_readobjects,&DH,varname_list, value_list))
			declare_error("Variable DH could not be retrieved from the json input file!");
		if (get_float_from_objectlist("TIME",number_readobjects,&TIME,varname_list, value_list))
			declare_error("Variable TIME could not be retrieved from the json input file!");
		if (get_float_from_objectlist("DT",number_readobjects,&DT,varname_list, value_list))
			declare_error("Variable DT could not be retrieved from the json input file!");

		/*=================================
		 	 section source parameters
		  =================================*/

// 		if (get_int_from_objectlist("SOURCE_TYPE",number_readobjects,&SOURCE_TYPE,varname_list, value_list))
// 			declare_error("Variable SOURCE_TYPE could not be retrieved from the json input file!");
// 		if (get_int_from_objectlist("SOURCE_SHAPE",number_readobjects,&SOURCE_SHAPE,varname_list, value_list))
// 			declare_error("Variable SOURCE_SHAPE could not be retrieved from the json input file!");
// 		else {
// 			if (SOURCE_SHAPE==3) {
// 				if (get_string_from_objectlist("SIGNAL_FILE",number_readobjects,SIGNAL_FILE,varname_list, value_list))
// 					declare_error("Variable SIGNAL_FILE could not be retrieved from the json input file!");
// 			}
// 		}
// 		if (get_int_from_objectlist("SRCREC",number_readobjects,&SRCREC,varname_list, value_list))
// 			declare_error("Variable SRCREC could not be retrieved from the json input file!");
// 		else {
// 			if (SRCREC==1) {
// 				if (get_string_from_objectlist("SOURCE_FILE",number_readobjects,SOURCE_FILE,varname_list, value_list))
// 					declare_error("Variable SOURCE_FILE could not be retrieved from the json input file!");
// 				if (get_int_from_objectlist("RUN_MULTIPLE_SHOTS",number_readobjects,&RUN_MULTIPLE_SHOTS,varname_list, value_list))
// 					declare_error("Variable RUN_MULTIPLE_SHOTS could not be retrieved from the json input file!");
// 			}
// 			if (SRCREC==2) {
// 				if (get_float_from_objectlist("PLANE_WAVE_DEPTH",number_readobjects,&PLANE_WAVE_DEPTH,varname_list, value_list))
// 					declare_error("Variable PLANE_WAVE_DEPTH could not be retrieved from the json input file!");
// 				else {
// 					if (PLANE_WAVE_DEPTH>0.0) {
// 						if (get_float_from_objectlist("PLANE_WAVE_ANGLE",number_readobjects,&PLANE_WAVE_ANGLE,varname_list, value_list))
// 							declare_error("Variable PLANE_WAVE_ANGLE could not be retrieved from the json input file!");
// 						if (get_float_from_objectlist("TS",number_readobjects,&TS,varname_list, value_list))
// 							declare_error("Variable TS could not be retrieved from the json input file!");
// 					}
// 				}
// 			}
// 		}


		/*=================================
		 	 section boundary parameters
		  =================================*/

		if (get_int_from_objectlist("FREE_SURF",number_readobjects,&FREE_SURF,varname_list, value_list))
			declare_error("Variable FREE_SURF could not be retrieved from the json input file!");
		if (get_int_from_objectlist("BOUNDARY",number_readobjects,&BOUNDARY,varname_list, value_list))
			declare_error("Variable BOUNDARY could not be retrieved from the json input file!");
		if (get_int_from_objectlist("FW",number_readobjects,&FW,varname_list, value_list))
			declare_error("Variable FW could not be retrieved from the json input file!");
		if (get_int_from_objectlist("ABS_TYPE",number_readobjects,&ABS_TYPE,varname_list, value_list))
                	declare_error("Variable ABS_TYPE could not be retrieved from the json input file!");

		if (ABS_TYPE==1) {
			if (get_float_from_objectlist("NPOWER",number_readobjects,&NPOWER,varname_list, value_list)) 
				declare_error("Variable NPOWER could not be retrieved from the json input file!");
			if (get_float_from_objectlist("K_MAX_CPML",number_readobjects,&K_MAX_CPML,varname_list, value_list)) 
				declare_error("Variable K_MAX_CPML could not be retrieved from the json input file!");
			if (get_float_from_objectlist("FPML",number_readobjects,&FPML,varname_list, value_list))
				declare_error("Variable FPML could not be retrieved from the json input file!");
			if (get_float_from_objectlist("VPPML",number_readobjects,&VPPML,varname_list, value_list))
				declare_error("Variable VPPML could not be retrieved from the json input file!");
                }
		if (ABS_TYPE==2) {
			if (get_float_from_objectlist("DAMPING",number_readobjects,&DAMPING,varname_list, value_list))
                        	declare_error("Variable DAMPING could not be retrieved from the json input file!");
			}


		/*=================================
			 section snapshot parameters
		=================================*/
		if (get_int_from_objectlist("SNAP",number_readobjects,&SNAP,varname_list, value_list))
			declare_error("Variable SNAP could not be retrieved from the json input file!");
		else {
			if (SNAP>0) {
				if (get_int_from_objectlist("SNAP_FORMAT",number_readobjects,&SNAP_FORMAT,varname_list, value_list))
					declare_error("Variable SNAP_FORMAT could not be retrieved from the json input file!");
				if (get_float_from_objectlist("TSNAP1",number_readobjects,&TSNAP1,varname_list, value_list))
					declare_error("Variable TSNAP1 could not be retrieved from the json input file!");
				if (get_float_from_objectlist("TSNAP2",number_readobjects,&TSNAP2,varname_list, value_list))
					declare_error("Variable TSNAP2 could not be retrieved from the json input file!");
				if (get_float_from_objectlist("TSNAPINC",number_readobjects,&TSNAPINC,varname_list, value_list))
					declare_error("Variable TSNAPINC could not be retrieved from the json input file!");
				if (get_string_from_objectlist("SNAP_FILE",number_readobjects,SNAP_FILE,varname_list, value_list))
					declare_error("Variable SNAP_FILE could not be retrieved from the json input file!");
			}
			if (SNAP==6){
                            TSNAP1=DT;
                            TSNAP2=TIME;
                            fprintf(fp, "\nVariable TSNAP1 is set to default value for freq.-domain wavefields: %f\n",TSNAP1);
                            fprintf(fp, "\nVariable TSNAP2 is set to default value for freq.-domain wavefields: %f\n",TSNAP2);
                        }
		}
		/* increments are read in any case, because they will be also used as increment for model output */
		if (get_int_from_objectlist("IDX",number_readobjects,&IDX,varname_list, value_list))
			declare_error("Variable IDX could not be retrieved from the json input file!");
		if (get_int_from_objectlist("IDY",number_readobjects,&IDY,varname_list, value_list))
			declare_error("Variable IDY could not be retrieved from the json input file!");
                
                if(SNAP==6){
                    if (get_string_from_objectlist("FREQ_FILE",number_readobjects,FREQ_FILE,varname_list, value_list)){
                        declare_error("\nVariable FREQ_FILE didn't set in the json input file!\n");
                    }
                }else if(SNAP==7){
                    if (get_int_from_objectlist("OMEGA1",number_readobjects,&OMEGA1,varname_list, value_list))
                        declare_error("Variable OMEGA1 could not be retrieved from the json input file!");
                    if (get_int_from_objectlist("OMEGA2",number_readobjects,&OMEGA2,varname_list, value_list))
                        declare_error("Variable OMEGA2 could not be retrieved from the json input file!");
                    if (get_int_from_objectlist("OMEGA_INC",number_readobjects,&OMEGA_INC,varname_list, value_list))
                        declare_error("Variable OMEGA_INC could not be retrieved from the json input file!");
                    if((OMEGA2<OMEGA1)||OMEGA1<1||OMEGA_INC<1) declare_error("Reset OMEGA1, OMEGA2, OMEGA_INC!");
                }

		/*=================================
			section seismogramm parameters
		=================================*/
                if (get_int_from_objectlist("SEISMO2",number_readobjects,&SEISMO2,varname_list, value_list))
                    SEISMO2=0;
		if (get_int_from_objectlist("SEISMO",number_readobjects,&SEISMO,varname_list, value_list))
			declare_error("Variable SEISMO could not be retrieved from the json input file!");
		else {
			if (SEISMO>0){
				if (get_string_from_objectlist("REC_FILE",number_readobjects,REC_FILE,varname_list, value_list))
					declare_error("Variable REC_FILE could not be retrieved from the json input file!");
				if (get_string_from_objectlist("SEIS_FILE",number_readobjects,SEIS_FILE,varname_list, value_list))
					declare_error("Variable SEIS_FILE could not be retrieved from the json input file!");
				if (get_int_from_objectlist("READREC",number_readobjects,&READREC,varname_list, value_list))
					declare_error("Variable READREC could not be retrieved from the json input file!");
				else {
					if (READREC==0) {
						if (get_float_from_objectlist("XREC1",number_readobjects,&XREC1,varname_list, value_list))
							declare_error("Variable XREC1 could not be retrieved from the json input file!");
						if (get_float_from_objectlist("XREC2",number_readobjects,&XREC2,varname_list, value_list))
							declare_error("Variable XREC2T could not be retrieved from the json input file!");
						if (get_float_from_objectlist("YREC1",number_readobjects,&YREC1,varname_list, value_list))
							declare_error("Variable YREC1 could not be retrieved from the json input file!");
						if (get_float_from_objectlist("YREC2",number_readobjects,&YREC2,varname_list, value_list))
							declare_error("Variable YREC2 could not be retrieved from the json input file!");
					}
				}
				if (get_int_from_objectlist("NDT",number_readobjects,&NDT,varname_list, value_list))
					declare_error("Variable NDT could not be retrieved from the json input file!");
				if (get_int_from_objectlist("SEIS_FORMAT",number_readobjects,&SEIS_FORMAT,varname_list, value_list))
					declare_error("Variable SEIS_FORMAT could not be retrieved from the json input file!");

				if (get_int_from_objectlist("REC_ARRAY",number_readobjects,&REC_ARRAY,varname_list, value_list))
					declare_error("Variable REC_ARRAY could not be retrieved from the json input file!");
				else {
					if (REC_ARRAY>0) {
						if (get_int_from_objectlist("DRX",number_readobjects,&DRX,varname_list, value_list))
							declare_error("Variable DRX could not be retrieved from the json input file!");
						if (get_float_from_objectlist("REC_ARRAY_DEPTH",number_readobjects,&REC_ARRAY_DEPTH,varname_list, value_list))
							declare_error("Variable REC_ARRAY_DEPTH could not be retrieved from the json input file!");
						if (get_float_from_objectlist("REC_ARRAY_DIST",number_readobjects,&REC_ARRAY_DIST,varname_list, value_list))
							declare_error("Variable REC_ARRAY_DIST could not be retrieved from the json input file!");
					}
				}
				if (get_float_from_objectlist("REFRECX",number_readobjects,&REFREC[1],varname_list, value_list))
					declare_error("Variable REFRECX could not be retrieved from the json input file!");
				if (get_float_from_objectlist("REFRECY",number_readobjects,&REFREC[2],varname_list, value_list))
					declare_error("Variable REFRECY could not be retrieved from the json input file!");
				if (get_float_from_objectlist("NGEOPH",number_readobjects,&NGEOPH,varname_list, value_list))
					declare_error("Variable NGEOPH could not be retrieved from the json input file!");
			}
		}

		/*=================================
			section general model and log parameters
		  =================================*/
		if (get_string_from_objectlist("MFILE",number_readobjects,MFILE,varname_list, value_list))
			declare_error("Variable MFILE could not be retrieved from the json input file!");
		(get_int_from_objectlist("WRITE_MODELFILES",number_readobjects,&WRITE_MODELFILES,varname_list, value_list));
		if (get_int_from_objectlist("LOG",number_readobjects,&LOG,varname_list, value_list))
			declare_error("Variable LOG could not be retrieved from the json input file!");
		if (get_int_from_objectlist("CHECKPTREAD",number_readobjects,&CHECKPTREAD,varname_list, value_list))
			declare_error("Variable CHECKPTREAD could not be retrieved from the json input file!");
		if (get_int_from_objectlist("CHECKPTWRITE",number_readobjects,&CHECKPTWRITE,varname_list, value_list))
			declare_error("Variable CHECKPTWRITE could not be retrieved from the json input file!");
		if (get_string_from_objectlist("LOG_FILE",number_readobjects,LOG_FILE,varname_list, value_list))
			declare_error("Variable LOG_FILE could not be retrieved from the json input file!");
		if (get_string_from_objectlist("CHECKPT_FILE",number_readobjects,CHECKPTFILE,varname_list, value_list))
			declare_error("Variable CHECKPT_FILE could not be retrieved from the json input file!");
		if (get_int_from_objectlist("READMOD",number_readobjects,&READMOD,varname_list, value_list))
			declare_error("Variable READMOD could not be retrieved from the json input file!");
		if (get_int_from_objectlist("OUT_TIMESTEP_INFO",number_readobjects,&OUTNTIMESTEPINFO,varname_list, value_list)) {}

		if (get_float_from_objectlist("TAU",number_readobjects,&TAU,varname_list, value_list))
			declare_error("Variable TAU could not be retrieved from the json input file!");
		if (get_int_from_objectlist("L",number_readobjects,&L,varname_list, value_list))
			declare_error("Variable L could not be retrieved from the json input file!");
        else {
            FL=vector(1,L);
            switch(L) {
			case 0:
				break;
			case 5:
				if (get_float_from_objectlist("FL5",number_readobjects,&FL[5],varname_list, value_list))
					declare_error("Variable FL5 could not be retrieved from the json input file!");
			case 4:
				if (get_float_from_objectlist("FL4",number_readobjects,&FL[4],varname_list, value_list))
					declare_error("Variable FL4 could not be retrieved from the json input file!");
			case 3:
				if (get_float_from_objectlist("FL3",number_readobjects,&FL[3],varname_list, value_list))
					declare_error("Variable FL3 could not be retrieved from the json input file!");
			case 2:
				if (get_float_from_objectlist("FL2",number_readobjects,&FL[2],varname_list, value_list))
					declare_error("Variable FL2 could not be retrieved from the json input file!");
			case 1:
				if (get_float_from_objectlist("FL1",number_readobjects,&FL[1],varname_list, value_list))
					declare_error("Variable FL1 could not be retrieved from the json input file!");
				break;
			default:
				declare_error("More than four relaxation Parameter (L>5) are not implemented yet!");
				break;
			}
		}

		/********************************************/
		/* Check files and directories if necessary */
		/********************************************/

// 		/* signal file */
// 		if (SOURCE_SHAPE == 3)
// 		{
// 			if (access(SIGNAL_FILE,0) != 0)
// 			{
// 				fprintf(fp, "\n==================================================================\n");
// 				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
// 				fprintf(fp, "        The signal file does not exist!\n");
// 				fprintf(fp, "        File name: <%s>", SIGNAL_FILE);
// 				fprintf(fp, "\n==================================================================\n");
// 				fserr = 1;
// 			}
// 			else if (access(SIGNAL_FILE,4) != 0)
// 			{
// 				fprintf(fp, "\n==================================================================\n");
// 				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
// 				fprintf(fp, "        The signal file does not have read access!\n");
// 				fprintf(fp, "        File name: <%s>", SIGNAL_FILE);
// 				fprintf(fp, "\n==================================================================\n");
// 				fserr = 1;
// 			}
// 		}
// 
// 		/* source file */
// 		if (SRCREC==1)
// 		{
// 			if (access(SOURCE_FILE,0) != 0)
// 			{
// 				fprintf(fp, "\n==================================================================\n");
// 				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
// 				fprintf(fp, "        The source file does not exist!\n");
// 				fprintf(fp, "        File name: <%s>", SOURCE_FILE);
// 				fprintf(fp, "\n==================================================================\n");
// 				fserr = 1;
// 			}
// 			else if (access(SOURCE_FILE,4) != 0)
// 			{
// 				fprintf(fp, "\n==================================================================\n");
// 				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
// 				fprintf(fp, "        The source file does not have read access!\n");
// 				fprintf(fp, "        File name: <%s>", SOURCE_FILE);
// 				fprintf(fp, "\n==================================================================\n");
// 				fserr = 1;
// 			}
// 		}

		/* model file -> this is checked in more detail in readmod.c!
		if (READMOD)
		{
			if (access(MFILE,0) != 0)
			{
				fprintf(fp, "\n==================================================================\n");
				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
				fprintf(fp, "        The model file does not exist!\n");
				fprintf(fp, "        File name: <%s>", MFILE);
				fprintf(fp, "\n==================================================================\n");
				fserr = 1;
			}
			else if (access(MFILE,4) != 0)
			{
				fprintf(fp, "\n==================================================================\n");
				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
				fprintf(fp, "        The model file does not have read access!\n");
				fprintf(fp, "        File name: <%s>", MFILE);
				fprintf(fp, "\n==================================================================\n");
				fserr = 1;
			}
		}*/

		/* receiver file */
		if (READREC)
		{
			if (access(REC_FILE,0) != 0)
			{
				fprintf(fp, "\n==================================================================\n");
				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
				fprintf(fp, "        The receiver file does not exist!\n");
				fprintf(fp, "        File name: <%s>", REC_FILE);
				fprintf(fp, "\n==================================================================\n");
				fserr = 1;
			}
			else if (access(REC_FILE,4) != 0)
			{
				fprintf(fp, "\n==================================================================\n");
				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
				fprintf(fp, "        The receiver file does not have read access!\n");
				fprintf(fp, "        File name: <%s>", REC_FILE);
				fprintf(fp, "\n==================================================================\n");
				fserr = 1;
			}
		}

		/* checkpoint file */
		if (CHECKPTREAD || CHECKPTWRITE)
		{
			if (access(CHECKPTFILE,0) != 0)
			{
				fprintf(fp, "\n==================================================================\n");
				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
				fprintf(fp, "        The checkpoint file does not exist!\n");
				fprintf(fp, "        File name: <%s>", CHECKPTFILE);
				fprintf(fp, "\n==================================================================\n");
				fserr = 1;
			}
			else if (access(CHECKPTFILE,6) != 0)
			{
				fprintf(fp, "\n==================================================================\n");
				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
				fprintf(fp, "        The checkpoint file does not have read and/or write access!\n");
				fprintf(fp, "        File name: <%s>", CHECKPTFILE);
				fprintf(fp, "\n==================================================================\n");
				fserr = 1;
			}
		}


		/********************************************/
		/* ERROR                                    */
		/********************************************/
		if (fserr)
		{
			fprintf(fp, "\n");
			sprintf(errormessage, "\n  in: <read_par_json.c> \n");
			declare_error(errormessage);
		}
		
                /*******************************************************/
                /* modified by Tao, in order to insert different source 
                    * for FWI, RTM or some other cases */
                /*******************************************************/
                /*set the way of adding source functions */
                if (get_int_from_objectlist("FLAG_SU_BIN",number_readobjects,&FLAG_SU_BIN,varname_list, value_list))
                    FLAG_SU_BIN=1;
                
                /* some parameters with respect to the boundary saving strategy */
                if(FLAG_SU_BIN==1||FLAG_SU_BIN==3){
                    i=0;
                    if (get_string_from_objectlist("SRC_F_SP",number_readobjects,SRC_F_SP,varname_list, value_list)){
                        fprintf(fp, "\nVariable SRC_F_SP didn't set in the json input file!\n");
                        i=i+1;
                    }
                    /*if (get_string_from_objectlist("SRC_F_VX",number_readobjects,SRC_F_VX,varname_list, value_list)){
                        fprintf(fp, "\nVariable SRC_F_VX didn't set in the json input file!\n");
                        i=i+1;
                    }
                    if (get_string_from_objectlist("SRC_F_VY",number_readobjects,SRC_F_VY,varname_list, value_list)){
                        fprintf(fp, "\nVariable SRC_F_VY didn't set in the json input file!\n");
                        i=i+1;
                    }
                    if(i==3) declare_error("Source file could not be retrieved from the json input file when FLAG_SU_BIN is equal to 1!");
                    if (get_int_from_objectlist("ADD_REPLACE_VX",number_readobjects,&ADD_REPLACE_VX,varname_list, value_list))
                        ADD_REPLACE_VX=0;
                    if (get_int_from_objectlist("ADD_REPLACE_VY",number_readobjects,&ADD_REPLACE_VY,varname_list, value_list))
                        ADD_REPLACE_VY=0;
                    if (get_int_from_objectlist("ADD_REPLACE_SP",number_readobjects,&ADD_REPLACE_SP,varname_list, value_list))
                        ADD_REPLACE_SP=0;*/  
                }
                
//                 if(FLAG_SU_BIN==2||FLAG_SU_BIN==3){
//                     if (get_int_from_objectlist("SRC_SNAP_TYPE",number_readobjects,&SRC_SNAP_TYPE,varname_list, value_list))
//                         SRC_SNAP_TYPE=1;
//                     if(SRC_SNAP_TYPE==2){
//                         if (get_float_from_objectlist("MEMORY_DEPTH",number_readobjects,&MEMORY_DEPTH,varname_list, value_list))
//                             declare_error("Variable MEMORY_DEPTH could not be retrieved from the json input file!");
//                     }else if(SRC_SNAP_TYPE==3){
//                         if (get_string_from_objectlist("SNAP_COORDS_FILE",number_readobjects,SNAP_COORDS_FILE,varname_list, value_list))
//                             declare_error("Variable SNAP_COORDS_FILE could not be retrieved from the json input file!");
//                     }
//                     if (get_string_from_objectlist("SRC_F_SNAP",number_readobjects,SRC_F_SNAP,varname_list, value_list))
//                             declare_error("Variable SRC_F_SNAP could not be retrieved from the json input file!");
//                     
//                     if (get_int_from_objectlist("ADD_REPLACE_SNAP",number_readobjects,&ADD_REPLACE_SNAP,varname_list, value_list))
//                         ADD_REPLACE_SNAP=3;
//                     
//                     if (get_int_from_objectlist("READ_LAST_SNAP",number_readobjects,&READ_LAST_SNAP,varname_list, value_list))
//                         READ_LAST_SNAP=0;
//                     
//                 }
                
//                 if(FLAG_SU_BIN==1&&SNAP==5){
//                     if (get_int_from_objectlist("SRC_SNAP_TYPE",number_readobjects,&SRC_SNAP_TYPE,varname_list, value_list))
//                         SRC_SNAP_TYPE=1;
//                     if(SRC_SNAP_TYPE==2){
//                         if (get_float_from_objectlist("MEMORY_DEPTH",number_readobjects,&MEMORY_DEPTH,varname_list, value_list))
//                             declare_error("Variable MEMORY_DEPTH could not be retrieved from the json input file!");
//                     }else if(SRC_SNAP_TYPE==3){
//                         if (get_string_from_objectlist("SNAP_COORDS_FILE",number_readobjects,SNAP_COORDS_FILE,varname_list, value_list))
//                             declare_error("Variable SNAP_COORDS_FILE could not be retrieved from the json input file!");
//                     }
//                     
//                     if (get_int_from_objectlist("ADD_REPLACE_SNAP",number_readobjects,&ADD_REPLACE_SNAP,varname_list, value_list))
//                         ADD_REPLACE_SNAP=3;
//                     if (get_int_from_objectlist("SAVE_LAST_SNAP",number_readobjects,&SAVE_LAST_SNAP,varname_list, value_list))
//                         SAVE_LAST_SNAP=0;
//                 }
                
//                 if(READREC==2){
//                     if (get_float_from_objectlist("SHOT_SPACING_X",number_readobjects,&SHOT_SPACING_X,varname_list, value_list))
//                         declare_error("Variable SHOT_SPACING_X could not be retrieved from the json input file!");
//                     if (get_float_from_objectlist("SHOT_SPACING_Y",number_readobjects,&SHOT_SPACING_Y,varname_list, value_list)){
//                         SHOT_SPACING_Y=0.0;
//                         fprintf(fp, "\nVariable SHOT_SPACING_Y is set to default value %f\n",SHOT_SPACING_Y);}
//                 }
                
                if (get_int_from_objectlist("FORWARD_ONLY",number_readobjects,&FORWARD_ONLY,varname_list, value_list))
                    FORWARD_ONLY=1;
                //if(FORWARD_ONLY!=1){
                    if (get_int_from_objectlist("INV_VP",number_readobjects,&INV_VP,varname_list, value_list))
                        INV_VP=0;
                    if (get_int_from_objectlist("INV_VS",number_readobjects,&INV_VS,varname_list, value_list))
                        INV_VS=0;
                    if (get_int_from_objectlist("INV_RHO",number_readobjects,&INV_RHO,varname_list, value_list))
                        INV_RHO=0;
                    if (get_float_from_objectlist("TSHIFT",number_readobjects,&TSHIFT,varname_list, value_list)){
                        TSHIFT = 0.0;
                        fprintf(fp, "\nVariable TSHIFT is set to default value %f\n",TSHIFT);
                    }
                    
                //}
                /* Maybe useful for illumination */
                if (get_int_from_objectlist("EPRECOND",number_readobjects,&EPRECOND,varname_list, value_list))
                    EPRECOND=0;
                if (get_int_from_objectlist("EPRECOND_PER_SHOT",number_readobjects,&EPRECOND_PER_SHOT,varname_list, value_list))
                    EPRECOND_PER_SHOT=0;
                if (get_float_from_objectlist("EPSILON",number_readobjects,&EPSILON,varname_list, value_list)){
                    EPSILON = 0.03;
                    fprintf(fp, "\nVariable EPSILON is set to default value %f\n",EPSILON);
                }
                
                
                
            
	} /* End of if(MYID==0) */
}
