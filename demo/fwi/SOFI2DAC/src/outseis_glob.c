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
 *   Write merged seismograms from all PEs to disk
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "segy.h"

void  outseis_glob(FILE *fp, FILE *fpdata, float **section,
		int **recpos, int **recpos_loc, int ntr, float Sx, 
                float Sy, int nsrc, int ns, int seis_form, int comp){

	/* declaration of extern variables */
	extern int NDT;
	extern float  DH, DT, REFREC[4];

	const float xshift=800.0, yshift=800.0;

	//fprintf(fp,"\n **Message from function outseis_glob \n\n");

	/* declaration of local variables */
	int i,j;
	segy tr;
	int itr ;
	float xr, yr, x, y;
// 	float Sx=0.0, Sy=0.0;

	/* if there are more than one source position/coordinate specified
	 * in SOURCE_FILE, only the first position is written into trace header fields,
	 * in case of RUN_MULTIPLE_SHOTS is activated the individual shot position is used */
// 	Sx=srcpos[1][ishot];
// 	Sy=srcpos[2][ishot];

	switch(seis_form){
	case 1 :
		for(itr=1;itr<=ntr;itr++){        /*SEGY (without file-header)*/
			xr=recpos[1][itr]*DH;
			yr=recpos[2][itr]*DH;
			x=xr-REFREC[1];
			y=yr-REFREC[2];

			/*tr.tracl=(int)recpos[3][itr];*/
                        tr.tracl=(int)itr;
                        tr.ep=0;     
                        tr.cdp=(int)0;
                        tr.trid=(short)1;           /* trace identification code: 1=seismic*/
                        tr.offset=(signed int)0;
                        tr.gelev=(signed int)iround(yr*1000.0);
                        tr.sdepth=(signed int)0;   /* source depth (positive) */
                        /* angle between receiver position and reference point
                            (sperical coordinate system: swdep=theta, gwdep=phi) */
                        tr.swdep=iround(1000.0);
                        /*tr.scalel=(signed short)-3;*/
                        tr.scalel=(signed short)-1000;
                        /*tr.scalco=(signed short)-3;*/
                        tr.scalco=(signed short)-1000;
                        tr.sx=(signed int)iround(Sx*1000.0);  /* X source coordinate */

                                /* group coordinates */
                        tr.gx=(signed int)iround(xr*1000.0);

                        tr.ns=(unsigned short)ns; /* number of samples in this trace */
                        tr.dt=(unsigned short)iround(((float)DT)*1.0e6); /* sample interval in micro-seconds */
                        /*tr.d1=(float)(TIME/ns2);*/        /* sample spacing for non-seismic data */
                        tr.d1=(float)tr.dt*1.0e-6;
                        
                        tr.tracr=0	;	/* trace sequence number within reel */

                        tr.fldr=0       ;	/* field record number */

                        tr.tracf=0      ;	/* trace number within field record */

                        tr.ep=0         ;	/* energy source point number */

                        tr.cdpt=0       ;	/* trace number within CDP ensemble */


                        tr.nvs=0       	;   /* number of vertically summed traces (see vscode
                                                in bhed structure) */

                        tr.nhs=0       	;   /* number of horizontally summed traces (see vscode
                                                in bhed structure) */

                        tr.duse=0     	;   /* data use:
                                                1 = production
                                                2 = test */

                        tr.selev=0     	; /* source elevation from sea level
                                                (above sea level is positive) */

                        tr.gdel=0      	; /* datum elevation at receiver group */

                        tr.sdel=0      	; /* datum elevation at source */

                        tr.gwdep=0     	; /* water depth at receiver group */

                        
                        tr.sy= (signed int)iround(Sy*1000.0)    ;   /* Y source coordinate */


                        tr.gy= (signed int)iround(yr*1000.0)    ;   /* Y group coordinate */

                        tr.counit=0    	;   /* coordinate units code:
                                                for previous four entries
                                                1 = length (meters or feet)
                                                2 = seconds of arc (in this case, the
                                                X values are longitude and the Y values
                                                are latitude, a positive value designates
                                                the number of seconds east of Greenwich
                                                or north of the equator */

                        tr.wevel=0     ;	/* weathering velocity */

                        tr.swevel=0    ;	/* subweathering velocity */

                        tr.sut=0       ;	/* uphole time at source */

                        tr.gut=0       ;	/* uphole time at receiver group */

                        tr.sstat=0     ;	/* source static correction */

                        tr.gstat=0     ;	/* group static correction */

                        tr.tstat=0     ;	/* total static applied */

                        tr.laga=0      ; /* lag time A, time in ms between end of 240-
                                                byte trace identification header and time
                                                break, positive if time break occurs after
                                                end of header, time break is defined as
                                                the initiation pulse which maybe recorded
                                                on an auxiliary trace or as otherwise
                                                specified by the recording system */

                        tr.lagb=0     	; /* lag time B, time in ms between the time break
                                                and the initiation time of the energy source,
                                                may be positive or negative */

                        tr.delrt=0     	; /* delay recording time, time in ms between
                                                initiation time of energy source and time
                                                when recording of data samples begins
                                                (for deep water work if recording does not
                                                start at zero time) */

                        tr.muts=0      ; /* mute time--start */

                        tr.mute=0      ; /* mute time--end */


                        tr.gain=0      ; /* gain type of field instruments code:
                                                1 = fixed
                                                2 = binary
                                                3 = floating point
                                                4 ---- N = optional use */

                        tr.igc=0       ; /* instrument gain constant */

                        tr.igi=0       ; /* instrument early or initial gain */

                        tr.corr=0      ; /* correlated:
                                                1 = no
                                                2 = yes */

                        tr.sfs=0       ; /* sweep frequency at start */

                        tr.sfe=0       ; /* sweep frequency at end */

                        tr.slen=0      ; /* sweep length in ms */

                        tr.styp=0      ; /* sweep type code:
                                                1 = linear
                                                2 = cos-squared
                                                3 = other */

                        tr.stas=0      	; /* sweep trace length at start in ms */

                        tr.stae=0      	; /* sweep trace length at end in ms */

                        tr.tatyp=0     	; /* taper type: 1=linear, 2=cos^2, 3=other */

                        tr.afilf=0     	; /* alias filter frequency if used */

                        tr.afils=0     	; /* alias filter slope */

                        tr.nofilf=0    	; /* notch filter frequency if used */

                        tr.nofils=0   	; /* notch filter slope */

                        tr.lcf=0      	; /* low cut frequency if used */

                        tr.hcf=0      	; /* high cut frequncy if used */

                        tr.lcs=0       	; /* low cut slope */

                        tr.hcs=0       	; /* high cut slope */

                        tr.year=0      	; /* year data recorded */

                        tr.day=0       	; /* day of year */

                        tr.hour=0     	; /* hour of day (24 hour clock) */

                        tr.minute=0    	; /* minute of hour */

                        tr.sec=0       	; /* second of minute */

                        tr.timbas=0    	; /* time basis code:
                                                1 = local
                                                2 = GMT
                                                3 = other */

                        tr.trwf=0      	; /* trace weighting factor, defined as 1/2^N
                                                volts for the least sigificant bit */

                        tr.grnors=0   	; /* geophone group number of roll switch
                                                position one */

                        tr.grnofr=0    	; /* geophone group number of trace one within
                                                original field record */

                        tr.grnlof=0    	; /* geophone group number of last trace within
                                                original field record */

                        tr.gaps=0      	;  /* gap size (total number of groups dropped) */

                        tr.otrav=0     	;  /* overtravel taper code:
                                                1 = down (or behind)
                                                2 = up (or ahead) */

                        /* local assignments */

                        tr.f1=0.0;	/* first sample location for non-seismic data */

                        tr.d2=0.0;	/* sample spacing between traces */

                        tr.f2=0.0;	/* first trace location */

                        tr.ungpow=0.0;	/* negative of power used for dynamic
                                                range compression */

                        tr.unscale=0.0;	/* reciprocal of scaling factor to normalize
                                                range */
                        tr.ntr=0      ;   /* number of traces */

                        tr.mark=0     ;


			for(j=1;j<=ns;j++) tr.data[j]=section[itr][j];

			fwrite(&tr,240,1,fpdata);
			fwrite(&tr.data[1],4,ns,fpdata);


		}
		break;


	case 2 :
		for(i=1;i<=ntr;i++){         /*ASCII ONE COLUMN*/
			for(j=1;j<=ns;j++) fprintf(fpdata,"%e\n", section[i][j]);
		}
		break;

	case 3 :                             /*BINARY */

		for(i=1;i<=ntr;i++)
			for(j=1;j<=ns;j++){
				fwrite(&section[i][j],sizeof(float),1,fpdata); }
		break;

	default :
		fprintf(fp," Don't know data format for seismograms !\n");
		fprintf(fp," No output written. ");
	}

	fclose(fpdata);
}
