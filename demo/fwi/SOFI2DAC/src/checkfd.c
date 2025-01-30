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

/*-------------------------------------------------------------
 *  Check FD-Grid for stability and grid dispersion.
 *  If the stability criterion is not fullfilled the program will
 *  terminate.                   
 *
 *  ----------------------------------------------------------*/

#include <limits.h>
#include "fd.h"

void checkfd(FILE *fp, float ** prho, float ** ppi, float ** ptaup, float *peta, float *hc){

	extern float DH, DT, TS;
	extern float XREC1, XREC2, YREC1, YREC2;
//         extern int NX, NY, MYID, PARAMETERIZATION, FW, L, NT, NDT, ACOUSTIC;
        extern int NX, NY, MYID, FW, L, NT, NDT;
	extern int READREC, NPROCX,NPROCY, SRCREC, FREE_SURF;
        extern int SEISMO, SEIS_FORMAT[6];
	extern char SOURCE_FILE[STRING_SIZE], REC_FILE[STRING_SIZE];

	/* local variables */

	float  c = 0.0, cmax_p=0.0, cmin_p=1e9, cmax_s=0.0, cmin_s=1e9, fmax, gamma;
	float  cmax=0.0, cmin=1e9, dtstab, dhstab, cmax_r, cmin_r;
	float sumu, sumpi, ws, ts, ppi_ref = 0.0, pu_ref;
	int nfw=iround(FW/DH);
	int i, j, k, l, ny1=1, nx, ny, nx_min, ny_min;
	float srec_minx=DH*NX*NPROCX+1, srec_miny=DH*NY*NPROCY+1;
	float srec_maxx=-1.0, srec_maxy=-1.0;
	
	
	ws=2.0*PI/TS; /*center frequency of source*/
	
	nx=NX; ny=NY; 

	/* low Q frame not yet applied as a absorbing boundary */
	/* if (!FREE_SURF) ny1=1+nfw;*/
	nfw=0;
	
	/* find maximum model phase velocity of P-waves at infinite
		 frequency within the whole model */
	if (L){  /*viscoelastic*/
		for (i=1+nfw;i<=(nx-nfw);i++){
			for (j=ny1;j<=(ny-nfw);j++){
			        
				
				sumpi=0.0;
				for (l=1;l<=L;l++){
					ts=DT/peta[l];
					sumpi=sumpi+((ws*ws*ts*ts*ptaup[j][i])/(1.0+ws*ws*ts*ts));
				}
				
// 				if(PARAMETERIZATION==1){
// 					ppi_ref=prho[j][i]*ppi[j][i]*ppi[j][i];}
// 				if(PARAMETERIZATION==3){
// 					ppi_ref=ppi[j][i]+2*pu[j][i];
                                        ppi_ref=ppi[j][i];
//                                 }
					
				/* minimum phase velocity of P waves */
				c=sqrt(ppi_ref/(prho[j][i]*(1.0+sumpi)));
				
				if (cmin_p>c) cmin_p=c;
				
				
				/* maximum phase velocity of shear waves */
				c=sqrt(ppi_ref*(1.0+L*ptaup[j][i])/(prho[j][i]*(1.0+sumpi)));
				
				if (cmax_p<c) cmax_p=c;
			}
		}
	}
	else{ /*L=0, elastic*/
		for (i=1+nfw;i<=(nx-nfw);i++){
			for (j=ny1;j<=(ny-nfw);j++){

// 				if(PARAMETERIZATION==3){
// 					c=sqrt((ppi[j][i]+2.0*pu[j][i])/prho[j][i]);
                                        c=sqrt(ppi[j][i]/prho[j][i]);
//                                 }

// 				if(PARAMETERIZATION==1){
// 					c=ppi[j][i];}

				if (cmax_p<c) cmax_p=c;
				if (cmin_p>c) cmin_p=c;
			}
		}

	}


	if (MYID==0){
		fprintf(fp,"\n\n\n **Message from checkfd (printed by PE %d):\n",MYID);
		fprintf(fp," Minimum and maximum P-wave velocities within subvolumes: \n ");
		fprintf(fp," MYID\t Vp_min(f=fc) \t Vp_max(f=inf) \n");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	fprintf(fp," %d \t %e \t %e \n", MYID, cmin_p, cmax_p);


        cmax=cmax_p;
        cmin=cmin_p;
	

	/* find global maximum for Vp and global minimum for Vs*/
	MPI_Allreduce(&cmax,&cmax_r,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	MPI_Allreduce(&cmin,&cmin_r,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);
	cmax=cmax_r;
	cmin=cmin_r;	

	fmax=2.0/TS;
	dhstab = (cmin/(hc[0]*fmax));
	gamma = fabs(hc[1]) + fabs(hc[2]) + fabs(hc[3]) + fabs(hc[4]) + fabs(hc[5]) + fabs(hc[6]);
	dtstab = DH/(sqrt(2)*gamma*cmax);
	/*dtstab=DH/(sqrt(2.0)*cmax);*/

	/* find global minimum for NX and NY */
	MPI_Allreduce(&NX,&nx_min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
	MPI_Allreduce(&NY,&ny_min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
	

	if (MYID == 0) {

	fprintf(fp," Global values for entire model: \n");

        fprintf(fp," Vp_max= %e m/s \t Vp_min=%e m/s \n\n", cmax,cmin);
	fprintf(fp,"\n\n ------------------ CHECK FOR GRID DISPERSION --------------------\n");
	fprintf(fp," To satisfactorily limit grid dispersion the number of gridpoints \n");

        fprintf(fp," per minimum wavelength (of P-waves) should be 6 (better more).\n");
	fprintf(fp," Here the minimum wavelength is assumed to be minimum model phase velocity \n");

        fprintf(fp," (of P-waves) at maximum frequency of the source\n");
	fprintf(fp," devided by maximum frequency of the source.\n");
	fprintf(fp," Maximum frequency of the source is approximately %8.2f Hz\n",2.0/TS);
        fprintf(fp," The minimum wavelength (of P-waves) in the following simulation will\n");
	fprintf(fp," be %e meter.\n", cmin/fmax);
	fprintf(fp," Thus, the recommended value for DH is %e meter.\n", dhstab);
	fprintf(fp," You have specified DH= %e meter.\n\n", DH);
	if (DH>dhstab)
		warning(" Grid dispersion will influence wave propagation, choose smaller grid spacing (DH).");
	
	fprintf(fp," \n\n ----------------------- CHECK FOR STABILITY ---------------------\n");
	fprintf(fp," The following simulation is stable provided that\n\n");
	fprintf(fp," \t p=cmax*DT/DH < 1/(sqrt(2)*gamma),\n\n");
	fprintf(fp," where cmax is the maximum phase velocity at infinite frequency\n");
	fprintf(fp," and gamma = sum(|FD coeff.|)\n");

	fprintf(fp," In the current simulation cmax is %8.2f m/s .\n\n",cmax);

	fprintf(fp," DT is the timestep and DH is the grid size.\n\n");
	fprintf(fp," In this simulation the stability limit for timestep DT is %e seconds .\n",dtstab);
	fprintf(fp," You have specified DT= %e s.\n", DT);
	if (DT>dtstab)
		declare_error(" The simulation will get unstable, choose smaller DT. ");
	else fprintf(fp," The simulation will be stable.\n");


        fprintf(fp," \n\n --------------------- CHECK FOR INPUT ERRORS ---------------------\n");

	if ((SEISMO)&& (MYID==0)){
		fprintf(fp," Checking the number of seismogram samples. \n");
		fprintf(fp,"    Number of timesteps %d.\n", NT);
		fprintf(fp,"    Seismogram sampling rate in timesteps %d.\n", NDT);
		fprintf(fp,"    Number of seismogram output samples %d.\n", NT/NDT);

		/* SU and SEG-Y allow 32767 samples, furthermore the exist programs allow for 65535 samples
		and pseudo-SEG-Y formats allow foralmost arbitrarily long traces. 
		For binary and textual output the limit is arbitrary. USHRT_MAX is the limit of an unsigned short specified in limits.h */

		if ((SEIS_FORMAT[0]==1) && (NT/NDT)>( USHRT_MAX )) {
			fprintf(fp," Maximum allowed number of samples per trace in SU format: %d \n", USHRT_MAX );
			declare_error(" Sorry. Too many samples per receiver! \n");
		}
	}



// 	if ((SEISMO>0) && (MYID==0)) {
// 		srec_minx=DH*NX*NPROCX+1, srec_miny=DH*NY*NPROCY+1;
// 		srec_maxx=-1.0, srec_maxy=-1.0;
// 		fprintf(fp,"\n Checking for receiver position(s) specified in input file.\n");
// 		fprintf(fp,"    Global grid size in m: %5.2f (x) : %5.2f (y) \n",NX*DH*NPROCX,NY*DH*NPROCY);
// 		if (FREE_SURF==0) fprintf(fp,"    Global grid size in m(-width of abs.boundary) : \n        %5.2f-%5.2f (x in m) : %5.2f-%5.2f (y in m) \n",FW*DH,NX*DH*NPROCX-FW*DH,FW*DH,NY*DH*NPROCY-FW*DH);
// 		if (FREE_SURF==1) fprintf(fp,"    Global grid size in m(-width of abs.boundary) : \n        %5.2f-%5.2f (x in m) : %5.2f-%5.2f (y in m) \n",FW*DH,NX*DH*NPROCX-FW*DH,DH,NY*DH*NPROCY-FW*DH);
// 
//                         /* find maximum and minimum source positions coordinate ---- from input file*/
//                         /*for usability reasons, "z" - as commonly used - denotes the depth (vertical direction),
//                   however, internally "y" is used for the vertical coordinate,
//                   we simply switch the "y" and "z" coordinate as read in the input file,
//                   therefore we determine the minimum/maximum position in y-direction by the ZREC1 variable and vice versa.
//                   this has to be considered for the receiver line coordinates specified in both the input file and separate source/receiver files*/
// 
// 		if (READREC==0 || READREC==2) {
// 			if (XREC1>XREC2) {
// 				srec_maxx=XREC1;
// 				srec_minx=XREC2;
// 			}
// 			else {
// 				srec_maxx=XREC2;
// 				srec_minx=XREC1;
// 			}
// 			if (YREC1>YREC2) {
// 				srec_maxy=YREC1;
// 				srec_miny=YREC2;
// 			}
// 			else {
// 				srec_maxy=YREC2;
// 				srec_miny=YREC1;
// 			}
// 			fprintf(fp,"    Number of receiver positions in input file : %i \n", ntr);
// 		}
// 		if (READREC==1) {
// 			/* find maximum and minimum source positions coordinate ---- from receiver file*/
// 			for (k=1;k<=ntr;k++){
// 				/* find maximum source positions coordinate*/
// 				if ((recpos[1][k]*DH)>srec_maxx) srec_maxx=recpos[1][k]*DH;
// 				if ((recpos[2][k]*DH)>srec_maxy) srec_maxy=recpos[2][k]*DH;
// 				/* find minimum source positions coordinate*/
// 				if ((recpos[1][k]*DH)<srec_minx) srec_minx=recpos[1][k]*DH;
// 				if ((recpos[2][k]*DH)<srec_miny) srec_miny=recpos[2][k]*DH;
// 			}
// 			fprintf(fp,"    Number of receiver positions in receiver file %s : %i \n", REC_FILE, ntr);
// 		}
// 
// 
// 
// 		fprintf(fp,"    Minimum receiver position coordinates : %5.2f (x) : %5.2f (y) \n",srec_minx,srec_miny);
// 		fprintf(fp,"    Maximum receiver position coordinates : %5.2f (x) : %5.2f (y) \n",srec_maxx,srec_maxy);
// 		/* checking if receiver coordinate of first receiver in line specified in input-file is inside the global grid */
// 		if (((srec_maxx<0.0) || (srec_maxy<0.0)) || ((srec_minx<0.0) || (srec_miny<0.0))) {
// 			declare_error("\n\n Coordinate of at least one receiver location is outside the global grid. \n\n");
// 		}
// 		if ((srec_maxx>NX*DH*NPROCX) || (srec_maxy>NY*DH*NPROCY)) {
// 			declare_error("\n\n Coordinate of at least one receiver location is outside the global grid. \n\n");
// 		}
// 		/* checking if receiver coordinate of first receiver in line specified in input-file is outside the Absorbing Boundary  */
// 		if ((srec_maxx<(FW*DH)) || (srec_minx<(FW*DH))) {
// 			/* this warning appears, when at least a single receiver is located in AB between 0 - FW+DX/DX/DZ ("inner boundary")*/
// 			warning("\n\n Coordinate of at least one receiver location is inside the Absorbing Boundary (warning 1, left boundary). \n\n");
// 		}
// 		if (srec_maxx>(NX*DH*NPROCX-FW*DH)) {
// 			/* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary")*/
// 			warning("\n\n Coordinate of at least one receiver location is inside the Absorbing Boundary (warning 2, right boundary). \n\n");
// 		}
// 		if (srec_maxy>(NY*DH*NPROCY-FW*DH)) {
// 			/* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary")*/
// 			warning("\n\n Coordinate of at least one receiver location is inside the Absorbing Boundary (warning 3, lower boundary). \n\n");
// 		}
// 		if ((srec_miny<(FW*DH)) && !(FREE_SURF)) {
// 			/* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary")*/
// 			warning("\n\n Coordinate of at least one receiver location is inside the Absorbing Boundary (warning 4, top boundary). \n\n");
// 		}
// 
// 		fprintf(fp," ... complete, receiver position specified in input file are located within the global grid.\n");
// 
// 	}
// 
// 
// 
// 
// 
// 	if ((SRCREC==1)&& (MYID==0)){
// 		srec_minx=DH*NX*NPROCX+1, srec_miny=DH*NY*NPROCY+1;
// 		srec_maxx=-1.0, srec_maxy=-1.0;
// 		fprintf(fp,"\n Checking for source position(s) specified in source file. \n");
// 		fprintf(fp,"    Global grid size in m: %5.2f (x) : %5.2f (y) \n",NX*DH*NPROCX,NY*DH*NPROCY);
// 		if (FREE_SURF==0) fprintf(fp,"    Global grid size in m (-width of abs.boundary) : \n        %5.2f-%5.2f (x in m) : %5.2f-%5.2f (y in m) \n",FW*DH,NX*DH*NPROCX-FW*DH,FW*DH,NY*DH*NPROCY-FW*DH);
// 		if (FREE_SURF==1) fprintf(fp,"    Global grid size in m(-width of abs.boundary) : \n        %5.2f-%5.2f (x in m) : %5.2f-%5.2f (y in m) \n",FW*DH,NX*DH*NPROCX-FW*DH,DH,NY*DH*NPROCY-FW*DH);
// 
// 
// 		for (k=1;k<=nsrc;k++){
// 			/* find maximum source positions coordinate*/
// 			if (srcpos[1][k]>srec_maxx) srec_maxx=srcpos[1][k];
// 			if (srcpos[2][k]>srec_maxy) srec_maxy=srcpos[2][k];
// 			/* find minimum source positions coordinate*/
// 			if (srcpos[1][k]<srec_minx) srec_minx=srcpos[1][k];
// 			if (srcpos[2][k]<srec_miny) srec_miny=srcpos[2][k];
// 		}
// 
// 		fprintf(fp,"    Number of source positions in source file %s. : %i.\n", SOURCE_FILE, nsrc);
// 		fprintf(fp,"    Minimum source position coordinates : %5.2f (x) : %5.2f (y) \n",srec_minx,srec_miny);
// 		fprintf(fp,"    Maximum source position coordinates : %5.2f (x) : %5.2f (y) \n",srec_maxx,srec_maxy);
// 		/* checking if receiver coordinate of first receiver in line specified in input-file is inside the global grid */
// 		if (((srec_maxx<0.0) || (srec_maxy<0.0)) || ((srec_minx<0.0) || (srec_miny<0.0))) {
// 			declare_error("\n\n Coordinate of at least one source location is outside the global grid. \n\n");
// 		}
// 		if ((srec_maxx>NX*DH*NPROCX) || (srec_maxy>NY*DH*NPROCY)) {
// 			declare_error("\n\n Coordinate of at least one source location is outside the global grid. \n\n");
// 		}
// 		/* checking if receiver coordinate of first receiver in line specified in input-file is outside the Absorbing Boundary  */
// 		if ((srec_maxx<(FW*DH)) || (srec_minx<(FW*DH))) {
// 			/* this warning appears, when at least a single receiver is located in AB between 0 - FW+DX/DX/DZ ("inner boundary")*/
// 			warning("\n\n Coordinate of at least one source location is inside the Absorbing Boundary (warning 1, left boundary). \n\n");
// 		}
// 		if (srec_maxx>(NX*DH*NPROCX-FW*DH)) {
// 			/* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary")*/
// 			warning("\n\n Coordinate of at least one source location is inside the Absorbing Boundary (warning 2, right boundary). \n\n");
// 		}
// 		if (srec_maxy>(NY*DH*NPROCY-FW*DH)) {
// 			/* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary")*/
// 			warning("\n\n Coordinate of at least one source location is inside the Absorbing Boundary (warning 3, lower boundary). \n\n");
// 		}
// 		if ((srec_miny<(FW*DH)) && !(FREE_SURF)) {
// 			/* this warning appears, when at least a single receiver is located in AB between NX/NY/NZ-FW+DX/DX/DZ and NX/NY/NZ ("outer boundary")*/
// 			warning("\n\n Coordinate of at least one source location is inside the Absorbing Boundary (warning 4, top boundary). \n\n");
// 		}
// 		fprintf(fp," ...complete, all source position(s) specified in source file are located within the global grid.\n");
// 
// 	}





	fprintf(fp,"\n\n ----------------------- ABSORBING BOUNDARY ------------------------\n");
        if((FW>nx_min)||(FW>ny_min)){
	  declare_error(" The width of the absorbing boundary is larger than one computational domain. Choose smaller FW or use less CPUs.");
	}

	fprintf(fp," Width (FW) of absorbing frame should be at least 10 gridpoints.\n");
	fprintf(fp," You have specified a width of %d gridpoints.\n",FW);
	if (FW<10) 
		warning(" Be aware of artificial reflections from grid boundaries ! \n");

	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
}

