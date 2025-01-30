/*-------------------------------------------------------------------
 * Copyright (C) 2018  For the list of authors, see file AUTHORS.
 *
 * This file is part of SEISPLATFORM.
 * 
 * SEISPLATFORM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published 
 * by the Free Software Foundation, version 2.0 of the License only.
 * 
 * SEISPLATFORM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with SEISPLATFORM. See file COPYING and/or 
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
-------------------------------------------------------------------*/

#include "comm.h"
/*
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
*/
#define SHIFT_IND 1
#define FREE_ARGUMENT char *

void declare_error(char err_text[]){
    //extern int MYID;
    
    //fprintf(stdout,"Message from PE %d\n",MYID);
    fprintf(stdout,"R U N - T I M E  E R R O R: \n");
    fprintf(stdout,"%s\n",err_text);
    fprintf(stdout,"...now exiting to system.\n");
    
    //MPI_Abort(MPI_COMM_WORLD, 1);
    exit(1);
}

void throw_error(char err_text[]){
    //extern int MYID;
    
    //fprintf(stdout,"Message from PE %d\n",MYID);
    fprintf(stdout,"R U N - T I M E  E R R O R: \n");
    fprintf(stdout,"%s\n",err_text);
    fprintf(stdout,"...now exiting to system.\n");
    
    //MPI_Abort(MPI_COMM_WORLD, 1);
    exit(1);
}

// void warning(char warn_text[]){
//     /* standard warnings handler */
//     fprintf(stdout,"W A R N I N G   M E S S A G E: \n");
//     fprintf(stdout,"%s\n",warn_text);
// }


double maximum(float **a, int nx, int ny){
    double maxi=0.0;
    int i, j;
    
    
    for (j=1;j<=ny;j++)
        for (i=1;i<=nx;i++)
            if (fabs((double) a[i][j])>maxi) maxi=fabs((double)a[i][j]);
    return maxi;
}

float minimum_m(float **mat, int nx, int ny)
{
    float minm;
    int i,j;
    
    minm = mat[1][1];
    for (i=1;i<=nx;i++)
        for (j=1;j<=ny;j++)
        {
            if (i*j==1) continue;
            if (mat[j][i] < minm)
            {
                minm = mat[j][i];
            }
        }
    return minm;
}

float maximum_m(float **mat, int nx, int ny)
{
    float maxm;
    int i,j;
    
    maxm = mat[1][1];
    for (i=1;i<=nx;i++)
        for (j=1;j<=ny;j++)
        {
            if (i*j==1) continue;
            if (mat[j][i] > maxm)
            {
                maxm = mat[j][i];
            }
        }
    return maxm;
}




float *vector(int ni, int nj){
    float *a;
    int k;
    
    a=(float *)malloc((size_t) ((nj-ni+1+SHIFT_IND)*sizeof(float)));
    if (!a) declare_error("util.c: allocation failure in function vector()");
    for (k=0;k<(nj-ni+1+SHIFT_IND);k++) a[k]=0.0;
    return a-ni+SHIFT_IND;
}


int *ivector(int ni, int nj){
    
    int *a;
    int k;
    
    a=(int *)malloc((size_t) ((nj-ni+1+SHIFT_IND)*sizeof(int)));
    if (!a) declare_error("util.c: allocation failure in function ivector()");
    for (k=0;k<(nj-ni+1+SHIFT_IND);k++) a[k]=0;
    return a-ni+SHIFT_IND;
}

unsigned short int *usvector(int ni, int nj){
    
    unsigned short int *a;
    int k;
    
    a=(unsigned short int *)malloc((size_t) ((nj-ni+1+SHIFT_IND)*sizeof(unsigned short int)));
    if (!a) declare_error("util.c: allocation failure in function usvector()");
    for (k=0;k<(nj-ni+1+SHIFT_IND);k++) a[k]=0;
    return a-ni+SHIFT_IND;
}


unsigned char *cvector(int ni, int nj){
    
    unsigned char *a;
    
    a=(unsigned char *)malloc((size_t) ((nj-ni+1+SHIFT_IND)*sizeof(unsigned char)));
    if (!a) declare_error("util.c: allocation failure in function cvector()");
    return a-ni+SHIFT_IND;
}


unsigned long *lvector(int ni, int nj){
    
    unsigned long *a;
    int k;
    
    a=(unsigned long *)malloc((size_t) ((nj-ni+1+SHIFT_IND)*sizeof(unsigned long)));
    if (!a) declare_error("util.c: allocation failure in function lvector()");
    for (k=0;k<(nj-ni+1+SHIFT_IND);k++) a[k]=0;
    return a-ni+SHIFT_IND;
}

double *dvector(int ni, int nj){
    
    double *a;
    int k;
    
    a=(double *)malloc((size_t) ((nj-ni+1+SHIFT_IND)*sizeof(double)));
    if (!a) declare_error("util.c: allocation failure in function dvector()");
    for (k=0;k<(nj-ni+1+SHIFT_IND);k++) a[k]=0.0;
    return a-ni+SHIFT_IND;
}

float **fmatrix(int mrl, int mrh, int mcl, int mch){
    
    int k,l, mrow=mrh-mrl+1,mcol=mch-mcl+1;
    float **ma;
    
    ma=(float **) malloc((size_t) ((mrow+SHIFT_IND)*sizeof(float*)));
    if (!ma) declare_error("util.c: allocation failure 1 in function fmatrix() ");
    ma += SHIFT_IND;
    ma -= mrl;
    
    ma[mrl]=(float *) malloc((size_t)((mrow*mcol+SHIFT_IND)*sizeof(float)));
    if (!ma[mrl]) declare_error("util.c: allocation failure 2 in function fmatrix() ");
    ma[mrl] += SHIFT_IND;
    ma[mrl] -= mcl;
    
    for (k=mrl+1;k<=mrh;k++) ma[k]=ma[k-1]+mcol;
    
    for (k=mrl;k<=mrh;k++)
        for (l=mcl;l<=mch;l++) ma[k][l]=0.0;
    
    return ma;
}

float **matrix(int mrl, int mrh, int mcl, int mch){
    
    int k,l, mrow=mrh-mrl+1,mcol=mch-mcl+1;
    float **ma;
    
    ma=(float **) malloc((size_t) ((mrow+SHIFT_IND)*sizeof(float*)));
    if (!ma) declare_error("util.c: allocation failure 1 in function matrix() ");
    ma += SHIFT_IND;
    ma -= mrl;
    
    ma[mrl]=(float *) malloc((size_t)((mrow*mcol+SHIFT_IND)*sizeof(float)));
    if (!ma[mrl]) declare_error("util.c: allocation failure 2 in function matrix() ");
    ma[mrl] += SHIFT_IND;
    ma[mrl] -= mcl;
    
    for (k=mrl+1;k<=mrh;k++) ma[k]=ma[k-1]+mcol;
    
    for (k=mrl;k<=mrh;k++)
        for (l=mcl;l<=mch;l++) ma[k][l]=0.0;
    
    return ma;
}


double **dmatrix(int mrl, int mrh, int mcl, int mch){
    
    int k,l, mrow=mrh-mrl+1,mcol=mch-mcl+1;
    double **ma;
    
    ma=(double **) malloc((size_t) ((mrow+SHIFT_IND)*sizeof(double*)));
    if (!ma) declare_error("util.c: allocation failure 1 in function matrix() ");
    ma += SHIFT_IND;
    ma -= mrl;
    
    ma[mrl]=(double *) malloc((size_t)((mrow*mcol+SHIFT_IND)*sizeof(double)));
    if (!ma[mrl]) declare_error("util.c: allocation failure 2 in function dmatrix() ");
    ma[mrl] += SHIFT_IND;
    ma[mrl] -= mcl;
    
    for (k=mrl+1;k<=mrh;k++) ma[k]=ma[k-1]+mcol;
    
    for (k=mrl;k<=mrh;k++)
        for (l=mcl;l<=mch;l++) ma[k][l]=0.0;
    
    return ma;
}

int **imatrix(int mrl, int mrh, int mcl, int mch){
    
    int k,l, mrow=mrh-mrl+1,mcol=mch-mcl+1;
    int **ma;
    
    ma=(int **) malloc((size_t) ((mrow+SHIFT_IND)*sizeof(int*)));
    if (!ma) declare_error("util.c: allocation failure 1 in function imatrix() ");
    ma += SHIFT_IND;
    ma -= mrl;
    
    ma[mrl]=(int *) malloc((size_t)((mrow*mcol+SHIFT_IND)*sizeof(int)));
    if (!ma[mrl]) declare_error("util.c: allocation failure 2 in function imatrix() ");
    ma[mrl] += SHIFT_IND;
    ma[mrl] -= mcl;
    
    for (k=mrl+1;k<=mrh;k++) ma[k]=ma[k-1]+mcol;
    
    for (k=mrl;k<=mrh;k++)
        for (l=mcl;l<=mch;l++) ma[k][l]=0;
    
    return ma;
}

unsigned short int **usmatrix(int mrl, int mrh, int mcl, int mch){
    
    int k,l, mrow=mrh-mrl+1,mcol=mch-mcl+1;
    unsigned short int **ma;
    
    ma=(unsigned short int **) malloc((size_t) ((mrow+SHIFT_IND)*sizeof(unsigned short int*)));
    if (!ma) declare_error("util.c: allocation failure 1 in function usmatrix() ");
    ma += SHIFT_IND;
    ma -= mrl;
    
    ma[mrl]=(unsigned short int *) malloc((size_t)((mrow*mcol+SHIFT_IND)*sizeof(unsigned short int)));
    if (!ma[mrl]) declare_error("util.c: allocation failure 2 in function usmatrix() ");
    ma[mrl] += SHIFT_IND;
    ma[mrl] -= mcl;
    
    for (k=mrl+1;k<=mrh;k++) ma[k]=ma[k-1]+mcol;
    
    for (k=mrl;k<=mrh;k++)
        for (l=mcl;l<=mch;l++) ma[k][l]=0;
    
    return ma;
}

float ***f3tensor(int mrl, int mrh, int mcl, int mch,int mdl, int mdh){
    
    int w,e,r, mrow=mrh-mrl+1,mcol=mch-mcl+1,mdep=mdh-mdl+1;
    float ***te;
    
    te=(float ***) malloc((size_t) ((mrow+SHIFT_IND)*sizeof(float**)));
    if (!te) declare_error("util.c: allocation failure 1 in function f3tensor() ");
    te += SHIFT_IND;
    te -= mrl;
    
    te[mrl]=(float **) malloc((size_t)((mrow*mcol+SHIFT_IND)*sizeof(float*)));
    if (!te[mrl]) declare_error("util.c: allocation failure 2 in function f3tensor() ");
    te[mrl] += SHIFT_IND;
    te[mrl] -= mcl;
    
    te[mrl][mcl]=(float *) malloc((size_t)((mrow*mcol*mdep+SHIFT_IND)*sizeof(float)));
    if (!te[mrl][mcl]) declare_error("util.c: allocation failure 3 in function f3tensor() ");
    te[mrl][mcl] += SHIFT_IND;
    te[mrl][mcl] -= mdl;
    
    for (e=mcl+1;e<=mch;e++) te[mrl][e]=te[mrl][e-1]+mdep;
    for (w=mrl+1;w<=mrh;w++){
        te[w]=te[w-1]+mcol;
        te[w][mcl]=te[w-1][mcl]+mcol*mdep;
        for (e=mcl+1;e<=mch;e++) te[w][e]=te[w][e-1]+mdep;
    }
    
    for (w=mrl;w<=mrh;w++)
        for (e=mcl;e<=mch;e++)
            for (r=mdl;r<=mdh;r++) te[w][e][r]=0.0;
    
    return te;
}


int ***i3tensor(int mrl, int mrh, int mcl, int mch,int mdl, int mdh){
    
    int w,e,r, mrow=mrh-mrl+1,mcol=mch-mcl+1,mdep=mdh-mdl+1;
    int ***te;
    
    te=(int ***) malloc((size_t) ((mrow+SHIFT_IND)*sizeof(int**)));
    if (!te) declare_error("util.c: allocation failure 1 in function i3tensor() ");
    te += SHIFT_IND;
    te -= mrl;
    
    te[mrl]=(int **) malloc((size_t)((mrow*mcol+SHIFT_IND)*sizeof(int*)));
    if (!te[mrl]) declare_error("util.c: allocation failure 2 in function i3tensor() ");
    te[mrl] += SHIFT_IND;
    te[mrl] -= mcl;
    
    te[mrl][mcl]=(int *) malloc((size_t)((mrow*mcol*mdep+SHIFT_IND)*sizeof(int)));
    if (!te[mrl][mcl]) declare_error("util.c: allocation failure 3 in function i3tensor() ");
    te[mrl][mcl] += SHIFT_IND;
    te[mrl][mcl] -= mdl;
    
    for (e=mcl+1;e<=mch;e++) te[mrl][e]=te[mrl][e-1]+mdep;
    for (w=mrl+1;w<=mrh;w++){
        te[w]=te[w-1]+mcol;
        te[w][mcl]=te[w-1][mcl]+mcol*mdep;
        for (e=mcl+1;e<=mch;e++) te[w][e]=te[w][e-1]+mdep;
    }
    
    for (w=mrl;w<=mrh;w++)
        for (e=mcl;e<=mch;e++)
            for (r=mdl;r<=mdh;r++) te[w][e][r]=0.0;
    
    
    return te;
}
float ****f4tensor(int nrl, int nrh, int ncl, int nch,int ndl, int ndh, int nvl, int nvh){
	/* allocate a float 4tensor with subscript range m[nrl..nrh][ncl..nch][ndl..ndh][nvl..nvh]
		   and intializing the matrix, e.g. m[nrl..nrh][ncl..nch][ndl..ndh][nvl..nvh]=0.0 */
	int i,j,d,v, nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1, nval=nvh-nvl+1;
	float ****t;

	/* allocate pointers to pointers to rows */
	t=(float ****) malloc((size_t) ((nrow+SHIFT_IND)*sizeof(float***)));
	if (!t) declare_error("allocation failure 1 in function f4tensor() ");
	t += SHIFT_IND;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(float ***) malloc((size_t)((nrow*ncol+SHIFT_IND)*sizeof(float**)));
	if (!t[nrl]) declare_error("allocation failure 2 in function f4tensor() ");
	t[nrl] += SHIFT_IND;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(float **) malloc((size_t)((nrow*ncol*ndep+SHIFT_IND)*sizeof(float*)));
	if (!t[nrl][ncl]) declare_error("allocation failure 3 in function f4tensor() ");
	t[nrl][ncl] += SHIFT_IND;
	t[nrl][ncl] -= ndl;

	/* allocate values and set pointers to them */
	t[nrl][ncl][ndl]=(float *) malloc((size_t)((nrow*ncol*ndep*nval+SHIFT_IND)*sizeof(float)));
	if (!t[nrl][ncl][ndl]) declare_error("allocation failure 4 in function f4tensor() ");
	t[nrl][ncl][ndl] += SHIFT_IND;
	t[nrl][ncl][ndl] -= nvl;

    for(i=nrl;i<=nrh;i++){
		if (i > nrl){
			t[i] = t[i-1] + ncol ;
		    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		    t[i][ncl][ndl] = t[i-1][ncl][ndl] + ncol*ndep*nval ;
		}
		for(j=ncl;j<=nch;j++){
			if (j > ncl){
				t[i][j]=t[i][j-1] + ndep ;
				t[i][j][ndl] = t[i][j-1][ndl] + ndep*nval ;
			}
			for(d=ndl;d<=ndh;d++){
				if (d > ndl) t[i][j][d] = t[i][j][d-1] + nval ;
			}
		}
	}

	/* initializing 4tensor */
	for (i=nrl;i<=nrh;i++)
		for (j=ncl;j<=nch;j++)
			for (d=ndl;d<=ndh;d++)
				for(v=nvl;v<=nvh;v++)  t[i][j][d][v]=0.0;
	
		/* return pointer to array of pointer to rows */
	return t;
}

char **chmatrix(int mrl, int mrh, int mcl, int mch){
    int k,l, mrow=mrh-mrl+1,mcol=mch-mcl+1;
    char **ma;
    
    ma=(char **) malloc((size_t) ((mrow+SHIFT_IND)*sizeof(char*)));
    if (!ma) throw_error("util.c: allocation failure 1 in function fmatrix() ");
    ma += SHIFT_IND;
    ma -= mrl;
    
    ma[mrl]=(char *) malloc((size_t)((mrow*mcol+SHIFT_IND)*sizeof(char)));
    if (!ma[mrl]) throw_error("util.c: allocation failure 2 in function fmatrix() ");
    ma[mrl] += SHIFT_IND;
    ma[mrl] -= mcl;
    
    for (k=mrl+1;k<=mrh;k++) ma[k]=ma[k-1]+mcol;
    
    for (k=mrl;k<=mrh;k++)
        for (l=mcl;l<=mch;l++) ma[k][l]='\0';
    
    return ma;
    
}

void free_vector(float *a, int ni, int nj){
    free((FREE_ARGUMENT) (a+ni-SHIFT_IND));
}

void free_ivector(int *a, int ni, int nj){
    free((FREE_ARGUMENT) (a+ni-SHIFT_IND));
}

void free_cvector(char *a, int ni, int nj){
    free((FREE_ARGUMENT) (a+ni-SHIFT_IND));
}

void free_dvector(double *a, int ni, int nj){
    free((FREE_ARGUMENT) (a+ni-SHIFT_IND));
}

void free_matrix(float **ma, int mrl, int mrh, int mcl, int mch){
    free((FREE_ARGUMENT) (ma[mrl]+mcl-SHIFT_IND));
    free((FREE_ARGUMENT) (ma+mrl-SHIFT_IND));
}

void free_imatrix(int **ma, int mrl, int mrh, int mcl, int mch){
    free((FREE_ARGUMENT) (ma[mrl]+mcl-SHIFT_IND));
    free((FREE_ARGUMENT) (ma+mrl-SHIFT_IND));
}

void free_usmatrix(unsigned short int **ma, int mrl, int mrh, int mcl, int mch){
    free((FREE_ARGUMENT) (ma[mrl]+mcl-SHIFT_IND));
    free((FREE_ARGUMENT) (ma+mrl-SHIFT_IND));
}

void free_f3tensor(float ***te, int mrl, int mrh, int mcl, int mch, int mdl, int mdh){
    free((FREE_ARGUMENT) (te[mrl][mcl]+mdl-SHIFT_IND));
    free((FREE_ARGUMENT) (te[mrl]+mcl-SHIFT_IND));
    free((FREE_ARGUMENT) (te+mrl-SHIFT_IND));
}

void free_i3tensor(int ***te, int mrl, int mrh, int mcl, int mch, int mdl, int mdh){
    free((FREE_ARGUMENT) (te[mrl][mcl]+mdl-SHIFT_IND));
    free((FREE_ARGUMENT) (te[mrl]+mcl-SHIFT_IND));
    free((FREE_ARGUMENT) (te+mrl-SHIFT_IND));
}

void free_f4tensor(float ****t, int nrl, int nrh, int ncl, int nch, int ndl, int ndh, int nvl,int nvh){
	/* free a float matrix allocated by f4tensor() */
	free((FREE_ARGUMENT) (t[nrl][ncl][ndl]+nvl-SHIFT_IND));
	free((FREE_ARGUMENT) (t[nrl][ncl]+ndl-SHIFT_IND));
	free((FREE_ARGUMENT) (t[nrl]+ncl-SHIFT_IND));
	free((FREE_ARGUMENT) (t+nrl-SHIFT_IND));
}

void free_chmatrix(char **ma, int mrl, int mrh, int mcl, int mch){
    free((FREE_ARGUMENT) (ma[mrl]+mcl-SHIFT_IND));
    free((FREE_ARGUMENT) (ma+mrl-SHIFT_IND));
}

void zero(float *A, int u_max){
    int u=0;
    while(++u<=u_max) *(A++)=0.0;
}

void normalize_data(float **data, int ntr, int ns){
    
    float *max_min_trace=NULL;
    float max, min, tmp=0.0;
    int i,j;
    
    
    max_min_trace = vector(1,ntr);
    for(i=1;i<=ntr;i++){
        
      		max=0.0;
      		min=0.0;
        
        for(j=2;j<=ns;j++){
            /* Looking for max and min */
            if(data[i][j]>max){max = data[i][j];}
            if(data[i][j]<min){min = data[i][j];}
            
            if (max>(fabs(min))){tmp=max;}
            else{tmp=fabs(min);}
        }
        
        max_min_trace[i]=tmp;
        
        /* set maximum_values=0.0 to 1.0 */
        if(max_min_trace[i]==0.0){max_min_trace[i]=1.0;}
        
        
        for(j=2;j<=ns;j++){
            data[i][j] = data[i][j]/max_min_trace[i];
        }
        
    }
}