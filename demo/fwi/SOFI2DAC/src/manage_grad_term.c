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

#include "fd.h"

void calc_correlation(float **fvx, float **fvy, float **fsp,
                      float **bvx, float **bvy, float **bsp,
                      float **cvxxyy, float **csp){
    /* [fvx, fvy, fsp] means the time derivative of vx, vy and pressure 
     * in source wavefield; and [bvx, bvy, bsp] means vx, vy and pressure 
     * in receiver wavefield
     */
    extern int NX1, NX2, NY1, NY2, IDX, IDY;
    int i, j;
    
    /* local variables */
    for (j=NY1;j<=NY2;j=j+IDY){
    for (i=NX1;i<=NX2;i=i+IDX){
        cvxxyy[j][i]+=fvx[j][i]*bvx[j][i]+fvy[j][i]*bvy[j][i];
        csp[j][i]+=fsp[j][i]*bsp[j][i];
    }}
}

void calc_gradient(float **cvxxyy, float **csp, 
                   float **prho, float **ppi){
/* 
 * in sofi2dac here, ppi=rho*vp*vp (bulk modulus), prho=rho (density)
 * in my IFOS2D p_mod->ppi means velocity and prho is same with sofi2dac
 */

    extern int NX1, NX2, NY1, NY2, IDX, IDY;
    extern float DT;
    /* local variables */
    int i, j;
    float Mvp, Mrho, lam, lam2, csp1, cvxxyy1;
    
    for (j=NY1;j<=NY2;j=j+IDY){
    for (i=NX1;i<=NX2;i=i+IDX){
        Mvp=sqrt(ppi[j][i]/prho[j][i]);
        Mrho=prho[j][i];
        lam=ppi[j][i]; //lam=Mrho[j][i]*Mvp[j][i]*Mvp[j][i];
        lam2=lam*lam;
        /* ----------------------------------------- */
        /* calculate gradient direction lamda mu rho */
        /* ----------------------------------------- */
        /* calculate gradient direction lamda */
        csp[j][i]=-DT*(0.25*csp[j][i]/lam2);
        /* calculate gradient direction rho */
        cvxxyy[j][i]=-DT*cvxxyy[j][i];
        /* ----------------------------------------- */
        /* calculate complete gradient direction     */
        /* ----------------------------------------- */
        csp1=csp[j][i];
        cvxxyy1=cvxxyy[j][i];
        /* calculate gradient direction vp */
        csp[j][i] = 2.0*Mvp*Mrho*csp1;
        /* calculate gradient direction rho */
        cvxxyy[j][i] = (Mvp*Mvp)*csp1+cvxxyy1;
    }}
}

void eprecond(float ** W, float ** vx, float ** vy){
    /*here we use the velocity compoments, you also can use acceleration */
    extern int NX1, NX2, NY1, NY2, IDX, IDY;
    int i, j;
    for (j=NY1;j<=NY2;j=j+IDY){
        for (i=NX1;i<=NX2;i=i+IDX){       
            W[j][i]+=((vx[j][i]*vx[j][i])+(vy[j][i]*vy[j][i]));
        }
    }
}

void eprecond1(float ** We, float ** Ws, float ** Wr){
    extern int NX1, NX2, NY1, NY2, IDX, IDY;
    int i, j;
    /* calculate energy weighting */
    /* energy weighted source and receiver contribution */
    for (j=NY1;j<=NY2;j=j+IDY){
        for (i=NX1;i<=NX2;i=i+IDX){        
            We[j][i]=sqrt(Ws[j][i]*Wr[j][i]);
        }
    }
}

void update_precond(float **eng, float prctile){
    
    extern int NXG, NYG, IDX, IDY;
    int i, j;
    /*update and apply preconditioning*/
    for (j=1;j<=NYG;j=j+IDY){
    for (i=1;i<=NXG;i=i+IDX){
        eng[j][i]+=prctile;
    }}
}

void apply_precond(float **eng, float **grad){
    
    extern int NXG, NYG, IDX, IDY;
    int i, j;
    /*update and apply preconditioning*/
    for (j=1;j<=NYG;j=j+IDY){
    for (i=1;i<=NXG;i=i+IDX){
        grad[j][i]=grad[j][i]/eng[j][i];
    }}
}

/* Return percentiles of the values in a matrix for the percentages p1 and p2 in the interval [0,100]. */
float global_percentiles(float p1, float ** fullmod) {
    /* btmp=global_percentiles(NXG, NYG, NX, NY, IDXI, IDYI, POS[1], POS[2], EPSILON_NUM_LOW, We, MYID_SHOT, MPI_COMM_SHOT); */
    extern int NXG, NYG, NX, NY, IDX, IDY, POS[3], MYID;
    int i,j,ij,ij1;
    int numx, numy, num;
    float *fullvec;
    float prctile;
    
    if(p1<0.0||p1>1.0){
        printf("p1: %f from function global_percentiles\n",p1);
        declare_error("p1 should be in (0.0 1.0)\n");
    }
    
    numx=(NXG-1)/IDX+1;
    numy=(NYG-1)/IDY+1;
    num = numx*numy;
    
    fullvec = vector(1,num);
    ij=0;
    for (j=1;j<=NYG;j=j+IDY){
    for (i=1;i<=NXG;i=i+IDX){
        ij=ij+1;
        fullvec[ij]= fullmod[j][i];
    }}
    
    quicksort2(fullvec,1,num);
    //////////
    ij1=find_min_positive_num2(1,num-1,fullvec);
    ij1=ij1-1+(int)(num*p1);
    //ij1=(int)(num*p1);
    //////////
    prctile=fullvec[ij1];
    
    free_vector(fullvec,1,num);
    printf("EPSILON_NUM_LOW=%e,  ---EPSILON_WE=%e---\n",p1,prctile);
    
    //printf("prctile=%e \n",prctile);
    
    return prctile;
};

/*  quickSort
 * This public-domain C implementation by Darel Rex Finley .
 * http://alienryderflex.com/quicksort */

void quicksort2(float *arr, int dummy, int elements) {
    
    /* number of levels: 64 is enough for 1.8e19 elements to be sorted */
    #define  MAX_LEVELS  64
    
    int	beg[MAX_LEVELS], end[MAX_LEVELS], i=0, L, R, swap ;
    float	piv;
    
    beg[0]=0;
    end[0]=elements;
    
    while (i>=0) {
        L=beg[i]; R=end[i]-1;
        if (L<R) {
            piv=arr[L];
            while (L<R) {
                while (arr[R]>=piv && L<R) R--;
                if (L<R) arr[L++]=arr[R];
                while (arr[L]<=piv && L<R) L++;
                if (L<R) arr[R--]=arr[L];
            }
            arr[L]=piv; beg[i+1]=L+1; end[i+1]=end[i]; end[i++]=L;
            if (end[i]-beg[i]>end[i-1]-beg[i-1]) {
                swap=beg[i]; beg[i]=beg[i-1]; beg[i-1]=swap;
                swap=end[i]; end[i]=end[i-1]; end[i-1]=swap;
            }
        }
        else	i--;
    }
}

int find_min_positive_num2(int n1,int n2,float *a){
    /*find the element position in which is the minimum positive value  */
    int i,mid,key,left,right;
    
    if((a[n1]<0.0)||(a[n2]<=0.0))   declare_error("element value of array should be more than 0\n");
    if(a[n1]>0.0){
        mid=n1;
        return mid;
    }
    left=n1,    right=n2;
    mid=(left+right)/2;
    while((left<right)&&(mid!=left)){
        if(a[mid]>0.0){
            left=left,    right=mid;
            mid=(left+right)/2;
        }else if(a[mid]==0.0){
            left=mid,    right=right;
            mid=(left+right)/2;
        }
    }
    return mid+1;
}

void savematIDXY(float **mod,char modfile[STRING_SIZE]){
    
    extern int NXG, NYG, NX, NY, IDX, IDY, POS[3], MYID;
    float **fullmod;
    FILE *fp;
    int i,j;

    fullmod = matrix(1,NYG,1,NXG);
    matrix_cat(NXG, NYG, NX, NY, POS[1], POS[2], mod, fullmod, MPI_COMM_WORLD);
    if(MYID==0){
        fp=fopen(modfile,"wb");
        for (i=1;i<=NXG;i=i+IDX){
        for (j=1;j<=NYG;j=j+IDY){
            fwrite(&fullmod[j][i],sizeof(float),1,fp);   
        }}
        fclose(fp);
    }
    free_matrix(fullmod, 1,NYG,1,NXG);
    MPI_Barrier(MPI_COMM_WORLD);
}

void savematIDXY2(float **fullmod,char modfile[STRING_SIZE]){
    
    extern int NXG, NYG, IDX, IDY;
    FILE *fp;
    int i,j;
    fp=fopen(modfile,"wb");
    for (i=1;i<=NXG;i=i+IDX){
    for (j=1;j<=NYG;j=j+IDY){
        fwrite(&fullmod[j][i],sizeof(float),1,fp);   
    }}
    fclose(fp);
}

void matrix_cat(int nxg, int nyg, int nx, int ny, int pos1, int pos2, 
                float **mod, float **fullmod, MPI_Comm mpi_comm_shot){
    /* !!!!just let nodes from one shot communictor(mpi_comm_shot) call this function!!!! */
    /*for example: matrix_cat(nxg, nyg, nx, ny, pos1, pos2, matrix1, fullmod, mpi_comm_shot); */
    
    int i,j,ii,jj,POS2NY,POS1NX;
    float **fullmod2;

    POS1NX=pos1*nx;
    POS2NY=pos2*ny;

    /* temporary global mod array for MPI-exchange */
    fullmod2 = matrix(1,nyg,1,nxg);
    
    /* loop over global model: copy model of local array	*/
    /* to appropriate locations in the global array		*/
    for(j=1;j<=ny;j++){
        jj=j+POS2NY;
        for(i=1;i<=nx;i++){
            ii=i+POS1NX;
            fullmod2[jj][ii] = mod[j][i];
        }
    }

    MPI_Allreduce(&fullmod2[1][1], &fullmod[1][1], nxg*nyg, MPI_FLOAT, MPI_SUM, mpi_comm_shot);
    
    free_matrix(fullmod2, 1,nyg,1,nxg);
    MPI_Barrier(mpi_comm_shot);
}

void matrix_scatter(int nxg, int nyg, int nx, int ny, int pos1, int pos2,
                    float **mod, float **fullmod){
    /* !!!!it should be called after fullmod is same at every node!!!! */
    /* distribute global spatial models on computational nodes */
    /* scattermodel(pi, fullmod); */

    int i,j,ii,jj,POS2NY,POS1NX;

    POS1NX=pos1*nx;
    POS2NY=pos2*ny;

    /* loop over global model: copy appropriate locations in the global array	*/
    /* to model of local array	*/
    for(j=1;j<=ny;j++){
        jj=j+POS2NY;
        for(i=1;i<=nx;i++){
            ii=i+POS1NX;
            mod[j][i] = fullmod[jj][ii];
        }
    }
}

void matrix_norm(int nx1, int nx2,int ny1, int ny2,int idx, int idy, float ** matrix){
    int i,j;
    float max=0.0;
    //global maximum
    max = fabs(matrix[ny1][nx1]);
    for (j=ny1;j<=ny2;j=j+idy){
    for (i=nx1;i<=nx2;i=i+idx){
        if(fabs(matrix[j][i]) > max) {
            max = fabs(matrix[j][i]);
        }
    }}
    max=max;
    if(max==0.0){
        printf("Please check input parameters of function matrix_norm ");
        exit(1);
    }
    for (j=ny1;j<=ny2;j=j+idy){
    for (i=nx1;i<=nx2;i=i+idx){
        matrix[j][i]=matrix[j][i]/max;
    }}
}

