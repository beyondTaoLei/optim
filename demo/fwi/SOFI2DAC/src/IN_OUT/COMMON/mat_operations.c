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
/* files to include */
#include "comm.h"

float M_max_abs(int nx1, int nx2,int ny1, int ny2, float **x){
    int i, j;
    float max_x;

    max_x=0.0;
    for (j=ny1;j<=ny2;j++){
    for (i=nx1;i<=nx2;i++){
        if(fabs(x[j][i])>max_x)    max_x=fabs(x[j][i]);
    }}
    
    return max_x;
}

void M_normalize(int nx1, int nx2,int ny1, int ny2, float **x){
    /* normalize matrix with its absolute maximum */
    int i, j;
    float max1;
    
    max1 = M_max_abs(nx1, nx2,ny1, ny2, x);
    if(max1!=0.0){
        for (j=ny1;j<=ny2;j++){
        for (i=nx1;i<=nx2;i++){
            x[j][i] = x[j][i]/max1;
        }}
    }
}

float M_normL2(int nx1, int nx2,int ny1, int ny2, float **x){
    /* The routine M_normL2 returns the Euclidian
       norm of a vector x of size (n2-n1+1)   */
    int i, j;
    float norm_x;

    norm_x=0.0;
    for (j=ny1;j<=ny2;j++){
    for (i=nx1;i<=nx2;i++){
        norm_x=norm_x+x[j][i]*x[j][i];
    }}
    norm_x=sqrt(norm_x);
    
    return norm_x;
    
}

/**********************************/
/**********************************/
/* operate matrixes element-wise.  */
/**********************************/
/**********************************/
void M_init(int nx1, int nx2,int ny1, int ny2, float **x){
    /* x=0.0 */
    int i, j;

    for (j=ny1;j<=ny2;j++){
    for (i=nx1;i<=nx2;i++){
        x[j][i]=0.0;
    }}
}

void M_copy(int nx1, int nx2,int ny1, int ny2, float **x, float **r){
    /* r=x */
    int i, j;

    for (j=ny1;j<=ny2;j++){
    for (i=nx1;i<=nx2;i++){
        r[j][i]=x[j][i];
    }}
}

void M_add(int nx1, int nx2,int ny1, int ny2, float **x, float **y,float **r){
    /* r=x+y */
    int i, j;

    for (j=ny1;j<=ny2;j++){
    for (i=nx1;i<=nx2;i++){
        r[j][i]=x[j][i]+y[j][i];
    }}
}

void M_subtract(int nx1, int nx2,int ny1, int ny2, float **x, float **y,float **r){
    /* r=x-y */
    int i, j;

    for (j=ny1;j<=ny2;j++){
    for (i=nx1;i<=nx2;i++){
        r[j][i]=x[j][i]-y[j][i];
    }}
}

void M_multiply(int nx1, int nx2,int ny1, int ny2, float **x, float **y,float **r){
    /* r=x*y */
    int i, j;

    for (j=ny1;j<=ny2;j++){
    for (i=nx1;i<=nx2;i++){
        r[j][i]=x[j][i]*y[j][i];
    }}
}

void M_divide(int nx1, int nx2,int ny1, int ny2, float **x, float **y,float **r){
    /* r=x/y */
    int i, j;

    for (j=ny1;j<=ny2;j++){
    for (i=nx1;i<=nx2;i++){
        r[j][i]=x[j][i]/y[j][i];
    }}
}

void M_square_sum(int nx1, int nx2,int ny1, int ny2, float **x, float **y,float **r){
    /* r=x*x+y*y */
    int i, j;

    for (j=ny1;j<=ny2;j++){
    for (i=nx1;i<=nx2;i++){
        r[j][i]=x[j][i]*x[j][i]+y[j][i]*y[j][i];
    }}
}

void M_square_sum_a(int nx1, int nx2,int ny1, int ny2, float **x, float **y,float **r){
    /* r+=x*x+y*y */
    int i, j;

    for (j=ny1;j<=ny2;j++){
    for (i=nx1;i<=nx2;i++){
        r[j][i]+=(x[j][i]*x[j][i]+y[j][i]*y[j][i]);
    }}
}

void M_AXPY(int nx1, int nx2,int ny1, int ny2, float a, float **x, float **y,float **r){
    /* r=a*x+y */
    int i, j;

    for (j=ny1;j<=ny2;j++){
    for (i=nx1;i<=nx2;i++){
        r[j][i]=a*x[j][i]+y[j][i];
    }}
}

void M_AXPb(int nx1, int nx2,int ny1, int ny2, float a, float **x,  float b, float **r){
    /* r=ax+b */
    int i, j;

    for (j=ny1;j<=ny2;j++){
    for (i=nx1;i<=nx2;i++){
        r[j][i]=a*x[j][i]+b;
    }}
}

void M_scaling(int nx1, int nx2,int ny1, int ny2, float **x, float a){
    /* x=x/a */
    int i, j;

    for (j=ny1;j<=ny2;j++){
    for (i=nx1;i<=nx2;i++){
        x[j][i]=x[j][i]/a;
    }}

}

void M_sqrt(int nx1, int nx2,int ny1, int ny2, float **x, float **r){
    /* r=sqrt(x) */
    int i, j;

    for (j=ny1;j<=ny2;j++){
    for (i=nx1;i<=nx2;i++){
        r[j][i]=sqrt(x[j][i]);
    }}
}

void M_multiply_sqrt(int nx1, int nx2,int ny1, int ny2, float **x, float **y, float **r){
    /* r=sqrt(x*y) */
    int i, j;

    for (j=ny1;j<=ny2;j++){
    for (i=nx1;i<=nx2;i++){
        r[j][i]=sqrt(x[j][i]*y[j][i]);
    }}
}

void M_multiply_a(int nx1, int nx2,int ny1, int ny2,float **a, float **b,float **cor){
    
    /* Local varibales */
    int i, j;
    
    for (j=ny1;j<=ny2;j++){
    for (i=nx1;i<=nx2;i++){
        cor[j][i] +=a[j][i]*b[j][i];
    }}
}

void M_multiply_a2(int nx1, int nx2,int ny1, int ny2, float **a, float **b, float **c, 
                   float **d, float **cor1, float **cor2){
    /* Local varibales */
    int i, j;
    
    for (j=ny1;j<=ny2;j++){
    for (i=nx1;i<=nx2;i++){
        cor1[j][i] +=(a[j][i]+c[j][i])*(b[j][i]+d[j][i]);
        cor2[j][i] +=(a[j][i]-c[j][i])*(b[j][i]-d[j][i]);
    }}
}
