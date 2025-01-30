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

/*  quickSort
 * This public-domain C implementation by Darel Rex Finley .
 * http://alienryderflex.com/quicksort */
void quicksort1(float *arr, int start_index, int elements) {
    /*
     arr is a vector;
     start_index means the start of index;
     elements: arr[elements-1] means the last element.
     */
    /* number of levels: 64 is enough for 1.8e19 elements to be sorted */
    #define  MAX_LEVELS  64
    
    int	beg[MAX_LEVELS], end[MAX_LEVELS], i=0, L, R, swap ;
    float	piv;
    
    beg[0]=start_index;
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

int find_min_positive_num1(int n1,int n2,float *a){
    /*find the element position in which is the minimum positive value  */
    int mid,left,right;
    
    if((a[n1]<0.0)||(a[n2]<=0.0)){
        printf("element value of array should be more than 0\n");
        exit(0);
    }
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

/* Return percentiles of the values in a matrix for the percentages 
 * p1 and p2 in the interval [0,100]. */
float global_percentile(float *fullvec, int n1, int n2, float p1) {
    int num, ij1, i, j;
    float prctile;
    float *fullvec2=NULL;
    
    if(p1<0.0||p1>1.0){
        printf("p1: %f from function global_percentile\n",p1);
        printf("p1 should be in (0.0 1.0)\n");
        exit(0);
    }
    
    num=n2-n1+1;
    //fullvec2 = vector(1,num);
    
    fullvec2=(float *)malloc((num+1)*sizeof(float));
    fullvec2[0]=0.0;
    
    j=1;
    for(i=n1;i<=n2;i++){
        fullvec2[j++]=fullvec[i];
    }
    
    
    /*
     * quicksort1(float *arr, int start_index, int elements) {
        arr is a vector;
        start_index means the start of index;
        elements: arr[elements-1] means the last element.
     */
    quicksort1(fullvec2,1,num+1);
    //////////
    ij1=find_min_positive_num1(1,num,fullvec2);
    ij1=ij1-1+(int)(num*p1);
    //ij1=(int)(num*p1);
    //////////
    prctile=fullvec2[ij1];
    //free_vector(fullvec2, 1,num);
    free(fullvec2);
    //printf("prctile=%e \n",prctile);
    return prctile;
   
}

int get_file_size(char filename[256]){
    /* return the file size by bytes */
    
    int size;
    FILE *fp;
    
    fp = fopen(filename,"r");
    if (fp==NULL) {
        printf(" Can't read file: %s\n",filename);
        exit(1);
    }
    fseek(fp,0L,SEEK_END);
    size=ftell(fp);
    fclose(fp);
    
    return size;
}

void V_write(char modfile[256],float *x_vector,int nn){
    /*for example: write_vector(modfile,x_global[offset],NXG*NYG); */
    FILE *fp;
    
    fp=fopen(modfile,"w");
    fwrite(x_vector,sizeof(float),nn,fp);
    fclose(fp);
}

void V_write_txt(char modfile[256],float *x_vector,int nn){
    /*for example: write_vector(modfile,x_global[offset],NXG*NYG); */
    int i;
    FILE *fp;
    
    fp=fopen(modfile,"w+");
    for (i=0;i<nn;i++){
        fprintf(fp,"%15e\n",x_vector[i]);
    }
    fclose(fp);
}

void V_read(char modfile[256],float *x_vector,int nn){
    /*for example: read_vector(modfile,x_global[offset],NXG*NYG); */
    FILE *fp;
    
    fp=fopen(modfile,"r");
    fread(x_vector,sizeof(float),nn,fp);
    fclose(fp);
}

float V_max_abs(int n1, int n2, float * x){
    int i;
    float max_x;

    max_x=0.0;
    for (i=n1;i<=n2;i++){
        if(fabs(x[i])>max_x)    max_x=fabs(x[i]);
    }
    
    return max_x;
}

int V_max_index(int n1, int n2, float * x){
    /* finds the index in which the element is maximum of whole array. */
    int i, index=n1;
    float max_x;

    max_x=x[n1];
    for (i=n1+1;i<=n2;i++){
        if(x[i]>=max_x){
            max_x=x[i];
            index=i;
        }
    }
    
    return index;
}

int V_max_abs_index(int n1, int n2, float * x){
    /* finds the index in which the element is maximum of whole array. */
    int i, index=n1;
    float max_x;

    max_x=fabs(x[n1]);
    for (i=n1+1;i<=n2;i++){
        if(fabs(x[i])>=max_x){
            max_x=fabs(x[i]);
            index=i;
        }
    }
    
    return index;
}

void V_normalize(int n1, int n2, float * x){
    /* normalize vector with its absolute maximum */
    int i;
    float max1;
    
    max1 = V_max_abs(n1, n2, x);
    if(max1!=0.0){
        for (i=n1;i<=n2;i++){
            x[i]=x[i]/max1;
        } 
    }
}

float V_normL2(int n1, int n2, float * x){
    /* The routine V_normL2 returns the Euclidian
       norm of a vector x of size (n2-n1+1)   */
    int i;
    float norm_x;

    norm_x=0.0;
    for (i=n1;i<=n2;i++){
        norm_x=norm_x+x[i]*x[i];
    }
    norm_x=sqrt(norm_x);
    
    return norm_x;
    
}

float V_dotproduct(int n1, int n2, float *x, float *y){
    /* calculates the dot product between two given vectors */
    int i;
    float sum1;

    sum1=0.0;
    for (i=n1;i<=n2;i++){
        sum1 +=x[i]*y[i];
    }
    
    return sum1;
}

/**********************************/
/**********************************/
/* operate vectors element-wise.  */
/**********************************/
/**********************************/
void V_init(int n1, int n2, float *x){
    /* x=0.0 */
    int i;

    for (i=n1;i<=n2;i++){
        x[i]=0.0;
    }
}

void V_copy(int n1, int n2, float *x, float *r){
    /* r=x */
    int i;

    for (i=n1;i<=n2;i++){
        r[i]=x[i];
    }
}

void V_pick(int n1, int n2, float *x, float *r){
    /* r[0:(n2-n1)]=x[n1:n2] */
    int i, j=0;
    
    for (i=n1;i<=n2;i++){
        r[j++]=x[i];
    }
}

void V_abs(int n1, int n2, float *x, float *r){
    /* r=abs(x) */
    int i;

    for (i=n1;i<=n2;i++){
        r[i]=fabs(x[i]);
    }
}

void V_add(int n1, int n2, float *x, float *y,float *r){
    /* r=x+y */
    int i;

    for (i=n1;i<=n2;i++){
        r[i]=x[i]+y[i];
    }
}

void V_subtract(int n1, int n2, float *x, float *y,float *r){
    /* r=x-y */
    int i;

    for (i=n1;i<=n2;i++){
        r[i]=x[i]-y[i];
    }
}

void V_multiply(int n1, int n2, float *x, float *y,float *r){
    /* r=x*y */
    int i;

    for (i=n1;i<=n2;i++){
        r[i]=x[i]*y[i];
    }
}

void V_divide(int n1, int n2, float *x, float *y,float *r){
    /* r=x/y */
    int i;

    for (i=n1;i<=n2;i++){
        r[i]=x[i]/y[i];
    }
}

void V_square_sum(int n1, int n2, float *x, float *y,float *r){
    /* r=x*x+y*y */
    int i;

    for (i=n1;i<=n2;i++){
        r[i]=x[i]*x[i]+y[i]*y[i];
    }
}

void V_square_sum_a(int n1, int n2, float *x, float *y,float *r){
    /* r+=x*x+y*y */
    int i;

    for (i=n1;i<=n2;i++){
        r[i]+=(x[i]*x[i]+y[i]*y[i]);
    }
}

void V_AXPY(int n1, int n2, float a, float * x, float * y,float *r){
    /* r=a*x+y */
    int i;

    for (i=n1;i<=n2;i++){
        r[i]=a*x[i]+y[i];
    }
}

void V_AXPb(int n1, int n2, float a, float * x,  float b, float *r){
    /* r=ax+b */
    int i;

    for (i=n1;i<=n2;i++){
        r[i]=a*x[i]+b;
    }
}

void V_scaling(int n1, int n2, float *x, float a){
    /* x=x/a */
    int i;

    for (i=n1;i<=n2;i++){
        x[i]=x[i]/a;
    }

}

void V_sqrt(int n1, int n2, float *x, float *r){
    /* r=sqrt(x) */
    int i;

    for (i=n1;i<=n2;i++){
        r[i]=sqrt(x[i]);
    }
}

void V_multiply_sqrt(int n1, int n2, float *x, float *y, float *r){
    /* r=sqrt(x*y) */
    int i;

    for (i=n1;i<=n2;i++){
        r[i]=sqrt(x[i]*y[i]);
    }
}

void V_multiply_a(int n1, int n2,float *a, float *b,float *cor){
    
    /* Local varibales */
    int i;
    
    for (i=n1;i<=n2;i++){
        cor[i] +=a[i]*b[i];
    }
}

void V_multiply_a2(int n1, int n2, float *a, float *b, float *c, 
                   float *d, float *cor1, float *cor2){
    /* Local varibales */
    int i;
    
    for (i=n1;i<=n2;i++){
        cor1[i] +=(a[i]+c[i])*(b[i]+d[i]);
        cor2[i] +=(a[i]-c[i])*(b[i]-d[i]);
    }
}

void V_grad(int n1, int n2, float *x, float *r){
    /* the gradient of an N-dimensional array */
    int i;

    for (i=n1+1;i<n2;i++){
        r[i] =0.5*(x[i+1]-x[i-1]);
    }
    r[n1]=x[n1+1]-x[n1];
    r[n2]=x[n2]-x[n2-1];
}

void V_integ(int n1, int n2, float *x, float *r){
    /* scaled integration */
    int i;
    float ints=0.0;

    for (i=n1;i<=n2;i++){
        ints += x[i];
        r[i] = ints;
    }
}

void V_shift(int n1, int n2, int lag, float *x, float *r){
    /* shift the array along the axis i-->i+lag*/
    int i;
    if(lag>0){ 
        //delay
        //for example: lag =2, then 1 2 3 4 5-->0 0 1 2 3
        for (i=n1;i<=n1+lag-1;i++){
            r[i] = 0.0;
        }
        for (i=n1+lag;i<=n2;i++){
            r[i] = x[i-lag];
        }
        
    }else{/* lag<=0 */
        //advance
        //for example: lag =-2, then 1 2 3 4 5-->3 4 5 0 0
        for (i=n1;i<=n2+lag;i++){
            r[i] = x[i-lag];
        }
        for (i=n2+lag+1;i<=n2;i++){
            r[i] = 0.0;
        }
    }
}

void V_reversal(int n1, int n2, float *x, float *r){
    /* reverse the array along the axis */
    //for example: 1 2 3 4 5--> 5 4 3 2 1
    int i, j;
    
    j=n2;
    for (i=n1; i<=n2; i++){
        r[j]=x[i];
        j--;
    }
}


/*
 *  gcc vector_operations.c -o vector_operations -lm && ./vector_operations
 */

// int main(int argc, char **argv){
//    
//     /* Local varibales */
//     int i, j,n1, n2;
//     float a[10],b[10],c[10],d[10],e[10],f[10],tmp;
//     float sub[5];
//     n1=0;
//     n2=9;
//     
//     
//     for (i=n1;i<=n2;i++){
//         a[i]=i*1.0;
//         b[i]=i*2.0;
//         c[i]=i*3.0;
//         d[i]=i*4.0;
//         e[i]=i*5.0;
//         f[i]=i*6.0;
//     }
    //a[1]=-1;a[4]=-20;
    
    //V_copy(n1, n2, a, c);
    //V_add(n1, n2, a, b, c);
    //V_subtract(n1, n2, a, b, c);
    //V_multiply(n1, n2, a, b, c);
    //V_divide(n1, n2, a, b, a);
    //V_square_sum(n1, n2, a, b, c);
    //V_square_sum_a(n1, n2, a, b, c);
//     V_AXPY(n1, n2, 5.0, a, b, c);
    //V_AXPb(n1, n2, 10.0, a,1, c);
    
    //V_scaling(n1, n2, a,10.0);
    //V_sqrt(n1, n2, a,c);
    //V_multiply_sqrt(n1, n2, a, b, c);
    //V_multiply_a(n1, n2, a, b, c);
    /*V_multiply_a2(n1, n2, a, b, c,d,e,f);*/
//     V_pick(3, 7, a, sub);
//     printf("\n");
//     for (i=n1;i<=n2;i++) printf("%5.1f ",c[i]);
//     for (i=n1;i<=n2;i++) printf("%5.1f ",a[i]);
//     printf("\n");
//     for (i=0;i<=4;i++) printf("%5.1f ",sub[i]);
//     printf("\n");
    
    //V_write("1.bin",e,10);
    
    //V_read("/home/tao/seis_software/proj01/precond/approx_hessian_it1_shot1.bin",a,1700*350);
    //tmp=global_percentile(a, n1, n2, 0.03);
    //printf("%e \n",tmp);
    //printf("%f \n",V_normL2(n1, n2, a));
    //i=find_min_positive_num1(n1,n2,a);
    //printf("%d \n",i);
    //for (i=n1;i<=n2;i++) printf("%5.1f ",f[i]);
    
// }

