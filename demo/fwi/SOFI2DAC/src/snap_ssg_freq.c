#include "fd.h"
float *read_freqs(char filename[256], int *nrow){
    int c;
    int nr, irow;
    float val, *vec = NULL;
    FILE *fp;
    
    fp=fopen(filename,"r");
    if (fp==NULL){
        printf(" File %s could not be opened !",filename);
        exit(0);
    }
    nr=0;
    /* The last line should be empty!!! */
    while ((c=fgetc(fp)) != EOF)
        if (c=='\n') ++(nr);
    rewind(fp);
    
    vec=vector(1, nr);
    for(irow=1;irow<=nr;irow++){
        fscanf(fp,"%f",&val);
        vec[irow]=val;
    }
    
    fclose(fp);
    printf(" Number of data found: %i\n",nr);

    *nrow= nr;
    
    return vec;
}
void allocate_mem_wfdf(int flag, st_wavefildfreq *wfd_f){
    
    extern int INV_RHO, NX, NY, NXI, NYI, FDORDER;
    extern char FREQ_FILE[STRING_SIZE];
    int nd, num_freqs;
    int i,j,k,l;
    nd = FDORDER/2;
    if(flag==1){/* allocate memory */
        wfd_f->freqs=read_freqs(FREQ_FILE, &num_freqs);
        wfd_f->num_freqs=num_freqs;
        wfd_f->psp_f =  f4tensor(1,num_freqs,1,2,1,NYI,1,NXI);
        for(l=1; l<=num_freqs; l++){
        for(k=1;k<=2;k++){
            for(j=1;j<=NYI;j++){
            for(i=1;i<=NXI;i++){
                wfd_f->psp_f[l][k][j][i]=0.0;
            }}
        }}
        if(INV_RHO==1){
            wfd_f->pvx_f =  f4tensor(1,num_freqs,1,2,1,NYI,1,NXI);
            wfd_f->pvy_f =  f4tensor(1,num_freqs,1,2,1,NYI,1,NXI);
            for(l=1; l<=num_freqs; l++){
            for(k=1;k<=2;k++){
                for(j=1;j<=NYI;j++){
                for(i=1;i<=NXI;i++){
                    wfd_f->pvx_f[l][k][j][i]=0.0;
                    wfd_f->pvy_f[l][k][j][i]=0.0;
                }}
            }}
        }
    }else{ /* deallocate memory */
        num_freqs=wfd_f->num_freqs;
        free_vector(wfd_f->freqs, 1, num_freqs);
        free_f4tensor(wfd_f->psp_f,1,num_freqs,1,2,1,NYI,1,NXI);
        if(INV_RHO==1){
            free_f4tensor(wfd_f->pvx_f,1,num_freqs,1,2,1,NYI,1,NXI);
            free_f4tensor(wfd_f->pvy_f,1,num_freqs,1,2,1,NYI,1,NXI);
        }
    }
}

int *generate_frequencies(int *num_freqs){
    
    extern int OMEGA1, OMEGA2, OMEGA_INC;
    int i, num;
    int *freqs =NULL;
    
    num=(int)((OMEGA2-OMEGA1)/OMEGA_INC+1);
    freqs=ivector(1, num);
    for(i=1;i<=num;i++) freqs[i]=OMEGA1+OMEGA_INC*(i-1);
    *num_freqs=num;
    return freqs;
}

void dft_forward(int nt, int num_freqs, int *freqs, float **vx,float ****Fvx){
/* modified based on discfourier.c in IFOS3D */
/*-------------------------------------------------------------------------
 * calculation of a discrete Fourier transfomation on the fly:
 * Fouriercomponents of forward or backpropagated wavefields are summed up for each frequency
 * S. Butezr 2013
 --------------------------------------------------------------------------*/
    /*
     INPUT:
     nt: from 1 to NT;
     num_freqs: the number of the frequencies of interest;
     freqs: the array of the frequency indexes, which element can be 
            from 1 to NT;
     vx: the 2D array of wavefields;
     Fvx: the wavefields in frequency domain.
     */
    extern int NX, NY, NT;
    int i, j, l;
    double alpha, trig1, trig2;

    for(l=1; l<=num_freqs; l++){
        /* for sig[1:NT] nt:[1:NT], jw:[1:NT], 
         * however you can specify any one frequency;*/
        alpha=2.*PI*(freqs[l]-1.)*(nt-1.)/NT;
        trig1=cos(alpha);
        trig2=-sin(alpha);
        
        for (j=1;j<=NY;j++){
        for (i=1;i<=NX;i++){
            Fvx[l][1][j][i]+=vx[j][i]*trig1;
            Fvx[l][2][j][i]+=vx[j][i]*trig2;
        }}
    }				
}

void dft_forward2(int nt, int num_freqs, float *freqs, float **vx,float ****Fvx){
/* modified based on discfourier.c in IFOS3D */
/*-------------------------------------------------------------------------
 * calculation of a discrete Fourier transfomation on the fly:
 * Fouriercomponents of forward or backpropagated wavefields are summed up for each frequency
 * S. Butezr 2013
 --------------------------------------------------------------------------*/
    /*
     INPUT:
     nt: from 1 to NT;
     num_freqs: the number of the frequencies of interest;
     freqs: the array of the frequency indexes, which element can be 
            from 1 to NT;
     vx: the 2D array of wavefields;
     Fvx: the wavefields in frequency domain.
     */
    extern int NX1, NX2, NY1, NY2, IDX, IDY, NT;
    extern float DT;
    int i, j, l, ii, jj;
    double alpha, trig1, trig2;

    for(l=1; l<=num_freqs; l++){
        /* for sig[1:NT] nt:[1:NT], jw:[1:NT], 
         * however you can specify any one frequency;*/
        alpha=2.*PI*freqs[l]*(nt-1.)*DT;
        trig1=cos(alpha);
        trig2=-sin(alpha);
        
        for(j=NY1;j<=NY2;j=j+IDY){
            jj=(j-NY1)/IDY+1;
            for(i=NX1;i<=NX2;i=i+IDX){
                ii=(i-NX1)/IDX+1;
                Fvx[l][1][jj][ii]+=vx[j][i]*trig1;
                Fvx[l][2][jj][ii]+=vx[j][i]*trig2;
        }}
    }				
}

void save_wavefields_freq(FILE *fp, st_wavefildfreq *wfd_f, int ishot){
    //SNAP=7 : values in vx, vy, and pressure field in the frequency domain
    extern int NX, NY, MYID, POS[3], IDX, IDY, INV_RHO;
    extern char SNAP_FILE[STRING_SIZE];
    
    int i,j, l, numx, numy, num;
    char snapfile_x[STRING_SIZE], snapfile_y[STRING_SIZE];
    char snapfile_p[STRING_SIZE];
    char snapfile_x_i[STRING_SIZE], snapfile_y_i[STRING_SIZE];
    char snapfile_p_i[STRING_SIZE];
    int num_freqs=wfd_f->num_freqs;
    float *freqs=wfd_f->freqs;
    float ****vx_f=wfd_f->pvx_f;
    float ****vy_f=wfd_f->pvy_f;
    float ****sp_f=wfd_f->psp_f;
    FILE *fp_vx, *fp_vy, *fp_p, *fp_vx_i, *fp_vy_i, *fp_p_i;
    
    numx=(NX-1)/IDX+1;
    numy=(NY-1)/IDY+1;
    num=numx*numy;
    
    sprintf(snapfile_p,  "%s.%s.shot%i.%i.%i.bin",SNAP_FILE,"sp",  ishot,POS[1],POS[2]);
    sprintf(snapfile_p_i,"%s.%s.shot%i.%i.%i.bin",SNAP_FILE,"sp.i",ishot,POS[1],POS[2]);
    fprintf(fp,"%s\n\n",snapfile_p);
    fp_p =fopen(snapfile_p,"w");
    fp_p_i =fopen(snapfile_p_i,"w");
    for(l=1; l<=num_freqs; l++){
        fseek(fp_p,  (l-1)*num*sizeof(float),SEEK_SET);
        fseek(fp_p_i,  (l-1)*num*sizeof(float),SEEK_SET);
        for (i=1;i<=NX;i+=IDX){
        for (j=1;j<=NY;j+=IDY){
            writedsk(fp_p,   sp_f[l][1][j][i],3);
            writedsk(fp_p_i, sp_f[l][2][j][i],3);
        }}
    }
    fclose(fp_p_i);
    fclose(fp_p);
    
    if(INV_RHO==1){
        sprintf(snapfile_x,  "%s.%s.shot%i.%i.%i.bin",SNAP_FILE,"vx",  ishot,POS[1],POS[2]);
        sprintf(snapfile_y,  "%s.%s.shot%i.%i.%i.bin",SNAP_FILE,"vy",  ishot,POS[1],POS[2]);
        sprintf(snapfile_x_i,"%s.%s.shot%i.%i.%i.bin",SNAP_FILE,"vx.i",ishot,POS[1],POS[2]);
        sprintf(snapfile_y_i,"%s.%s.shot%i.%i.%i.bin",SNAP_FILE,"vy.i",ishot,POS[1],POS[2]);
        fprintf(fp,"%s\n", snapfile_x);
        fprintf(fp,"%s\n", snapfile_y);
        fp_vx=fopen(snapfile_x,"w");
        fp_vy=fopen(snapfile_y,"w");
        fp_vx_i=fopen(snapfile_x_i,"w");
        fp_vy_i=fopen(snapfile_y_i,"w");
        
        for(l=1; l<=num_freqs; l++){
            fseek(fp_vx, (l-1)*num*sizeof(float),SEEK_SET);
            fseek(fp_vy, (l-1)*num*sizeof(float),SEEK_SET);
            fseek(fp_vx_i, (l-1)*num*sizeof(float),SEEK_SET);
            fseek(fp_vy_i, (l-1)*num*sizeof(float),SEEK_SET);
            for (i=1;i<=NX;i+=IDX){
            for (j=1;j<=NY;j+=IDY){
                writedsk(fp_vx,  vx_f[l][1][j][i],3);
                writedsk(fp_vy,  vy_f[l][1][j][i],3);
                writedsk(fp_vx_i,vx_f[l][2][j][i],3);
                writedsk(fp_vy_i,vy_f[l][2][j][i],3);
            }}
        }
        fclose(fp_vx);
        fclose(fp_vy);
        fclose(fp_vx_i);
        fclose(fp_vy_i);
    }
}