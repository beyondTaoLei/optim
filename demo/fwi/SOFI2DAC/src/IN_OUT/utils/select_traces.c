#include "par.h"
#include "comm.h"
/*
  gcc -o select_traces select_traces.c -I../include -L../lib -L/home/tao/seis_software/SAC/lib -lIN_OUT -lsac -lsacio -lm -lfftw3
  ./select_traces in=data/Tdiff.bin out=data/out2.dat s=/home/tao/seis_software/FWI_PRJ_Latest/snap/SR/S_POS.dat r=/home/tao/seis_software/FWI_PRJ_Latest/snap/SR/R_POS.dat tc=50 perc=0.1
 */
/* Global variables */
int xargc; 
char **xargv;

int main(int argc, char *argv[]){
    
    /* IN and OUT */
    int verbose;
    float tc, perc;
    char *in, *out, *s, *r;
    
    /* Define variables */
    int nsrc, ntr;
    int isrc, itr, index, i, j;
    float **mat_s=NULL, **mat_r=NULL, **data=NULL, **data2=NULL;
    FILE *fp;
    
    /* initialize getpar */
    initargs(argc,argv);
    //verbose=0; getparint("verbose",&verbose);
    if(!getparstring("in", &in)){
        printf("command as: ./select_traces in=data/Tdiff.bin out=data/out2.dat s=/home/tao/seis_software/FWI_PRJ_Latest/snap/SR/S_POS.dat r=/home/tao/seis_software/FWI_PRJ_Latest/snap/SR/R_POS.dat tc=50 perc=0.1\n");
        throw_error("parameter [in] must be specified!\n");}
    if(!getparstring("s", &s)) throw_error("parameter [s] must be specified!\n");
    if(!getparstring("r", &r)) throw_error("parameter [r] must be specified!\n");
    if(!getparstring("out", &out))  out= "out.dat";
    if(!getparfloat("tc",&tc)) throw_error("tc must be specified!\n");
    if(!getparfloat("perc",&perc)) throw_error("perc must be specified!\n");
    
    mat_s=read_text(0, 2, &nsrc, s);
    mat_r=read_text(0, 2, &ntr, r);
    data = matrix(1, nsrc, 1, ntr);
    data2 = matrix(1, 6, 1, nsrc*ntr);
    V_read(in,&data[1][1],nsrc*ntr);
    
    printf("ntr=%d nsrc=%d \n",ntr, nsrc);
    

    
    for(isrc=1;isrc<=nsrc;isrc++){
        for(itr=1;itr<=ntr;itr++){
            index=(isrc-1)*ntr+itr;
            /* global index no. */
            data2[1][index]= isrc;
            data2[2][index]= itr;
            data2[3][index]= mat_s[isrc][2];
            data2[4][index]= mat_r[itr][2];
            data2[5][index]= data[isrc][itr];
            if(fabs(data2[5][index])< tc*perc) data2[6][index]=1.;
            //printf("%d ", (int)data2[index][6]);
        }
    }
    for(i=1;i<=nsrc*ntr;i++){
        if(data2[6][i]==0.){
            for(j=i+1;j<=nsrc*ntr;j++){
                if(data2[3][i] == data2[4][j]&&data2[4][i] == data2[3][j]){
                    data2[6][j]=1.;
                    break;
                }
            }
        }
    }
    
    fp=fopen(out,"w+");
    for(index=1;index<=nsrc*ntr;index++){
        //printf("%d ", (int)data2[6][index]);
        if((int)data2[6][index]==0) 
            fprintf(fp, "%8d %8d %.4e\n", (int)data2[1][index],(int)data2[2][index],data2[5][index]);
    }
    fclose(fp);
    
//     out="data/out4.bin";
//     V_write(out,&data2[6][1],nsrc*ntr);

    free_matrix(mat_s, 1, nsrc, 1, 2);
    free_matrix(mat_r, 1, ntr, 1, 2);
    free_matrix(data, 1, nsrc, 1, ntr);
    free_matrix(data2, 1, 6, 1, nsrc*ntr);
    
    return 0;
}