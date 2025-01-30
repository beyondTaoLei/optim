#include "par.h"
#include "comm.h"
/*
  gcc -o bin2mat2d bin2mat2d.c -I../include -L../lib -L/home/tao/seis_software/SAC/lib -lIN_OUT -lsac -lsacio -lm -lfftw3
  ./bin2mat2d n1=36 in=../bin/data/Tdiff.bin out=../bin/data/out.dat
 */
/* Global variables */
int xargc; 
char **xargv;

int main(int argc, char *argv[]){
    
    /* IN and OUT */
    int verbose, n1, n2, flag;
    float xmin, ymin, dh;
    char *in, *out;
    
    /* Define variables */
    int ix, iy, size;
    float x,y;
    float **data=NULL;
    float **sig=NULL;
    char fout[256];
    FILE *fp;
    
    /* initialize getpar */
    initargs(argc,argv);
    //verbose=0; getparint("verbose",&verbose);
    if(!getparstring("in", &in)){
        printf("command as: ./bin2mat2d n1=36 in=../bin/data/Tdiff.bin out=../bin/data/out.dat\n");
        throw_error("parameter [in] must be specified!\n");}
    if (!getparstring("out", &out))  out= "out.dat";
    
    if(!getparint("n1",&n1)) throw_error("parameter n1 must be specified!\n");
    if(!getparint("n2",&n2)){
        fp = fopen(in,"r");
        if (fp==NULL) {
            printf(" It Was not able to read BIN file: %s\n",in);
            exit(1);
        }
        fseek(fp,0L,SEEK_END);
        size=ftell(fp);
        n2=size/4/n1;
        fclose(fp);
    }
    printf("n1=%d n2=%d \n",n1, n2);
    /*output the true coordinates, one should set flag=1 */
    flag=0; getparint("flag",&flag);
    if(flag==1){
        xmin=0.0; getparfloat("xmin",&xmin);
        ymin=0.0; getparfloat("ymin",&ymin);
        dh=1.0; getparfloat("dh",&dh);
        printf("xmin=%f ymin=%f dh=%f \n",xmin, ymin, dh);
    }
    
    data = matrix(1, n2, 1, n1);
    V_read(in,&data[1][1],n2*n1);
    
    fp=fopen(out,"w+");
    if(flag==0){
        for(ix=1;ix<=n2;ix++){
            for(iy=1;iy<=n1;iy++){
                fprintf(fp, "%8d %8d %.4e\n",ix,iy,data[ix][iy]);
            }
        }
    }else{
        for(ix=1;ix<=n2;ix++){
            x=xmin+(ix-1)*dh;
            for(iy=1;iy<=n1;iy++){
                y=ymin+(iy-1)*dh;
                fprintf(fp, "%f %f %.4e\n",x,y,data[ix][iy]);
            }
        }
    }
    fclose(fp);

    
    free_matrix(data, 1, n2, 1, n1);

    return 0;
}