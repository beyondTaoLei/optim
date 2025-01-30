#include <sacio.h>
#include "comm.h"
/*
 gcc -c bin2sac.c -I/home/tao/seis_software/seisplatform/FD/SOFI2D/src/IN_OUT/include -I/home/tao/seis_software/SAC/include
 gcc -o bin2sac bin2sac.o -L/home/tao/seis_software/seisplatform/FD/SOFI2D/src/IN_OUT/lib -L/home/tao/seis_software/SAC/lib -lIN_OUT -lsac -lsacio -lm
 ./bin2sac ./aa/toy_example_vy.su.shot1 shot1/hh
 */


int main(int argc, char *argv[]){
  
    /* Define variables to be passed to wsac1() */
    int nerr;//max
    //float yfunc[ MAX ], x, beg, del;
    float beg, del;
    //char kname[ 10 ];
    
    int ns;
    float dt;
    float *sig=NULL;
    long size=0;
    char file_sac[256];
    char filein[256]; 
    FILE *fp;
    /////////////////////
    if(argc !=4){
        printf("The comamnd should be as: ./bin2sac 0.1 signal.bin signal.sac\n");
        exit(1);
    }else{
        sscanf(argv[1], "%f", &dt);
        sscanf(argv[2], "%s", filein);
        sscanf(argv[3], "%s", file_sac);
        
    }
    fp = fopen(filein,"r");
    if (fp==NULL) {
        printf(" It Was not able to read BIN file: %s\n",filein);
        exit(1);
    }
    fseek(fp,0L,SEEK_END);
    size=ftell(fp);
    ns=size/4;
    
    /* Define the beginning time, time sampling. */
    beg = dt;
    del = dt;
    
    sig= vector(1, ns);
    V_read(filein,&sig[1],ns);
    wsac1(file_sac, &sig[1], &ns, &beg, &del, &nerr, strlen(file_sac)) ;
    

    /* Write the SAC file kname
        - kname holds the name of the file to be written
        - yfunc Input Amplitude data
        - max number of points to be writtne
        - beg Beginning Time of the data
        - del Time Sampling of the series
        - nerr Error return Flag
        - strlen(kname) Length of the character array kname
    */
    //wsac1 (kname, yfunc, &max, &beg, &del, &nerr, strlen( kname )) ;

    /* Check the Error status
        - 0 on Success
        - Non-Zero on Error
    */
    if(nerr != 0) {
        fprintf(stderr, "Error writing SAC File: %s\n", file_sac);
        exit(-1);
    }

    exit(0);
}