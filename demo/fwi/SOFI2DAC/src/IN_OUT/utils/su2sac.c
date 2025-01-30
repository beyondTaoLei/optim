#include <sacio.h>
#include "par.h"
#include "source.h"

#define FALSE   0
#define TRUE    1
/* Global variables */
int xargc; 
char **xargv;

int main(int argc, char *argv[]){
    
    /* IN and OUT */
    char *fin, *fout;
    
    /* Define variables */
    int ntr, ns, itr, it;
    float dt, Sx, Sy;
    float xdummy[5];
    float **pos=NULL, **data=NULL;
    float *sig=NULL;
    char kevnm[256], kstnm[256];
    char fout2[256];

    float evla, evlo, stla, stlo;
    int nerr;
    float b, delta;
    int true  = TRUE;
    
    /* initialize getpar */
    initargs(argc,argv);
    if(!getparstring("fin", &fin)){
        printf("command as: ./su2sac fin=data/toy_example_p.su fout=out\n");
        throw_error("parameter [fin] must be specified!\n");}
    if (!getparstring("fout", &fout))  fout= "out";
    
    pos=read_su_header_all(fin, &ntr, &ns, &dt, &Sx, &Sy);
    data = read_su(fin, &ntr, &ns);
    sig = vector(1, ns);

    evla   = Sx;
    evlo   = Sy;
    b      = dt;
    delta  = dt;
    
    newhdr () ;
    strcpy(kevnm, "Event Name");
    setihv("IFTYPE", "ITIME", &nerr , strlen("IFTYPE"), strlen("ITIME"));
    setihv("IZTYPE", "IB",    &nerr , strlen("IZTYPE"), strlen("IB"));
    setkhv("KEVNM",  kevnm, &nerr, strlen("KEVNM"), SAC_STRING_LENGTH);
    setlhv("LEVEN",  &true,   &nerr , strlen("LEVEN"));

    setfhv("EVLA",   &evla,     &nerr, strlen("EVLA"));
    setfhv("EVLO",   &evlo,     &nerr, strlen("EVLO"));
    setfhv("B",      &b,      &nerr , strlen("B"));
    setfhv("DELTA",  &delta,  &nerr , strlen("DELTA")) ;
    setnhv("NPTS",   &ns,    &nerr, strlen("NPTS"));

    for(itr=1;itr<=ntr;itr++){
        for(it=1;it<=ns;it++){
            sig[it]= data[itr][it];
        }
        stla=pos[1][itr];
        stlo=pos[2][itr];
        
        /* update the file name for each trace */
        /* Define the file to be written    */
        if(ntr>1){
            sprintf(fout2, "%s_%d.sac",fout, itr);
        }else{
            sprintf(fout2, "%s.sac",fout);
        }
        sprintf(kstnm, "tr%d", itr);
        setkhv ( "KSTNM", kstnm, &nerr, strlen("KSTNM"), SAC_STRING_LENGTH);
        setfhv ( "STLA" , &stla , &nerr , strlen("STLA"));
        setfhv ( "STLO" , &stlo , &nerr , strlen("STLO"));
        wsac0(fout2, xdummy, &sig[1], &nerr, strlen(fout2));
    }
    
    /* deallocate memory */
    free_matrix(pos, 1, 3, 1, ntr);
    free_matrix(data, 1, ntr, 1, ns);
    free_vector(sig, 1, ns);
                    
    return 0;     
}