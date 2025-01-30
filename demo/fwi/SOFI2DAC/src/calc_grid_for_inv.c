
#include "fd.h"
//call this function if(FORWARD_ONLY==0)
/* save every IDX and IDY spatial point during the forward modelling */
void calc_grid_for_inv(void){

    extern int MYID, NXG, NYG, NX, NY, NT, POS[3], IDX, IDY;
    extern int NX1, NX2, NY1, NY2, NXI, NYI;
//     extern float TSNAPINC;
    /* local variables */
    int i,ii;
    
    /* save every IDX and IDY spatial point during the forward modelling */
    //IDX=2;
    //IDY=2;
    //NTDTINV=ceil((float)NT/(float)(DTINV));		/* round towards next higher integer value */
//     NTDTINV=(NT-1)/DTINV+1;
//     IDT=iround(TSNAPINC/DT);
    
    for(i=1;i<=NX;i++){
        ii=i+POS[1]*NX;
        if(ii%IDX==1){
            NX1=i;
            break;
        }
    }
    for(i=NX;i>=1;i--){
        ii=i+POS[1]*NX;
        if(ii%IDX==1){
            NX2=i;
            break;
        }
    }
    for(i=1;i<=NY;i++){
        ii=i+POS[2]*NY;
        if(ii%IDY==1){
            NY1=i;
            break;
        }
    }
    for(i=NY;i>=1;i--){
        ii=i+POS[2]*NY;
        if(ii%IDY==1){
            NY2=i;
            break;
        }
    }
    if(IDX==1){
        NX1=1;
        NX2=NX;}
    if(IDY==1){
        NY1=1;
        NY2=NY;}
    
    NXI=(NX2-NX1)/IDX+1;
    NYI=(NY2-NY1)/IDY+1;
//     NXGI=(NXG-1)/IDX+1;
//     NYGI=(NYG-1)/IDY+1;
    if(MYID==0){
        printf("NXGI, NYGI=%d %d\n", (NXG-1)/IDX+1, (NYG-1)/IDY+1);
    }
    
//     MPI_Barrier(MPI_COMM_WORLD);
    printf("MYID=%d,POS=%d %d, NX1=%d,NX2=%d,NY1=%d,NY2=%d\n",MYID,POS[1],POS[2],NX1,NX2,NY1,NY2);
    MPI_Barrier(MPI_COMM_WORLD);

}
