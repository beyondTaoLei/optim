#include "fd.h"
void extract_current_sources(char srcname[256], st_acquisition *acq){
    
    extern int MYID, IENDX, IENDY,NT, NXG, NYG;
    extern int FREE_SURF, FW, FDORDER, FDORDER_TIME, POS[3];
    extern float DT, DH, REFREC[4];
    
    /* local variables */
    int nsrc_glob, nt, nsrc_out, nsrc_loc;
    float dt, Sx, Sy;
    int **srcpos_all=NULL, *repeated_tmp=NULL;
    float **signals_all=NULL, **srcpos_tmp3d=NULL;
    
    /*1. read source from su-format file*/
    //Note: we will add the signal at receiver position tr.gx tr.gy
    if(MYID==0){
        srcpos_tmp3d = read_su_header_all(srcname, &nsrc_glob, &nt, &dt, &Sx, &Sy);
        signals_all = read_su(srcname, &nsrc_glob, &nt);
        /*2. map position to gridpoint*/
        srcpos_all=coord_grid(DH, REFREC, srcpos_tmp3d, &repeated_tmp, nsrc_glob, &nsrc_out);
        if(nsrc_glob!=nsrc_out||nt!=NT||dt!=DT){
            if(nsrc_glob>nsrc_out)  free_ivector(repeated_tmp,1,nsrc_glob-nsrc_out);
            printf("%d %d %d %d %f %f\n",nsrc_glob, nsrc_out,nt,NT,dt,DT);
            printf("Please make sure the source positions are non-redundant, NT and DT are correct.\n");
            exit(0);
        }
        free_matrix(srcpos_tmp3d,1,3,1,nsrc_glob);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&nsrc_glob,1,MPI_INT,0,MPI_COMM_WORLD);
    if(MYID!=0){
        srcpos_all = imatrix(1,3,1,nsrc_glob);
        signals_all = matrix(1,nsrc_glob,1,NT);
    }
    MPI_Bcast(&Sx,1,MPI_FLOAT,0,MPI_COMM_WORLD);
    MPI_Bcast(&Sy,1,MPI_FLOAT,0,MPI_COMM_WORLD);
    MPI_Bcast(&srcpos_all[1][1],nsrc_glob*3,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&signals_all[1][1],nsrc_glob*NT,MPI_FLOAT,0,MPI_COMM_WORLD);
    /*3. generate local positions and local signals*/
    /* find this single source positions on subdomains */ 
    acq->srcpos_loc = splitpos(IENDX, IENDY, POS, srcpos_all, &nsrc_loc, nsrc_glob);
    acq->nsrc_loc=nsrc_loc;
    acq->Sx=Sx;
    acq->Sy=Sy;
    acq->signals_loc = matrix(1,nsrc_loc,1,NT);
    seis_glob2loc(nsrc_loc, NT, signals_all, acq->signals_loc, acq->srcpos_loc);
    free_imatrix(srcpos_all,1,3,1,nsrc_glob);
    free_matrix(signals_all,1,nsrc_glob,1,NT);
    MPI_Barrier(MPI_COMM_WORLD);
}
    
    
    
    
    
    