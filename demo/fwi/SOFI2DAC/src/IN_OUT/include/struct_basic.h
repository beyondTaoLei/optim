
typedef struct Source{
    /* the values from json file */
    /*flag_su_bin:
    1 su 
    2 bin 
    3 for FWI case, set geo_srci with 2, and set geo_srca with 1
    
    src_snap_type:
    1: boundary saving
    2: the specified depth
    3: the specified coordinates
    4: the whole model space
    
    add_replace:
    1: tmp = compnt+b*sig
    2: tmp = b*sig;
    3: tmp = sig;
    0 other: tmp = compnt;
    return tmp;
    */
    int flag_su_bin, src_snap_type, read_last_snap;
    int add_replace_vx, add_replace_vy, add_replace_sp, add_replace_snap;
    int add_replace_sxx, add_replace_sxy, add_replace_syy;
    float memory_depth;
    char src_f_vx[STRING_SIZE], src_f_vy[STRING_SIZE];
    char src_f_sp[STRING_SIZE], src_f_snap[STRING_SIZE];
    char snap_coords_file[STRING_SIZE];
    
    //int   **srcpos
    int **srcpos_loc;
    int nsrc_glob, nsrc_loc;
    
    /* the file of saving the snapshots in time order */
    FILE *fp_vx,  *fp_vy, *fp_sp;
    FILE *fp_sxx, *fp_sxy, *fp_syy;
    float *snap_all; /* [nsrc] */
    
    /*just avariable in several node*/  
    /*the dimension of signal is [nsrc_loc][nt] */
    float **signal_vx,  **signal_vy,  **signal_sp;
    float **signal_sxx, **signal_sxy, **signal_syy;
    
    /*the dimension of snap is [nsrc_loc] */
    float *snap_vx,  *snap_vy,  *snap_sp;
    float *snap_sxx, *snap_sxy, *snap_syy;
    
}st_source;

typedef struct Boundary{
    //int   **srcpos
    int **pos_loc;
    int nsnap_glob, nsnap_loc;
    
    /* the file of saving the snapshots in time order */
    FILE *fp_vx,  *fp_vy, *fp_sp;
    FILE *fp_sxx, *fp_sxy, *fp_syy;

}st_bdr;

typedef struct Wavefield{
    
    int ng_loc;
    /* the wavefield compenents */
    float *fvx, *fvy, *fsxx, *fsyy, *fsxy, *fsp;
    float *evx, *evy;

}st_wfd;

typedef struct Gradient{
    
    /* flags which decide whether the parameter will be inversed */
    int inv_vp, inv_vs, inv_rho;
    int ng_loc, ng_glob;
    /* the coordinates of inversed grids */
    int *winx, *winy;
    int **win_loc;/* i/x: [1][ij], j/y: [2][ij] */
    /* model parameters */
    float *Mvp, *Mvs, *Mrho;
    
    /* zero-lag cross correlation */
    float *cvx, *cvy, *csxx, *csyy, *csxy, *csp;
    /* the gradient compenents */
    float *grad_lam, *grad_mu, *grad_rho_s;
    float *grad_vp, *grad_vs, *grad_rho;
    
    /* preconditioning */
    int eprecond, eprecond_per_shot;
    float epsilon, *Ws, *Wr, *We;

}st_grad;