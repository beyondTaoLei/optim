/* files to include */
#include "comm.h"
#include "segy.h"

int **coord_grid(float dh, float *refrec, float **coord, int **repeated_tr,
                 int ntr_prev, int *ntr);
float insert_src(float compnt, float sig, float b,int flag);
int **read_coordinates(int *ntr, float dh, float *refrec, char rec_file[256]);
void seis_glob2loc(int ntr,int ns, float **fulldate, float **section, int ** recpos_loc);
void seis_glob2loc_vec(int ntr, float *fulldate, float *section, int ** recpos_loc);
int source_subfix(char source_file[256]);
int **splitpos(int iendx, int iendy, int *pos, int **recpos_g,int *ntr_loc, int ntr_glob);
void read_su_pos(char srcname[256], float dh, float *refrec, int nt_in,
                 float dt_in, int ***srcpos, int *nsrc_glob);
int **evaluate_wfd_coords(float dh, float *refrec, int fdorder, int free_surf,
                          int fw, int nxg, int nyg, float memory_depth,
                          int source_snap_type, int *nsrc_glob);
void extract_snap_from_su(float **signal_sp, float *snap_sp, int nsrc_loc, int it);

void read_wfd_file(FILE *fp1, float *wfd, int num, int hin);
void write_wfd_file(FILE *fp1, float *wfd, int num, int hin);
void read_su_info(char filename[256], int *ntr, int *ns, float *dt);
float **read_su_header(char filename[256], int *ntr, int *ns, float *dt);
float **read_su_header_all(char filename[256], int *ntr, int *ns, float *dt, float *Sx, float *Sy);
float **read_su(char filename[256], int *ntr, int *ns);
void write_su(char filename[256], float Sx, float Sy, float **coord, float **seis, int ntr, int ns, float dt);

//void update_shot_num(char filein[256],int ishot,char fileout[256]);
void write_snap_boundary(FILE *fp1, float **sp, int **pos_loc, int num, int hin);

float **sources(char src_file[256], int *nsrc);
float ** wavelet(float ** srcpos_loc, int nsrc, int src_shape, int nt, float dt);

// sac file
float *read_sac(char filename[256], int *ns, float *dt);



