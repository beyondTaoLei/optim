/*------------------------------------------------------------------------
 *  fd.h - include file for viscoelastic FD programs          
 *  See COPYING file for copying and redistribution conditions.
 *  ---------------------------------------------------------------------*/

/* files to include */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#define iround(x) ((int)(floor)(x+0.5))
#define min(x,y) ((x<y)?x:y)    
#define max(x,y) ((x<y)?y:x)
#define fsign(x) ((x<0.0)?(-1):1)    

#define PI (3.141592653589793)
#define NPAR 120
#define STRING_SIZE 256
#define STRING_SIZE2 256
#define REQUEST_COUNT 4

#include "source.h"
#include "grad.h"
#include "ssg.h"
#include "struct_basic.h"
typedef struct Acquisition{
    //Acquisition
    int nsrc_loc;
    float Sx, Sy;
    int **srcpos_loc;
    float **signals_loc;
}st_acquisition;

typedef struct Wavefild_freq{
    int num_freqs;
    float *freqs;
    float ****pvx_f, ****pvy_f, ****psp_f;
}st_wavefildfreq;

/* declaration of functions */
void av_rho(float **rho, float **rip, float **rjp);
void av_tau(float **taus, float **tausipjp);
void	catseis(float **data, float **fulldata, int *recswitch, int ntr_glob, int ns);
void checkfd(FILE *fp, float ** prho, float ** ppi, float ** ptaup, float *peta, float *hc);
void comm_ini(float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
float ** buffertop_to_bot, float ** bufferbot_to_top, 
MPI_Request *req_send, MPI_Request *req_rec);

void cpml_update_p_x ( int i, int j, float  * vxx, float ** psi_vxx,
                       float * K_x, float * a_x, float * b_x);
void cpml_update_p_y ( int i, int j, float * vyy, float ** psi_vyy,
                       float * K_y, float * a_y, float * b_y);
void cpml_update_v_x ( int i, int j,float  * sp_x, float ** psi_sxx_x,
                       float * K_x_half, float * a_x_half, float * b_x_half);
void cpml_update_v_y ( int i, int j,float * sp_y, float ** psi_syy_y,
                       float * K_y_half, float * a_y_half, float * b_y_half);
void exchange_p(int nd, float ** sp, float ** bufferlef_to_rig, float ** bufferrig_to_lef,
		float ** buffertop_to_bot, float ** bufferbot_to_top,
		MPI_Request *req_send, MPI_Request *req_rec);
void exchange_par(void);
void exchange_v(int nd, float ** vx, float ** vy,
		float ** bufferlef_to_rig, float ** bufferrig_to_lef,
		float ** buffertop_to_bot, float ** bufferbot_to_top,
		MPI_Request * req_send, MPI_Request * req_rec);
float **holbergcoeff_all();
float *holbergcoeff();
void info(FILE *fp);
void initproc(void);
int read_objects_from_intputfile(FILE *fp, char *input_file,char ** varname_list,char ** value_list);
void print_objectlist_screen(FILE *fp, int number_readobject,char ** varname_list,char ** value_list);
int count_occure_charinstring(char stringline[STRING_SIZE2], char teststring[]);
void copy_str2str_uptochar(char string_in[STRING_SIZE2], char string_out[STRING_SIZE2], char teststring[]);
int get_int_from_objectlist(char string_in[STRING_SIZE2], int number_readobject, int * int_buffer,
                            char ** varname_list,char ** value_list);
int get_float_from_objectlist(char string_in[STRING_SIZE2], int number_readobject, float * double_buffer,
                              char ** varname_list,char ** value_list);
int get_string_from_objectlist(char string_in[STRING_SIZE2], int number_readobject, char string_buffer[STRING_SIZE2],
                               char ** varname_list,char ** value_list);
int is_string_blankspace(char string_in[STRING_SIZE2]);
void remove_blankspaces_around_string(char string_in[STRING_SIZE2] );
void add_object_tolist(char string_name[STRING_SIZE2],char string_value[STRING_SIZE2], int * number_readobject
                       , char ** varname_list,char ** value_list);

void prepare_boundary_saving(int ishot, int nshot1, st_bdr *bdr);
void write_snap_boundary_file(int ishot, int nt, float **psp, float **pvx, float **pvy, st_bdr *bdr);
void finish_boundary_file(st_bdr *bdr);
void free_boundary_struct_out_shot(st_bdr *bdr);

//void set_wfd_pointer_ac(float *pvx, float *pvy, float *psp, st_wfd *wfd);
void allocate_mem_wfd(int flag, int ng_loc, st_wfd *wfd);
void allocate_mem_wfd_e(int flag, int ng_loc, st_wfd *wfd);
void allocate_mem_grad(int flag, int ng_loc, st_grad *grad);
void trans_para(float **rho, float **pi, st_grad *grad);
void project_wfd(float **vx, float **vy, float **sp, int **win_loc, st_wfd *wfd);
void calc_correlation(float **fvx, float **fvy, float **fsp,
                      float **bvx, float **bvy, float **bsp,
                      float **cvxxyy, float **csp);
void calc_gradient(float **cvxxyy, float **csp, 
                   float **prho, float **ppi);
void savematIDXY(float **mod,char modfile[STRING_SIZE]);
void savematIDXY2(float **fullmod,char modfile[STRING_SIZE]);
void eprecond1(float ** We, float ** Ws, float ** Wr);
void update_precond(float **eng, float prctile);
void apply_precond(float **eng, float **grad);
float global_percentiles(float p1, float ** fullmod);
void quicksort2(float *arr, int dummy, int elements);
int find_min_positive_num2(int n1,int n2,float *a);
void matrix_cat(int nxg, int nyg, int nx, int ny, int pos1, int pos2, 
                float **mod, float **fullmod, MPI_Comm mpi_comm_shot);
void matrix_scatter(int nxg, int nyg, int nx, int ny, int pos1, int pos2,
                    float **mod, float **fullmod);
void matrix_norm(int nx1, int nx2,int ny1, int ny2,int idx, int idy, float ** matrix);
void allocate_mem_seis(int flag, float ***sectionvx, float ***sectionvy, float ***sectionp, 
                       float ***sectiondiv, int ntr, int ns);
void allocate_mem_seis2(int flag, float ***sectionvx, float ***sectionvy, float ***sectionp, 
                       int ntr, int ns);
void save_seis(FILE *fp, float **sectionvx, float **sectionvy, 
               float **sectionp, float **sectiondiv,
               float ** seismo_fulldata, int * recswitch, int **recpos, 
               int **recpos_loc, int ntr_glob, int ns, int ishot,
               float Sx, float Sy);
void save_seis2(FILE *fp, float **sectionvx, float **sectionvy, float **sectionp, 
                float ** seismo_fulldata, int * recswitch, int **recpos, 
                int **recpos_loc, int ntr_glob, int ns, int ishot, float Sx, float Sy);

int **create_window(int nx1, int nx2, int ny1, int ny2, int idx, int idy, int *icount);
int **create_window_inv(int *icount);
void catvector_mpi(float *v, int **win, float *g, int n1, int n2, int ng_glob);

void psource(float **sp, int **srcpos_loc, float **signals, int nt, int nsrc, int sw);
void source_boundary(float **sp, int **srcpos_loc, float **signals, int nt, int nsrc);
void source_last_snap(int ishot, float **psp, int flag);

void init_struct_src_elmts(int flag, st_source *geo_src);
void prepare_src_position(int ishot, int nshot1, st_source *geo_src);
void prepare_sig_file(int ishot, st_source *geo_src);
void extract_snap_time(int nt, int nt_snap, st_source *geo_src);
void insert_particle_velocity(float **pvx, float **pvy, float  **prip, float **prjp, st_source *geo_src);
void insert_velocity_last_snap(int ishot, float **pvx, float **pvy, st_source *geo_src);
void insert_stress(float **psp, st_source *geo_src);
void insert_stress_last_snap(int ishot, float **psp, st_source *geo_src);
void finish_src_file(int ns, st_source *geo_src);
void free_src_struct_out_shot(int nx, int ny, st_source *geo_src);


void matcopy_acoustic(float ** rho, float ** pi);
void matcopy_viscac(float ** rho, float ** pi, float ** taup);

void write_matrix_disk(float ** local_matrix,char path_name[STRING_SIZE]);
float global_maximum(float ** gradiant_1);
float average_matrix(float ** matrix);
float matrix_product(float ** matrix1, float **matrix2);
float ** get_global_from_local_matrix(float ** local_matrix);
void get_local_from_global_matrix(float ** global_matrix,float ** local_matrix);


void merge(int nsnap, int type);
void mergemod(char modfile[STRING_SIZE], int format);

void note(FILE *fp);
void  outseis_glob(FILE *fp, FILE *fpdata, float **section,
		int **recpos, int **recpos_loc, int ntr, float Sx, 
                float Sy, int nsrc, int ns, int seis_form, int comp);
void PML_pro(float * d_x, float * K_x, float * alpha_prime_x, float * a_x, float * b_x, 
            float * d_x_half, float * K_x_half, float * alpha_prime_x_half, float * a_x_half, float * b_x_half,
            float * d_y, float * K_y, float * alpha_prime_y, float * a_y, float * b_y, 
            float * d_y_half, float * K_y_half, float * alpha_prime_y_half, float * a_y_half, float * b_y_half);
void prepare_update_p(float *etajm, float *peta, float **ppi, float **prho, float **ptaup, float **g, float *bjm, float *cjm, float ***e);

void read_par_json(FILE *fp, char *fileinp);
float readdsk(FILE *fp_in, int format);
void readmod_acoustic(float  **  rho, float **  pi);
void readmod_viscac(float  **  rho, float **  pi, float ** taup, float * eta);
int **receiver(FILE *fp, int *ntr);
int **receiver_boundary(FILE *fp, int *ntr);

void saveseis_glob(FILE *fp, float **sectiondata, int  **recpos, int  **recpos_loc, int ntr, float Sx, 
                float Sy, int ishot,int ns, int sectiondatatype);
void saveseis_glob2(FILE *fp, float **sectiondata, int  **recpos, int  **recpos_loc, int ntr, float Sx, 
                float Sy, int ishot,int ns, int sectiondatatype);
void seismo_ssg(int lsamp, int ntr, int **recpos, float **sectionvx,
                float **sectionvy, float **sectionp, float **sectiondiv,
                float **vx, float **vy, float **sp, float **pi, float *hc);
void seismo_ssg2(int lsamp, int ntr, int **recpos, float **sectionvx,
                 float **sectionvy, float **sectionp,float **vx, 
                 float **vy, float **sp);

float *read_freqs(char filename[256], int *nrow);
void allocate_mem_wfdf(int flag, st_wavefildfreq *wfd_f);


void snap(FILE *fp,int nt, int nsnap, float **vx, float **vy, float **sxx,
          float **syy, float **sp, float **u, float **pi, float *hc, int ishot, st_wavefildfreq *wfd_f);
void snap_last(FILE *fp, float **vx, float **vy,float **sp, int ishot);
void calc_grid_for_inv(void);
void eprecond(float ** W, float ** vx, float ** vy);
int *generate_frequencies(int *num_freqs);
void dft_forward(int nt, int num_freqs, int *freqs, float **vx,float ****Fvx);
void dft_forward2(int nt, int num_freqs, float *freqs, float **vx,float ****Fvx);
void save_wavefields_freq(FILE *fp, st_wavefildfreq *wfd_f, int ishot);
// void zero_wavefields_freq(int num_freqs, float ****vx_f, float ****vy_f, 
//                           float ****sp_f);
// void calc_wavefields_freq(float ****vx_f, float ****vy_f, float ****sp_f, 
//                           float **vx, float **vy, float **sp, int *freqs,
//                           int nt, int num_freqs);
void extract_current_sources(char srcname[256], st_acquisition *acq);
// void split_su_pos_loc(char srcname[256], int myid, int myid1, 
//                       MPI_Comm mpi_comm_no, float dh, float *refrec, 
//                       int iendx, int iendy, int *POS1, int nt_in, 
//                       float dt_in, int ***srcpos_loc, int *nsrc_loc);
// void split_su_sig_loc(char srcname[256], int myid, int myid1, 
//                       MPI_Comm mpi_comm_no, int **srcpos_loc, 
//                       int ntr_loc, float ***signals);
// void split_wfd_pos_loc(int myid, int myid1, MPI_Comm mpi_comm_no, 
//                         int iendx, int iendy, int *POS1, float dh, 
//                         float *refrec, int fdorder, int free_surf, 
//                         int fw, int nxg, int nyg, float memory_depth, 
//                         int source_snap_type, char snap_coords_file[256],
//                         int *nsrc_glob, int ***srcpos_loc, int *nsrc_loc);
// void split_wfd_sig_loc(FILE *fp_sp, int myid, int myid1, 
//                        MPI_Comm mpi_comm_no, float *snap_all, 
//                        float *snap_loc, int ** recpos_loc, 
//                        int nsrc_glob, int nsrc_loc, int it);
int **splitrec(int **recpos,int *ntr_loc, int ntr, int *recswitch);
float **splitsrc(float **srcpos,int *nsrc_loc, int nsrc);


/* SSG_fd.c */
// void SSG_fd_v(float *vxx, float *vxy, float *vyy, float *vyx, 
//            float **vx, float **vy, int i, int j, int fdoh, float *hc);
// void SSG_fd_s(float *sxx_x, float *sxy_x, float *sxy_y, float *syy_y, 
//            float ** sxx, float ** syy, float ** sxy, int i, int j, 
//            int fdoh, float *hc);
// void SSG_fd_s_ac(float *sp_x, float *sp_y, float **sp, int i, int j, 
//                int fdoh, float *hc);
// void SSG_fd_v_ac(float *vxx, float *vyy, float **vx, float **vy, 
//                  int i, int j, int fdoh, float *hc);
// float calc_vxx(float **vx, int i, int j, int fdoh, float *hc);
// float calc_vyy(float **vy, int i, int j, int fdoh, float *hc);
// float calc_vxy(float **vx, int i, int j, int fdoh, float *hc);
// float calc_vyx(float **vy, int i, int j, int fdoh, float *hc);
// float calc_sxx_x(float **sxx, int i, int j, int fdoh, float *hc);
// float calc_sxy_x(float **sxy, int i, int j, int fdoh, float *hc);
// float calc_sxy_y(float **sxy, int i, int j, int fdoh, float *hc);
// float calc_syy_y(float **syy, int i, int j, int fdoh, float *hc);
// float calc_sp_x(float **sp, int i, int j, int fdoh, float *hc);
// float calc_sp_y(float **sp, int i, int j, int fdoh, float *hc);


void subgrid_bounds ( int nx1, int nx2, int ny1, int ny2, int * gx, int * gy );
void surface_acoustic_PML(int ndepth, float ** sp);
double time_count(double time1, double *time2);

// void update_p_interior(int *gx, int *gy, float **rho, float **pi, float *hc, 
//                        float **vx, float **vy, float **sp, float **u);
void update_p_PML(int *gx, int *gy, float **rho, float **pi, float *hc,
                  float **vx, float **vy, float **sp, float **u,
                  float **psi_vxx, float **psi_vyy,
                  float *K_x, float *a_x, float *b_x, 
                  float *K_y, float *a_y, float *b_y, 
                  int nx2, int ny2, float flag);
void update_p_visc_PML(int nx1, int nx2, int ny1, int ny2, float ** vx, float ** vy, float ** sp, float ** pi, float **rho, float *hc, int infoout,
		       float ***p, float **g, float *bjm, float *cjm, float ***e,  
		       float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
		       float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
		       float ** psi_vxx, float ** psi_vyy, float ** psi_vxy, float ** psi_vyx);
// void update_v_interior(int *gx, int *gy, float  **rip, float **rjp, float *hc,
//                        float **vx, float **vy, float **uvx, float **uvy, float **sp);
void update_v_PML(int *gx, int *gy, float  **rip, float **rjp, float *hc,
                  float **vx, float **vy, float **uvx, float **uvy, float **sp,
                  float **psi_sxx_x, float **psi_syy_y,
                  float *K_x_half, float *a_x_half, float *b_x_half, 
                  float *K_y_half, float *a_y_half, float *b_y_half,
                  int nx2, int ny2, float flag);

void update_p_interior(int *gx, int *gy, float **rho, float **pi, float *hc, 
                       float **vx, float **vy, float **sp, float **u, float flag);
// void update_p_interior_back(int *gx, int *gy, float **rho, float **pi, float *hc, 
//                        float **vx, float **vy, float **sp, float **u);
void update_v_interior(int *gx, int *gy, float  **rip, float **rjp, float *hc,
                       float **vx, float **vy, float **uvx, float **uvy, float **sp, float flag);
// void update_v_interior_back(int *gx, int *gy, float  **rip, float **rjp, float *hc,
//                        float **vx, float **vy, float **uvx, float **uvy, float **sp);
// void wavefield_update_p_b(int i, int j, float **rho, float **pi,
//                         float vxx, float vyy,  float **sp, float **u);
// void wavefield_update_v_b(int i, int j, float **rip, float **rjp, 
//                         float sp_x, float sp_y, float **vx, float **vy, 
//                         float **uvx, float **uvy);


void declare_error(char err_text[]);
void warning(char warn_text[]);
double maximum(float **a, int nx, int ny);
float minimum_m(float **mat, int nx, int ny);
float maximum_m(float **mat, int nx, int ny);
float *vector(int ni, int nj);
int *ivector(int ni, int nj);
unsigned short int *usvector(int ni, int nj);
unsigned char *cvector(int ni, int nj);
unsigned long *lvector(int ni, int nj);
double *dvector(int ni, int nj);
float **fmatrix(int mrl, int mrh, int mcl, int mch);
float **matrix(int mrl, int mrh, int mcl, int mch);
double **dmatrix(int mrl, int mrh, int mcl, int mch);
int **imatrix(int mrl, int mrh, int mcl, int mch);
unsigned short int **usmatrix(int mrl, int mrh, int mcl, int mch);
float ***f3tensor(int mrl, int mrh, int mcl, int mch,int mdl, int mdh);
int ***i3tensor(int mrl, int mrh, int mcl, int mch,int mdl, int mdh);
float ****f4tensor(int nrl, int nrh, int ncl, int nch,int ndl, int ndh, int nvl, int nvh);
void free_vector(float *a, int ni, int nj);
void free_ivector(int *a, int ni, int nj);
void free_cvector(char *a, int ni, int nj);
void free_dvector(double *a, int ni, int nj);
void free_matrix(float **ma, int mrl, int mrh, int mcl, int mch);
void free_imatrix(int **ma, int mrl, int mrh, int mcl, int mch);
void free_usmatrix(unsigned short int **ma, int mrl, int mrh, int mcl, int mch);
void free_f3tensor(float ***te, int mrl, int mrh, int mcl, int mch, int mdl, int mdh);
void free_i3tensor(int ***te, int mrl, int mrh, int mcl, int mch, int mdl, int mdh);
void free_f4tensor(float ****t, int nrl, int nrh, int ncl, int nch, int ndl, int ndh, int nvl,int nvh);
void zero(float *A, int u_max);
void normalize_data(float **data, int ntr, int ns);
// void quicksort(float *arr, int dummy, int elements);


// void wavefield_update_p(int i, int j, float **rho, float **pi,
//                         float vxx, float vyy,  float **sp, float **u);
// void wavefield_update_p_new(int i, int j, float **rho, float **pi,
//                         float vxx, float vyy,  float **sp, float **u,
//                         float flag);
// void wavefield_update_v(int i, int j, float **rip, float **rjp, 
//                         float sp_x, float sp_y, float **vx, float **vy, 
//                         float **uvx, float **uvy);
// void wavefield_update_v_new(int i, int j, float **rip, float **rjp, 
//                         float sp_x, float sp_y, float **vx, float **vy, 
//                         float **uvx, float **uvy, float flag);
void write_par(FILE *fp);
void writedsk(FILE *fp_out, float amp, int format);
void writemod(char modfile[STRING_SIZE], float ** array, int format);
void zero_acoustic(int ny1, int ny2, int nx1, int nx2, float **vx, float **vy, float **sp);
void zero_fdveps_ac(int ny1, int ny2, int nx1, int nx2, float ** vx, float ** vy, float ** sp, float ** vxp1, float ** vyp1,
                 float ** psi_sxx_x, float ** psi_sxy_x, float ** psi_vxx, float ** psi_vyx, float ** psi_syy_y, float ** psi_sxy_y, float ** psi_vyy, float ** psi_vxy,
                 float ** psi_vxxs);
void zero_fdveps_viscac(int ny1, int ny2, int nx1, int nx2, float ** vx, float ** vy, float ** sp, float ** vxp1, float ** vyp1,
			float ** psi_sxx_x, float ** psi_sxy_x, float ** psi_vxx, float ** psi_vyx, float ** psi_syy_y,
			float ** psi_sxy_y, float ** psi_vyy, float ** psi_vxy, float ** psi_vxxs, float ***pp);
