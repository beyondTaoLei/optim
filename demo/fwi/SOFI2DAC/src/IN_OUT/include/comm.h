/* files to include */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>

#define iround(x) ((int)(floor)(x+0.5))
#define min(x,y) ((x<y)?x:y)
#define max(x,y) ((x<y)?y:x)
#define fsign(x) ((x<0.0)?(-1):1)

#define PI (3.141592653589793)
#define STRING_SIZE 256
#define STRING_SIZE2 256

// util.c
/* very great utilities  from sofi2d software */
//void declare_error(char err_text[]);
//void warning(char warn_text[]);
void throw_error(char err_text[]);
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
char **chmatrix(int mrl, int mrh, int mcl, int mch);
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
void free_chmatrix(char **ma, int mrl, int mrh, int mcl, int mch);
void zero(float *A, int u_max);
void normalize_data(float **data, int ntr, int ns);


/* vector operations */
//void quicksort(float *arr, int start_index, int elements);
void quicksort1(float *arr, int start_index, int elements);
//int find_min_positive_num(int n1,int n2,float *a);
int find_min_positive_num1(int n1,int n2,float *a);
float global_percentile(float *fullvec, int n1, int n2, float p1);
int get_file_size(char filename[256]);
void V_write(char modfile[256],float *x_vector,int nn);
void V_write_txt(char modfile[256],float *x_vector,int nn);
void V_read(char modfile[256],float *x_vector,int nn);
float V_max_abs(int n1, int n2, float * x);
int V_max_index(int n1, int n2, float * x);
int V_max_abs_index(int n1, int n2, float * x);
void V_normalize(int n1, int n2, float * x);
float V_normL2(int n1, int n2, float * x);
float V_dotproduct(int n1, int n2, float *x, float *y);
/* operate vectors element-wise.  */
void V_copy(int n1, int n2, float *x, float *r);
void V_pick(int n1, int n2, float *x, float *r);
void V_abs(int n1, int n2, float *x, float *r);
void V_add(int n1, int n2, float *x, float *y,float *r);
void V_subtract(int n1, int n2, float *x, float *y,float *r);
void V_multiply(int n1, int n2, float *x, float *y,float *r);
void V_divide(int n1, int n2, float *x, float *y,float *r);
void V_square_sum(int n1, int n2, float *x, float *y,float *r);
void V_square_sum_a(int n1, int n2, float *x, float *y,float *r);
void V_AXPY(int n1, int n2, float a, float * x, float * y,float *r);
void V_AXPb(int n1, int n2, float a, float * x,  float b, float *r);
void V_scaling(int n1, int n2, float *x, float a);
void V_sqrt(int n1, int n2, float *x, float *r);
void V_multiply_sqrt(int n1, int n2, float *x, float *y, float *r);
void V_multiply_a(int n1, int n2,float *a, float *b,float *cor);
void V_multiply_a2(int n1, int n2, float *a, float *b, float *c, 
                   float *d, float *cor1, float *cor2);
void V_grad(int n1, int n2, float *x, float *r);
void V_integ(int n1, int n2, float *x, float *r);
void V_shift(int n1, int n2, int lag, float *x, float *r);
void V_reversal(int n1, int n2, float *x, float *r);


/* matrix operations */
float M_max_abs(int nx1, int nx2,int ny1, int ny2, float **x);
void M_normalize(int nx1, int nx2,int ny1, int ny2, float **x);
float M_normL2(int nx1, int nx2,int ny1, int ny2, float **x);
/* operate matrixes element-wise.  */
void M_init(int nx1, int nx2,int ny1, int ny2, float **x);
void M_copy(int nx1, int nx2,int ny1, int ny2, float **x, float **r);
void M_add(int nx1, int nx2,int ny1, int ny2, float **x, float **y,float **r);
void M_subtract(int nx1, int nx2,int ny1, int ny2, float **x, float **y,float **r);
void M_multiply(int nx1, int nx2,int ny1, int ny2, float **x, float **y,float **r);
void M_divide(int nx1, int nx2,int ny1, int ny2, float **x, float **y,float **r);
void M_square_sum(int nx1, int nx2,int ny1, int ny2, float **x, float **y,float **r);
void M_square_sum_a(int nx1, int nx2,int ny1, int ny2, float **x, float **y,float **r);
void M_AXPY(int nx1, int nx2,int ny1, int ny2, float a, float **x, float **y,float **r);
void M_AXPb(int nx1, int nx2,int ny1, int ny2, float a, float **x,  float b, float **r);
void M_scaling(int nx1, int nx2,int ny1, int ny2, float **x, float a);
void M_sqrt(int nx1, int nx2,int ny1, int ny2, float **x, float **r);
void M_multiply_sqrt(int nx1, int nx2,int ny1, int ny2, float **x, float **y, float **r);
void M_multiply_a(int nx1, int nx2,int ny1, int ny2,float **a, float **b,float **cor);
void M_multiply_a2(int nx1, int nx2,int ny1, int ny2, float **a, float **b, float **c, 
                   float **d, float **cor1, float **cor2);


// project.c
// int count_elements(int nx1, int nx2, int ny1, int ny2, int idx, int idy);
// void create_win(int nx1, int nx2, int ny1, int ny2, int idx, int idy, 
//                 int *winx, int *winy);
void project2vec(float **m, float *v, int **win, int n1, int n2);
void project2vec2(float **p, float *v, int nx1, int nx2, int ny1,
                      int ny2, int idx, int idy);
void project2arr(float **p, float *v, int nx1, int nx2, int ny1,
                      int ny2, int idx, int idy);

// xcorr.c
void xcorr(float * temp_TS, float * temp_TS1, float * temp_xcorr, int maxlag, int ns);

//envelope.c
//void envelope(float *signals, float *env, int ns);
void hilbert(float *signals, float *hilb, int ns);
void envelope_hilbert(float *signals, float *hilb, float *env, int ns);
void phase(float *signals, float *hilb, float *phase, int ns);
void phase2(float *signals, float *hilb, float *phase, int ns);
// void hilbert2(float *signals, float *hilb, int ns);
void phase_diff(float *signals_o, float *hilb_o, float *signals_s, 
                float *hilb_s, float *phase_diff, int ns);

//windows.c
// void time_window(float *time_win, int nlen);
float *generate_freqs(int ns, float dt, int fac, int *nfreq);
void gauss_window(int nfreq, float *freqs, float *win, float alfa, float fc);
void cosine_window(int nfreq, float *freqs, float *win, float *fc);
void cosine_window2(int nfreq, float *freqs, float *win, float fc);
void Taper(float *data, int type, int nlen, float width);

//interpolation.c
float linint(float x0, float y0, float x1, float y1, float x);
float linint_set(float *x, float *y, float xx, int n);

//freqdomain_filt_vector.c
void freqdomain_filt_vector(float *data, float *f_win, int nfreq, int ns);

//read_text.c
float **read_text(int flag, int ncol, int *nrow, char filename[256]);

//fft_use.c
void fft_amp_phase(float *signals, float *phase, int ns);
void fft_real_imag(float *signals, float *imag, int ns, int flag);

//vector_operations2.c
void V_der(int ns, float dt, float *x, float *r);
void V_int(int ns, float dt, float *x, float *r);

//json_parser.c
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

//string_opts.c
char *strrpc(char *str,char *oldstr,char *newstr);
void strrpc2(char str[256],char *oldstr,char *newstr, char str_out[256]);
int StringToInteger(char *p);
char **read_string(int *nrow, char filename[256]);
char **read_string1(int *nrow, char filename[256]);

//update_shot_num.c
void update_shot_num(char filein[256],int ishot,char fileout[256]);
void update_shot_num2(char filein[256],int ishot,char fileout[256]);
int extract_shot_num(char filein[256]);

//intersection.c
void SelectSort(int *nums,int numsSize);
int* intersection(int* nums1, int nums1Size, int* nums2, int nums2Size, int* returnSize);
