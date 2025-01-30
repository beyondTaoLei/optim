/* files to include */
#include "comm.h"

/* SSG_fd.c */
void SSG_fd_v(float *vxx, float *vxy, float *vyy, float *vyx, 
           float **vx, float **vy, int i, int j, int fdoh, float *hc);
void SSG_fd_s(float *sxx_x, float *sxy_x, float *sxy_y, float *syy_y, 
           float ** sxx, float ** syy, float ** sxy, int i, int j, 
           int fdoh, float *hc);
void SSG_fd_s_ac(float *sp_x, float *sp_y, float **sp, int i, int j, 
               int fdoh, float *hc);
void SSG_fd_v_ac(float *vxx, float *vyy, float **vx, float **vy, 
                 int i, int j, int fdoh, float *hc);
float calc_vxx(float **vx, int i, int j, int fdoh, float *hc);
float calc_vyy(float **vy, int i, int j, int fdoh, float *hc);
float calc_vxy(float **vx, int i, int j, int fdoh, float *hc);
float calc_vyx(float **vy, int i, int j, int fdoh, float *hc);
float calc_sxx_x(float **sxx, int i, int j, int fdoh, float *hc);
float calc_sxy_x(float **sxy, int i, int j, int fdoh, float *hc);
float calc_sxy_y(float **sxy, int i, int j, int fdoh, float *hc);
float calc_syy_y(float **syy, int i, int j, int fdoh, float *hc);
float calc_sp_x(float **sp, int i, int j, int fdoh, float *hc);
float calc_sp_y(float **sp, int i, int j, int fdoh, float *hc);

void wavefield_update_p(int i, int j, float **rho, float **pi,
                        float vxx, float vyy,  float **sp, float **u,
                        float dt, float dh, float flag);
void wavefield_update_v(int i, int j, float **rip, float **rjp, 
                        float sp_x, float sp_y, float **vx, float **vy, 
                        float **uvx, float **uvy, float dt, float dh, float flag);

void wavefield_update_s_el ( int i, int j,float   vxx, float  vyx,float vxy,float  vyy, float **sxy,
                   float **sxx, float ** syy, float ** pi, float ** u, float ** uipjp );
void wavefield_update_s_el_new(int i, int j, float **pi, float **u, float **uipjp,
                               float vxx, float vyx, float vxy, float vyy, 
                               float **sxy, float **sxx, float **syy, 
                               float **uxy, float **ux, float ** uy,
                               float dt, float dh, float flag);
void wavefield_update_v_el ( int i, int j,float   sxx_x, float  sxy_x,float sxy_y,float  syy_y, float **vx,
                          float **vy, float ** rip, float ** rjp );
void wavefield_update_v_el_new(int i, int j, float **rip, float **rjp, 
                            float sxx_x, float sxy_x, float sxy_y, float syy_y, 
                            float **vx, float **vy, float **uvx, float **uvy, 
                            float dt, float dh, float flag);

// void wavefield_update_v_el_new2(int i, int j, float rip, float rjp, 
//                             float sxx_x, float sxy_x, float sxy_y, float syy_y, 
//                             float *vx, float *vy, float *uvx, float *uvy, 
//                             float dt, float dh, float flag);
// void wavefield_update_s_el_new2(int i, int j, float pi, float u, float uipjp,
//                                float vxx, float vyx, float vxy, float vyy, 
//                                float *sxy, float *sxx, float *syy, 
//                                float *uxy, float *ux, float * uy,
//                                float dt, float dh, float flag);

