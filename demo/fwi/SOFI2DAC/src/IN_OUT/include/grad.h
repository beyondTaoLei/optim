/* files to include */
#include "comm.h"

void calc_correlation_time(int n1,int n2,int inv_vp,int inv_vs,int inv_rho,
        float *fvx, float *fvy, float *fsxx, float *fsyy, float *fsxy,
        float *bvx, float *bvy, float *bsxx, float *bsyy, float *bsxy, 
        float *cvx, float *cvy, float *csxx, float *csyy, float *csxy);
void gradient(int n1, int n2, float dt, int inv_vp, int inv_vs, 
              int inv_rho, float *Mrho, float *Mvs, float *Mvp,
              float *cvx,float *cvy,float *csxx,float *csyy,float *csxy,
              float *grad_lam, float *grad_mu, float *grad_rho_s,
              float *grad_vp, float *grad_vs, float *grad_rho);
// void eprecond1(int n1, int n2, float p1, int eprecond, float *Ws, float *Wr, float *We);
void save_grad(int n1, int n2, int inv_vp, int inv_vs, int inv_rho, 
               float *grad_vp, float *grad_vs, float *grad_rho, int ishot);

void calc_correlation_ac_time(int n1,int n2,int inv_vp,int inv_rho,
        float *fvx, float *fvy, float *fsp, float *bvx, float *bvy, 
        float *bsp, float *cvx, float *cvy, float *csp);
void gradient_ac(int n1, int n2, float dt, int inv_vp, int inv_rho, 
                 float *Mrho, float *Mvp, float *cvx,float *cvy,float *csp,
                 float *grad_lam, float *grad_rho_s,float *grad_vp, float *grad_rho);
void save_grad_ac(int n1, int n2, int inv_vp, int inv_rho, 
               float *grad_vp, float *grad_rho, int ishot);

