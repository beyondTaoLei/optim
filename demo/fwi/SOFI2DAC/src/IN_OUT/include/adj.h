/* files to include */
#include "comm.h"
void EnvelopeImage(float *data, float *TPoint, int NumCtrT, int ns, 
                    float dt, float StaDist, float **env);
float misfit_L2(float *obs, float *syn, float *res, float tstart, 
                float tend, int ns, float dt);
float misfit_WF(float *obs, float *syn, float *res, float tstart, 
                        float tend, int ns, float dt);
float misfit_IP(float *obs, float *syn, float *res, int ns, float dt, float eps);
float misfit_ENV(float *obs, float *syn, float *res, float tstart, 
                 float tend, int ns, float dt, float eps);
//float misfit_TT(float *obs, float *syn, float *res, float tstart, 
//                        float tend, int ns, float dt, float *delay);
// float misfit_TT_kernel(float *obs, float *syn, float *res, float tstart, 
//                         float tend, int ns, float dt);
float misfit_TT(float *obs, float *syn, float *res, float tstart, 
                float tend, int ns, float dt, int flag, float *delay);
float misfit_TT_multi(float *obs, float *syn, float *res, float tstart, 
                float tend, int ns, float dt, float **freqs, int nfreq, int order);
float misfit_TT_multi_debug(float *obs, float *syn, float *res, float tstart, 
                float tend, int ns, float dt, float **freqs, int nfreq, int order,float *lags);
// float misfit_TT_kernel(float *obs, float *syn, float *res, float tstart, 
//                         float tend, int ns, float dt, float *delay);
float misfit_AMP(float *obs, float *syn, float *res, float tstart, 
                        float tend, int ns, float dt);
float misfit_Rytov(float *obs, float *syn, float *res, float tstart, 
                float tend, int ns, float dt, int flag, float *delay);
float misfit_CC(float *obs, float *syn, float *res, float tstart, 
                float tend, int ns, float dt, float *delay);

// void misfit_PH_Freq(float *obs, float *syn, float *res, float tstart, float tend, 
//                 int ns, float dt);



