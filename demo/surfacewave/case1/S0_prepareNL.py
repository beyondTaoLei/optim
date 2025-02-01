#!/usr/bin/env python3
"""
prepares input json file for surface wave dispersion inversion problem
Usage:
    
"""

import os
import sys
import json
import numpy as np
import pandas as pd

# Input paras
fdir    =sys.argv[1] # optim
# make directory
os.makedirs(fdir, exist_ok=True)
os.makedirs('fig', exist_ok=True)
os.makedirs('log', exist_ok=True)
os.makedirs('job', exist_ok=True)

optim={
    'optimroot'             :'/home/tao/Downloads/optim', #[IN] software path
    'Forward part'          :'comment',
    'mpiexep'               :'mpirun',
    'algo'                  :'dunkin',                  #[IN] algorithm for disperison calculation
    'wavetype'              :'rayleigh',                #[IN] 'love' or 'rayleigh'
    'maxmode'               :0,                         #[IN] code, fund. 0, first high mode 1
    'ncores'                :8,                         #[IN] number of cores you can apply for modelling
    'fstn'                  :'list/stn.csv',            #[IN] the measure coordinates
    'fperiods'              :'list/periods.csv',        #[IN] measured periods in second, starting with low periods
    'fmodfd'                :'model/inv/inv.vs',        #[IN] the current model in work memory
    'nxyzfd'                :[1, 1300, 400],            #[IN] X-, Y- and Z-dir grid number of input model
    'oxyzfd'                :[0.0, 0.0, 0.0],           #[IN] X-, Y- and Z-dir original coordinate of input model in meter
    'dxyzfd'                :[500.0, 500.0, 500.0],     #[IN] X-, Y- and Z-dir spacing of input model in meter
    'Gradient part'         :'comment',
    'smooth_size'           :[0.0, 30000.0, 9000.0],    #[IN] smooth size for descent direction, in meter (such as, 6*d3g:6*d2g:6*d1g)
    'precond'               :1,                         #[IN] preconditional operator for gradient.
    'EPSILON'               :0.05,                      #[IN] preconditional gradient is grad/(diagHessian+diagHessian.mean(axis=0).max()*EPSILON).
    'Optimization'          :'comment',
    'NLCG'                  :'comment',
    'powell'                :1,                         #[IN] wheter to check Powell restart condition(1-yes, 0-no)
    'l-BFGS'                :'comment',
    'l'                     :5,                         #[IN] the maximum storing number of the pairs {sk,yk}(~3-7)
    'TN'                    :'comment',
    'conv_CG'               :0,                         #[IN] initial status
    'niter_max_CG'          :12,                        #[IN] maximum number of inner conjugate gradient loop
    'flag_singularity'      :1,                         #[IN] 1: singularity test; 0: no test
    'flag_negative'         :2,                         #[IN] 1: negative curvature test; 2: strong negative curvature test; 0: no test
    'flag_truncation'       :1,                         #[IN] 1: truncation test based on residuals; 2: test based on parabolic function, not work now; 0: no test
    'flag_min_res'          :1,                         #[IN] 1: record minimum residual when applying the inner CG agrithm; 0: no record
    'pro_sing'              :1.0E-10,                   #[IN] threshold \varsigma for singularity test
    'pro_curv'              :1.0E-10,                   #[IN] threshold \delta for negative curvature test (flag_negative = 1)
    'pro_trun'              :0.5,                       #[IN] threshold c_q for truncation test
    'eta'                   :0.0,                       #[IN] if ||residuals|| <= eta * ||gradient||, skip the inner CG loop
    'eta_update'            :0,                         #[IN] 1: update the eta during inversion cycle; 0: not update
    'Steplength'            :'comment',
    'try_old_sl'            :1,                         #[IN] whether to use steplength from last iter.(1-yes, 0-no)
    'max_m0'                :4700.0,                    #[IN] the maximum of initial model
    'eps_scale'             :0.03,                      #[IN] max_m0*eps_scale as one unit to test the steplength
    'alpha_lb'              :0.1,                       #[IN] lower limit of steplength
    'alpha_ub'              :10,                        #[IN] upper limit of steplength
    'mod_limits'            :1,                         #[IN] whether to check the bound of model value, 1: yes, 0: no
    'mod_lb'                :2500,                      #[IN] lower limit of model value if mod_limits=1
    'mod_ub'                :5000,                      #[IN] upper limit of model value if mod_limits=1
    'Inversion loop'        :'comment',
    'abort'                 :0,                         #[IN] intial stage phase
    'iter0'                 :1,                         #[IN] the initial iteration at current inversion stage/phase
    'niter_min'             :30,                        #[IN] minimum iterations
    'niter_max'             :40,                        #[IN] maximum iterations, for future design
    'pro'                   :0.01,                      #[IN] (Ei2 - Ei) / Ei2
    'pro2'                  :0.01,                      #[IN] k2 / k1 
    'pro3'                  :0.01,                      #[IN] average2 / average1 
    'time_sleep'            :10,                        #[IN] time interval for checking whether current step has completed, in second
}
optim['demoroot'] = os.path.dirname(os.path.dirname(__file__)) #.../surfacewave
optim['prjroot'] = os.getcwd() # project path
optim['fdir'] = fdir
optim['odir'] = 'obs'                               #[IN] work directory for observed data
optim['wdir'] = 'syn'                               #[IN] work directory for synthetic data
optim['fmod'] = optim['fmodfd']                     # grid of FD = grid of inversion
optim['fcost'] = 'list/fcost.dat'                   #[OUT] misfit function
optim['fgrad'] = 'list/grad.bin'                    #[OUT] the gradient at point fmod
optim['fhess'] = 'list/hess.bin'                    #[OUT] the Hessian vector product at point [fmod, fvctr]
optim['fvctr'] = 'list/vctr.bin'                    #[OUT] the vector pk in LCG algorithm
# inverse grids
optim['n3'] = optim['nxyzfd'][0]
optim['n2'] = optim['nxyzfd'][1]
optim['n1'] = optim['nxyzfd'][2]
optim['d3'] = optim['dxyzfd'][0]
optim['d2'] = optim['dxyzfd'][1]
optim['d1'] = optim['dxyzfd'][2]
optim['n'] = optim['n3']*optim['n2']*optim['n1'] #only for TRN now

# write work infomation
foptim = os.path.join(fdir,'optim.json')
fout = open(foptim,'w')
fout.write(json.dumps(optim,indent=4))
fout.close()
print('*'*60+'\n*Output file: '+foptim)

print('Finished...', __file__)
#####################