#!/usr/bin/env python3
"""
generate initial model and input json file for Rosenbrock problem
Usage:
    cd /home/tao/Nutstore_Files/works/optim/demo/rosenbrock/prj
    python ../S0_prepare.py optim_PSTD
    python ../S0_prepare.py optim_PNLCG
    python ../S0_prepare.py optim_LBFGS
    python ../S0_prepare.py optim_TRN
"""

import os
import sys
import json
import numpy as np

#Input paras
fdir    =sys.argv[1] # optim_PSTD, optim_PNLCG, optim_LBFGS, optim_TRN 
# make directory
os.makedirs(fdir, exist_ok=True)
os.makedirs('FD', exist_ok=True)
os.makedirs('fig', exist_ok=True)
try:
    os.remove(os.path.join(fdir, 'step.json'))
except OSError:
    pass

# initial model
a=np.array([0.25, 0.25], np.float32)
a.tofile('FD/mod.bin')

optim={
    'optimroot'         :'/home/tao/Nutstore_Files/works/optim', #[IN] software path
    'demoroot'          :'/home/tao/Nutstore_Files/works/optim/demo/rosenbrock', #[IN] current demo path
    'Optimization'      :'comment',
    'fmod'              :'FD/mod.bin',      #[IN] the current model in work memory
    'fcost'             :'FD/misfit.dat',   #[OUT] the misfit at point fmod
    'fgrad'             :'FD/grad.bin',     #[OUT] the gradient at point fmod
    'fhess'             :'FD/hess.bin',     #[OUT] the Hessian vector product at point [fmod, fvctr]
    'fvctr'             :'FD/vctr.bin',     #[OUT] the vector pk in LCG algorithm
    'PNLCG'             :'comment',
    'powell'            :0,                 #[IN] wheter to check Powell restart condition(1-yes, 0-no)
    'l-BFGS'            :'comment',
    'l'                 :5,                 #[IN] the maximum storing number of the pairs {sk,yk}(~3-7)
    'TN'                :'comment',
    'conv_CG'           :0,                 #[IN] initial status for debug, set to 0
    'niter_max_CG'      :12,                #[IN] maximum number of inner conjugate gradient loop
    'flag_singularity'  :1,                 #[IN] 1: singularity test; 0: no test
    'flag_negative'     :2,                 #[IN] 1: negative curvature test; 2: strong negative curvature test; 0: no test
    'flag_truncation'   :1,                 #[IN] 1: truncation test based on residuals; 2: test based on parabolic function, not work now; 0: no test
    'flag_min_res'      :1,                 #[IN] 1: record minimum residual when applying the inner CG agrithm; 0: no record
    'pro_sing'          :1.0E-10,           #[IN] threshold \varsigma for singularity test
    'pro_curv'          :1.0E-10,           #[IN] threshold \delta for negative curvature test (flag_negative = 1)
    'pro_trun'          :0.5,               #[IN] threshold c_q for truncation test
    'eta'               :0.9,               #[IN] if ||residuals|| <= eta * ||gradient||, skip the inner CG loop
    'eta_update'        :0,                 #[IN] 1: update the eta during inversion cycle; 0: not update
    'steplength'        :'comment',
    'try_old_sl'        :1,                 #[IN] whether to use steplength from last iter.(1-yes, 0-no)
    'eps_scale'         :0.1,               #[IN] max_m0*eps_scale as one unit to test the steplength
    'max_m0'            :0.25,              #[IN] the maximum of initial model
    'alpha_lb'          :0.01,              #[IN] lower limit of steplength
    'alpha_ub'          :5,                 #[IN] upper limit of steplength
    'mod_limits'        :0,                 #[IN] 1: check the bound of model value;0 no limit for model value
    'mod_lb'            :1000,              #[IN] lower limit of model value
    'mod_ub'            :5000,              #[IN] upper limit of model value
    'Inversion loop'    :'comment',
    'iter0'             :1,                 #[IN] the initial iteration at current inversion stage/phase
    'niter_max'         :200,               #[IN] maximum iterations, for future design
}
mod = np.fromfile(optim['fmod'], np.float32)
optim['n'] = mod.shape[0] #only for TRN now
optim['fdir'] = fdir
# write work infomation
foptim = os.path.join(fdir,'optim.json')
fout = open(foptim,'w')
fout.write(json.dumps(optim,indent=4))
fout.close()
print("*"*60+"\n*Output file: "+foptim)

print("Finished...", __file__)
#####################