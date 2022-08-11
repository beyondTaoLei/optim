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
    'fmod'              :'FD/mod.bin',      #[IN] the current model in work memory
    'fcost'             :'FD/misfit.dat',   #[OUT] the misfit at point fmod
    'fgrad'             :'FD/grad.bin',     #[OUT] the gradient at point fmod
    'fhess'             :'FD/hess.bin',     #[OUT] the Hessian vector product at point [fmod, fvctr]
    'fvctr'             :'FD/vctr.bin',     #[OUT] the vector pk in LCG algorithm
    'iter0'             :1,                 #[IN] the initial iteration at current inversion stage/phase
    'niter_max'         :200,               #[IN] maximum iterations, for future design
    'conv'              :0.1,               #[IN] PRO, for future design
    'steplength'    :'comment',
    'try_old_sl'        :1,                 #[IN] whether to use steplength from last iter.(1-yes, 0-no)
    'eps_scale'         :0.1,               #[IN] max_m0*eps_scale as one unit to test the steplength
    'max_m0'            :0.25,              #[IN] the maximum of initial model
    'alpha_lb'          :0.01,              #[IN] lower limit of steplength
    'alpha_ub'          :5,                 #[IN] upper limit of steplength
    'mod_limits'        :0,                 #[IN] 1: check the bound of model value;0 no limit for model value
    'mod_lb'            :1000,              #[IN] lower limit of model value
    'mod_ub'            :5000,              #[IN] upper limit of model value
    'PNLCG'         :'comment',
    'powell'            :0,                 #[IN] wheter to check Powell restart condition(1-yes, 0-no)
    'l-BFGS'        :'comment',
    'l'                 :5,                 #[IN] the maximum storing number of the pairs {sk,yk}(~3-7)
    'TN'            :'comment',     
    'niter_max_CG'      :12,                #[IN] maximum number of inner conjugate gradient loop
    'flag_singularity'  :1,                 #[IN] singularity test, please check the code
    'flag_negative'     :2,                 #[IN] negative curvature test, please check the code
    'flag_truncation'   :1,                 #[IN] truncation test, please check the code
    'flag_min_res'      :1,                 #[IN] record minimum residual when applying the inner CG agrithm
    'pro_sing'          :1.0E-10,           #[IN] for singularity test
    'pro_curv'          :1.0E-10,           #[IN] for negative curvature test (if flag_negative == 1)
    'pro_trun'          :0.5,               #[IN] for truncation test
    'eta'               :0.9,               #[IN] check the code
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