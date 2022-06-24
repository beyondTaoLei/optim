#!/usr/bin/env python3
"""
generate initial model and input json file for Rosenbrock problem
Usage:
    python S0_prepare.py optim_PSTD
    python S0_prepare.py optim_PNLCG
    python S0_prepare.py optim_LBFGS
    python S0_prepare.py optim_TRN
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
    'fmod'              :'FD/mod.bin',
    'fcost'             :'FD/misfit.dat',
    'fgrad'             :'FD/grad.bin',
    'fhess'             :'FD/hess.bin',
    'fvctr'             :'FD/vctr.bin',
    'iter0'             :1,         #the initial iteration at current inversion stage/phase
    'niter_max'         :200,       #maximum iterations, for future design
    'conv'              :0.1,       #PRO, for future design
    'steplength'    :'comment',
    'try_old_sl'        :1,         #whether to use steplength from last iter.(1-yes, 0-no)
    'eps_scale'         :0.1,       #max_m0*eps_scale as one unit to test the steplength
    'max_m0'            :0.25,      #the maximum of initial model
    'alpha_lb'          :0.01,      #lower limit of steplength
    'alpha_ub'          :5,         #upper limit of steplength
    'PNLCG'         :'comment',
    'powell'            :0,         #wheter to check Powell restart condition(1-yes, 0-no)
    'l-BFGS'        :'comment',
    'l'                 :5,         #the maximum storing number of the pairs {sk,yk}(~3-7)
    'TN'            :'comment',     
    'niter_max_CG'      :12,        #maximum number of inner conjugate gradient loop
    'flag_singularity'  :1,         #singularity test, please check the code
    'flag_negative'     :2,         #negative curvature test, please check the code
    'flag_truncation'   :1,         #truncation test, please check the code
    'flag_min_res'      :1,         #record minimum residual when applying the inner CG agrithm
    'pro_sing'          :1.0E-10,   #for singularity test
    '#pro_curv'         :1.0E-10,   #no useful
    'pro_trun'          :0.5,       #for truncation test
    'eta'               :0.9,       #check the code
    'mod_limits'        :0,         #1: check the bound of model value;0 no limit for model value
    'mod_lb'            :1000,      #lower limit of model value
    'mod_ub'            :5000       #upper limit of model value
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