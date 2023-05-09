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

# Input paras
fdir    =sys.argv[1] # optim_PSTD, optim_PNLCG, optim_LBFGS, optim_TRN 
# make directory
os.makedirs(fdir, exist_ok=True)
os.makedirs('fig', exist_ok=True)
os.makedirs('job', exist_ok=True)
os.makedirs('log', exist_ok=True)
try:
    os.remove(os.path.join(fdir, 'step.json'))
except OSError:
    pass

optim={
    "fperiods"              :"list/periods.csv",  #in second, sorted starting with low periods
    "fobstime"              :"list/obs_time.csv",
    "algo"                  :"dunkin",
    "wavetype"              :"rayleigh",            #'love' or 'rayleigh'
    "maxmode"               :0,                     #code, fund. 0, first high mode 1
    "wdir"                  :"step1",
    "ox"                    :0.0,                   #in meter
    "oz"                    :0.0,                   #in meter
    "dx"                    :200.0,                 #in meter
    "dz"                    :100.0,                 #in meter
    "nx"                    :120,
    "nz"                    :60,
    "mx"                    :104,
    "mz"                    :123,
    "nxs"                   :30,
    "nxe"                   :89,
    "nzs"                   :0,
    "nze"                   :59,
    "mxs"                   :22,
    "mxe"                   :81,
    "mzs"                   :0,
    "mze"                   :59,
    "pertb"                 :0.001,
    "Gradient part"         :"comment",
    'smooth_size'           :900.0,                #in meter(6*dx)
    'Optimization'          :'comment',
    'fmod'                  :'model/inv/inv.vs',    #vp/vs/rho: in m/s,m/s,kg/m3
    'fmodref'               :'model/inv/inv_ref.vs',   #[IN] reference model mref, just define the file name
    'fmod2'                 :'model/inv/inv.p',    # electrical resistivity
    'fcost'                 :'step1/misfit.dat',
    'fgrad'                 :'step1/grad.bin',
    'iter0'                 :1,         #the initial iteration at current inversion stage/phase
    'niter_max'             :30,        #maximum iterations, for future design
    'conv'                  :0.1,       #PRO, for future design
    'Steplength'            :'comment',
    'try_old_sl'            :1,         #whether to use steplength from last iter.(1-yes, 0-no)
    'eps_scale'             :0.04,     #max_m0*eps_scale as one unit to test the steplength
    'max_m0'                :3500.0,    #the maximum of initial model
    'alpha_lb'              :0.1,       #lower limit of steplength
    'alpha_ub'              :10,        #upper limit of steplength
    'PNLCG'                 :'comment',
    'powell'                :1,         #wheter to check Powell restart condition(1-yes, 0-no)
    'LBFGS'                 :'comment',
    'l'                     :5,         #the maximum storing number of the pairs {sk,yk}(~3-7)
    'TRN'                   :'comment',     
    'niter_max_CG'          :12,        #maximum number of inner conjugate gradient loop
    'flag_singularity'      :1,         #singularity test, please check the code
    'flag_negative'         :2,         #negative curvature test, please check the code
    'flag_truncation'       :1,         #truncation test, please check the code
    'flag_min_res'          :1,         #record minimum residual when applying the inner CG agrithm
    'pro_sing'              :1.0E-10,   #for singularity test
    '#pro_curv'             :1.0E-10,   #no useful
    'pro_trun'              :0.5,       #for truncation test
    'eta'                   :0.9,       #check the code
    'mod_limits'            :1,         #1: check the bound of model value;0 no limit for model value
    'mod_lb'                :1200,      #lower limit of model value
    'mod_ub'                :4000       #upper limit of model value
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

foptim = os.path.join(fdir,'mod_info.dat')
fp=open(foptim, 'w')
fp.write('%d %d\n' % (optim['nx'], optim['nz']))
fp.write('%d %d\n' % (optim['mod_lb'], optim['mod_ub']))
fp.close()

print("Finished...", __file__)
#####################
