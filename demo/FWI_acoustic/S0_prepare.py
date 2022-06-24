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
    'submit'                :1, #1: only generate the lsf command, 0: run inversion workflow
    'platform'              :'lsf', #lsf, 
    'ncores'                :40,                #how many CPU cores in one node
    'Path definition'       :'comment',
    'prjroot'               :'/data/ess-leit/Marmousi2_new',    #change the name according to your project
    'optimroot'             :'/work/ess-leit/seisplatform/mod_OPTIM',
    'mpiexec'               :'/share/intel/2018u4/compilers_and_libraries_2018.5.274/linux/mpi/intel64/bin/mpiexec',
    'FDEXE'                 :'/work/ess-leit/seisplatform/bin/sofi2Dac',
    'ADEXE'                 :'/work/ess-leit/seisplatform/bin/sofi2Dac_grad',
    'Wave simulation'       :'comment',
    'NPROCX'                :4,
    'NPROCY'                :1,
    'NX'                    :1700,
    'NY'                    :350,
    'DH'                    :10.0,
    'FREE_SURF'             :1,
    'FW'                    :15,
    'IDX'                   :1,
    'IDY'                   :1,
    'fevn_info'             :'list/eventinfo.csv',
    'Frequency filtering'   :'comment',
    'fc'                    :2.5,       #low cut-off frequency
    'order'                 :6,         #the order of the Butterworth lowpass filter
    'zero_phase'            :0,         #0, forward; 1, filter trace in both forward and reverse direciton 
    "Misfit definition"     :'comment',
    'lnorm'                 :1,         #1: L2 norm
    'Gradient'              :'comment',
    'taper'                 :1,
    'ftaper'                :'taper/taper_all.bin',
    'smooth_size'           :5,
    'eps_scale'             :0.002,
    'Optimization'          :'comment',
    'fdio'                  :'fdio',
    'fmod'                  :'model/inv/inv.vp',
    'fcost'                 :'fdio/comm/misfit.dat',
    'fgrad'                 :'fdio/comm/grad.bin',
    'fhess'                 :'fdio/comm/hess.bin',
    'fvctr'                 :'fdio/comm/vctr.bin',
    'iter0'                 :1,         #the initial iteration at current inversion stage/phase
    'niter_max'             :30,        #maximum iterations, for future design
    'conv'                  :0.1,       #PRO, for future design
    'Steplength'            :'comment',
    'try_old_sl'            :1,         #whether to use steplength from last iter.(1-yes, 0-no)
    'eps_scale'             :0.002,     #max_m0*eps_scale as one unit to test the steplength
    'max_m0'                :4700.0,    #the maximum of initial model
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
    'mod_lb'                :1000,      #lower limit of model value
    'mod_ub'                :5000       #upper limit of model value
}
mod = np.fromfile(optim['fmod'], np.float32)
optim['n'] = mod.shape[0] #only for TRN now
optim['fdir'] = fdir
optim['NP'] = optim['NPROCX']*optim['NPROCY']
# write work infomation
foptim = os.path.join(fdir,'optim.json')
fout = open(foptim,'w')
fout.write(json.dumps(optim,indent=4))
fout.close()
print("*"*60+"\n*Output file: "+foptim)

print("Finished...", __file__)
#####################