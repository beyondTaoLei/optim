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
import pandas as pd

# Input paras
fdir    =sys.argv[1] # optim
# make directory
os.makedirs(fdir, exist_ok=True)
os.makedirs('fig', exist_ok=True)
os.makedirs('log', exist_ok=True)
os.makedirs('job', exist_ok=True)

optim={
    'submit'                :1,                         #[IN] 1: run inverison; 0: only generate inversion commands
    'optimroot'             :'/home/tao/Nutstore_private/works/software/optim', #[IN] software path
    'Forward part'          :'comment',
    'mpiexec'               :'mpirun',
    'ncores'                :8,                         #[IN] number of cores you can apply for modelling
    'fcmp'                  :'list/cmp.csv',            #[IN] the measure coordinates (CMP)
    'fmodfd'                :'model/inv/inv.vs',        #[IN] the current model in work memory
    'fperiods'              :'list/periods.csv',        #[IN] measured periods in second, starting with low periods
    'algo'                  :'dunkin',                  #[IN] algorithm for disperison calculation
    'wavetype'              :'rayleigh',                #[IN] 'love' or 'rayleigh'
    'maxmode'               :0,                         #[IN] code, fund. 0, first high mode 1
    'wdir'                  :'syn',                     #[IN] work directory for synthetic data
    'n3g'                   :1,                         #[IN] X- grid number of forward modeling
    'n2g'                   :1400,                      #[IN] Y- grid number of forward modeling 
    'n1g'                   :100,                       #[IN] Z- grid number of forward modeling
    'o3g'                   :0.0,                       #[IN] X- original coordinate of model in meter
    'o2g'                   :0.0,                       #[IN] Y- original coordinate of model in meter
    'o1g'                   :0.0,                       #[IN] Z- original coordinate of model in meter
    'd3g'                   :100.0,                     #[IN] X- spacing in meter
    'd2g'                   :100.0,                     #[IN] Y- spacing in meter
    'd1g'                   :100.0,                     #[IN] Z- spacing in meter
    'Gradient part'         :'comment',
    'odir'                  :'obs',                     #[IN] work directory for observed data
    'fcost'                 :'list/fcost.dat',          #[OUT] the misfit at point fmod
    'fgrad'                 :'list/grad.bin',           #[OUT] the gradient at point fmod
    'fhess'                 :'list/hess.bin',           #[OUT] the Hessian vector product at point [fmod, fvctr]
    'fvctr'                 :'list/vctr.bin',           #[OUT] the vector pk in LCG algorithm
    'smooth_size'           :'0.0:1000.0:500.0',        #[IN] smooth size for descent direction, in meter (such as, 6*d3g:6*d2g:6*d1g)
    'percond'               :1,                         #[IN] preconditional operator for gradient.
    'EPSILON'               :0.05,                      #[IN] preconditional gradient is grad/(diagHessian+diagHessian.mean(axis=0).max()*EPSILON).
    'Optimization'          :'comment',
    'PNLCG'                 :'comment',
    'powell'                :1,                         #[IN] wheter to check Powell restart condition(1-yes, 0-no)
    'l-BFGS'                :'comment',
    'l'                     :5,                         #[IN] the maximum storing number of the pairs {sk,yk}(~3-7)
    'TN'                    :'comment',
    'conv_CG'               :0,                         #[IN] initial status
    'niter_max_CG'          :12,                        #[IN] maximum number of inner conjugate gradient loop
    'flag_singularity'      :1,                         #[IN] singularity test, please check the code
    'flag_negative'         :2,                         #[IN] negative curvature test, please check the code
    'flag_truncation'       :1,                         #[IN] truncation test, please check the code
    'flag_min_res'          :1,                         #[IN] record minimum residual when applying the inner CG agrithm
    'pro_sing'              :1.0E-10,                   #[IN] for singularity test
    'pro_curv'              :1.0E-10,                   #[IN] for negative curvature test (if flag_negative == 1)
    'pro_trun'              :0.5,                       #[IN] for truncation test
    'eta'                   :0.9,                       #[IN] check the code
    'Steplength'            :'comment',
    'try_old_sl'            :1,                         #[IN] whether to use steplength from last iter.(1-yes, 0-no)
    'eps_scale'             :0.005,                      #[IN] max_m0*eps_scale as one unit to test the steplength
    'max_m0'                :4000.0,                    #[IN] the maximum of initial model
    'alpha_lb'              :0.1,                       #[IN] lower limit of steplength
    'alpha_ub'              :10,                        #[IN] upper limit of steplength
    'mod_limits'            :1,                         #[IN] 1: check the bound of model value;0 no limit for model value
    'mod_lb'                :2000,                      #[IN] lower limit of model value in m/s
    'mod_ub'                :4000,                      #[IN] upper limit of model value in m/s
    'Inversion loop'        :'comment',
    'iter0'                 :1,                         #[IN] the initial iteration at current inversion stage/phase
    'niter_max'             :30,                        #[IN] maximum iterations, for future design
}
optim['prjroot'] = os.getcwd() # project path
optim['fdir'] = fdir
optim['fmod'] = optim['fmodfd'] # grid of FD = grid of inversion
optim['n3'] = optim['n3g']
optim['n2'] = optim['n2g']
optim['n1'] = optim['n1g']
optim['d3'] = optim['d3g']
optim['d2'] = optim['d2g']
optim['d1'] = optim['d1g']
optim['n'] = optim['n3']*optim['n2']*optim['n1'] #only for TRN now

#read indexes of CMP
df = pd.read_csv(optim['fcmp'])
optim['ntraces'] = len(df['offset'])

df = pd.read_csv(optim['fperiods'])
periods = np.array(df['T0'])
optim['nper'] = len(periods)

# write work infomation
foptim = os.path.join(fdir,'optim.json')
fout = open(foptim,'w')
fout.write(json.dumps(optim,indent=4))
fout.close()
print('*'*60+'\n*Output file: '+foptim)

print('Finished...', __file__)
#####################