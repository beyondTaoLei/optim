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
fdir    =sys.argv[1] # optim
# make directory
os.makedirs(fdir, exist_ok=True)
os.makedirs('fig', exist_ok=True)
os.makedirs('log', exist_ok=True)
os.makedirs('job', exist_ok=True)

optim={
    'submit'                :1,                         #[IN] 1: run inverison; 0: only generate inversion commands
    'optimroot'             :'/home/tao/Nutstore_Files/works/seisplatform/mod_OPTIM', #[IN] software path
    'prjroot'               :'/home/tao/Nutstore_Files/works/seisplatform/mod_OPTIM/demo/surface_phv_lsqr/prj',
    'Forward part'          :'comment',
    'fmodfd'                :'model/inv/inv.vs',        #[IN] current model for forward modelling
    'fperiods'              :'list/periods.csv',        #[IN] measured periods in second, starting with low periods
    'algo'                  :'dunkin',                  #[IN] algorithm for disperison calculation
    'wavetype'              :'rayleigh',                #[IN] 'love' or 'rayleigh'
    'maxmode'               :0,                         #[IN] code, fund. 0, first high mode 1
    'wdir'                  :'syn',                     #[IN] work directory for synthetic data
    'n2g'                   :1300,                      #[IN] X- grid number of forward modeling 
    'n1g'                   :400,                       #[IN] Y- grid number of forward modeling
    'o2g'                   :0.0,                       #[IN] X- original coordinate of model in meter
    'o1g'                   :0.0,                       #[IN] Y- original coordinate of model in meter
    'd2g'                   :500.0,                     #[IN] X- spacing in meter
    'd1g'                   :500.0,                     #[IN] Y- spacing in meter
    'Diff sensmat part'     :'comment',
    'odir'                  :'obs',                     #[IN] work directory for observed data
    'fcost'                 :'list/fcost.dat',          #[OUT] misfit function
    'fdiff'                 :'list/diff.dat',           #[OUT] differnece between observed and synthetic data for each pairs
    'fsensmat'              :'list/sensmat.npz',        #[OUT] sensitivity matrix [G]
    'Descent part lsqr'     :'comment',
    'dim'                   :2,                         #[IN] dimension for current study
    'modeltype'             :3,                         #[IN] 1:smoothest, 2:flattest, 3:smallest perturbation
    'descent_type'          :1,                         #[IN] desecent with format of absolute value(1), percent vaule(2)
    'fmodref'               :'model/inv/ref.vs',        #[IN] reference model in forward modelling gird
    'fwgtd'                 :'list/wgtd.npz',           #[IN] data weighting
    'fwgtm'                 :'list/wgtm.npz',           #[IN] model weighting
    'niter_lsqr'            :100,                       #[IN] maximum number of lsqr
    'lambda0'               :100.0,                     #[IN] initial lambda for L-curve test, suggest which is not more than 1.0E18
    'step'                  :4.0,                       #[IN] the ratio between the neighbor damping factor
    'testmax'               :20,                        #[IN] the maximum test times
    'Update part'           :'comment',
    'eps_scale'             :0.03,                      #[IN] step length in proportion
    'smooth_size'           :'3000.0:3000.0',           #[IN] smooth size for descent direction, in meter (such as, 6*dxg:6*dyg)
    'mod_limits'            :1,                         #[IN] 1: check the bound of model value;0 no limit for model value
    'mod_lb'                :2500,                      #[IN] lower limit of model value in m/s
    'mod_ub'                :5000,                      #[IN] upper limit of model value in m/s
    'Inversion loop'        :'comment',
    'iter0'                 :1,                         #[IN] the initial iteration at current inversion stage/phase
    'niter_max'             :30,                        #[IN] maximum iterations, for future design
}
optim['fdir'] = fdir
optim['fmod'] = optim['fmodfd'] # grid of FD = grid of inversion
optim['n2'] = optim['n2g']
optim['n1'] = optim['n1g']
optim['d1'] = optim['d1g']
optim['d2'] = optim['d2g']
# write work infomation
foptim = os.path.join(fdir,'optim.json')
fout = open(foptim,'w')
fout.write(json.dumps(optim,indent=4))
fout.close()
print('*'*60+'\n*Output file: '+foptim)

print('Finished...', __file__)
#####################