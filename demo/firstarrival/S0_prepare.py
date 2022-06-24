#!/usr/bin/env python3
"""
Prepare configure parameters for tomography
Usage:
    python S0_prepare.py optim
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
    'prjroot'               :'/home/tao/Nutstore_Files/works/seisplatform/mod_OPTIM/demo/firstarrival/prj', #[IN] project path
    'Forward part'          :'comment',
    'fmodfd'                :'model/inv/inv.vp',        #[IN] current model for forward modelling
    'n2g'                   :680,                       #[IN] X- grid number of forward modeling              
    'n1g'                   :155,                       #[IN] Y- grid number of forward modeling
    'o2g'                   :0.0,                       #[IN] X- original coordinate of model in meter
    'o1g'                   :0.0,                       #[IN] Y- original coordinate of model in meter
    'd2g'                   :25.0,                      #[IN] X- spacing in meter
    'd1g'                   :25.0,                      #[IN] Y- spacing in meter
    'ncores'                :8,                         #[IN] number of cores you can apply for modelling
    'Diff sensmat part'     :'comment',
    'inc2'                  :2,                         #[IN] X- grid sampling for sensitivity matrix
    'inc1'                  :2,                         #[IN] Y- grid sampling for sensitivity matrix
    'fevn_info'             :'list/eventinfo.csv',      #[IN] event list file
    'fcost'                 :'list/fcost.dat',          #[OUT] misfit function
    'fdiff'                 :'list/diff.dat',           #[OUT] differnece between observed and synthetic data for each pairs
    'fsensmat'              :'list/sensmat.npz',        #[OUT] sensitivity matrix [G]
    'Descent part lsqr'     :'comment',
    'dim'                   :2,                         #[IN] dimension for current study
    'modeltype'             :3,                         #[IN] 1:smoothest, 2:flattest, 3:smallest perturbation
    'descent_type'          :1,                         #[IN] desecent with format of absolute value(1), percent vaule(2)
    'fmodfdref'             :'model/inv/ref.vp',        #[IN] reference model in forward modelling gird
    'fmod'                  :'model/inv/inv_tomo.slowness',   #[OUT] current model [x] for tomography inversion
    'fmodref'               :'model/inv/ref_tomo.slowness',   #[OUT] reference model mref, just define the file name
    'fwgtd'                 :'list/wgtd.npz',           #[IN] data weighting
    'fwgtm'                 :'list/wgtm.npz',           #[IN] model weighting
    'niter_inner_max'       :100,                       #[IN] maximum number of lsqr
    'lambda0'               :10000.0,                   #[IN] initial lambda for L-curve test, suggest which is not more than 1.0E18
    'step'                  :2.0,                       #[IN] the ratio between the neighbor damping factor
    'testmax'               :20,                        #[IN] the maximum test times
    'Update part'           :'comment',
    'eps_scale'             :0.02,                      #[IN] step length in proportion
    'smooth_size'           :'150.0:150.0',             #[IN] smooth size for descent direction, in meter (such as, 6*dxg:6*dyg)
    'mod_limits'            :1,                         #[IN] 1: check the bound of model value;0 no limit for model value
    'mod_lb'                :1500,                      #[IN] lower limit of model value in m/s
    'mod_ub'                :4800,                      #[IN] upper limit of model value in m/s
    'Inversion loop'        :'comment',
    'iter0'                 :1,                         #[IN] the initial iteration at current inversion stage/phase
    'niter_max'             :30,                        #[IN] maximum iterations, for future design
}
optim['fdir'] = fdir
optim['n1'] = len(np.arange(optim['n1g'])[::optim['inc1']])
optim['n2'] = len(np.arange(optim['n2g'])[::optim['inc2']])
optim['d1'] = optim['d1g'] * optim['inc1']
optim['d2'] = optim['d2g'] * optim['inc2']
if optim['modeltype'] != 2:
    optim['numreg'] = 1
else:
    optim['numreg'] = optim['dim']
# write work infomation
foptim = os.path.join(fdir,'optim.json')
fout = open(foptim,'w')
fout.write(json.dumps(optim,indent=4))
fout.close()
print('*'*60+'\n*Output file: '+foptim)

print('Finished...', __file__)