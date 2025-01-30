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
os.makedirs('list', exist_ok=True)
os.makedirs('fig', exist_ok=True)
os.makedirs('log', exist_ok=True)
os.makedirs('job', exist_ok=True)

optim={
    'Program path'          :'comment',
    'optimroot'             :'/mnt/c/Users/90541/Documents/work/software/optim', #[IN] software path
    'Forward part'          :'comment',
    'ncores'                :8,                         #[IN] number of cores you can apply for modelling
    'fsrc'                  :'list/src.csv',            #[IN] event list file
    'fmodfd'                :'model/inv/inv.vp',        #[IN] current model for forward modelling
    'nxyfd'                 :[680, 155],                #[IN] X- and Y-dir. grid number of forward modeling
    'oxyfd'                 :[0.0, 0.0],                #[IN] X- and Y-dir. original coordinate of model in meter
    'dxyfd'                 :[25.0, 25.0],              #[IN] X- and Y-dir. spacing in meter
    'incxy'                 :[2, 2],                    #[IN] X- and Y-dir. grid sampling for sensitivity matrix
    'Descent part'          :'comment',
    'fmodfdref'             :'model/inv/ref.vp',        #[IN] reference model in forward modelling gird
    'fwgtd'                 :'list/wgtd.npz',           #[IN] data weighting
    'fwgtm'                 :'list/wgtm.npz',           #[IN] model weighting
    'modeltype'             :3,                         #[IN] 1:smoothest, 2:flattest, 3:smallest perturbation
    'descent_type'          :1,                         #[IN] desecent with format of absolute value(1), percent vaule(2)
    'precond'               :1,                         #[IN] preconditional operator for conjugate gradient.
    'niter_inner_max'       :100,                       #[IN] maximum number of lsqr (or cg)
    'lambda0'               :156.25,                    #[IN] lambda for regularization
    'Update part'           :'comment',
    'eps_scale'             :0.02,                      #[IN] step length in proportion
    'smooth_size'           :[150.0, 150.0],            #[IN] smooth size for descent direction, in meter (such as, 6*dxg:6*dyg)
    'mod_limits'            :1,                         #[IN] whether to check the bound of model value, 1: yes, 0: no
    'mod_lb'                :1500,                      #[IN] lower limit of model value if mod_limits=1
    'mod_ub'                :4800,                      #[IN] upper limit of model value if mod_limits=1
    'Inversion loop'        :'comment',
    'abort'                 :0,                         #[IN] intial stage phase
    'iter0'                 :1,                         #[IN] the initial iteration at current inversion stage/phase
    'niter_min'             :20,                        #[IN] minimum iterations
    'niter_max'             :30,                        #[IN] maximum iterations, for future design
    'pro'                   :0.01,                      #[IN] (Ei2 - Ei) / Ei2
    'pro2'                  :0.01,                      #[IN] k2 / k1 
    'pro3'                  :0.01,                      #[IN] average2 / average1 
}
optim['demoroot'] = os.path.dirname(os.path.dirname(__file__)) #.../firstarrival
optim['prjroot'] = os.getcwd() # project path
optim['fdir'] = fdir
optim['dim'] = 2  #[IN] dimension for current study
optim['fmod'] = 'model/inv/inv_tomo.slowness'       #[OUT] current model [x] for tomography inversion
optim['fmodref'] = 'model/inv/ref_tomo.slowness'    #[OUT] reference model mref, just define the file name
optim['fcost'] = 'list/fcost.dat'                   #[OUT] misfit function
optim['fdiff'] = 'list/diff.dat'                    #[OUT] differnece between observed and synthetic data for each pairs
optim['fsensmat'] = 'list/sensmat.npz'              #[OUT] sensitivity matrix [G]
# grids
nx, ny      =optim['nxyfd']
dx, dy      =optim['dxyfd']
incx, incy  =optim['incxy']
optim['n2'] = len(np.arange(nx)[::incx])
optim['n1'] = len(np.arange(ny)[::incy])
optim['d2'] = dx * incx
optim['d1'] = dy * incy
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