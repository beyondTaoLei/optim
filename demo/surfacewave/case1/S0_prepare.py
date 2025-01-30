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
    'incxyz'                :[1, 10, 1],                #[IN] X-, Y- and Z-dir grid sampling for sensitivity matrix
    'Descent part'          :'comment',
    'fmodfdref'             :'model/inv/ref.vs',        #[IN] reference model in input model grid
    'fwgtd'                 :'list/wgtd.npz',           #[IN] data weighting
    'fwgtm'                 :'list/wgtm.npz',           #[IN] model weighting
    'modeltype'             :3,                         #[IN] 1:smoothest, 2:flattest, 3:smallest perturbation
    'descent_type'          :1,                         #[IN] desecent with format of absolute value(1), percent vaule(2)
    'precond'               :0,                         #[IN] preconditional operator for conjugate gradient.
    'niter_inner_max'       :100,                       #[IN] maximum number of lsqr
    'lambda0'               :6.0e-3,                    #[IN] initial lambda for L-curve test, suggest which is not more than 1.0E18
    'Update part'           :'comment',
    'eps_scale'             :0.03,                      #[IN] step length in proportion
    'smooth_size'           :[0.0, 30000.0, 9000.0],    #[IN] smooth size for descent direction, in meter (such as, 6*d3g:6*d2g:6*d1g)
    'mod_limits'            :1,                         #[IN] whether to check the bound of model value, 1: yes, 0: no
    'mod_lb'                :2500,                      #[IN] lower limit of model value if mod_limits=1
    'mod_ub'                :5000,                      #[IN] upper limit of model value if mod_limits=1
    'Inversion loop'        :'comment',
    'abort'                 :0,                         #[IN] intial stage phase
    'iter0'                 :1,                         #[IN] the initial iteration at current inversion stage/phase
    'niter_min'             :20,                        #[IN] minimum iterations
    'niter_max'             :30,                        #[IN] maximum iterations, for future design
    'pro'                   :0.01,                      #[IN] (Ei2 - Ei) / Ei2
    'pro2'                  :0.01,                      #[IN] k2 / k1 
    'pro3'                  :0.01,                      #[IN] average2 / average1 
}

optim['demoroot'] = os.path.dirname(os.path.dirname(__file__)) #.../surfacewave
optim['prjroot'] = os.getcwd() # project path
optim['fdir'] = fdir
optim['odir'] = 'obs'                               #[IN] work directory for observed data
optim['wdir'] = 'syn'                               #[IN] work directory for synthetic data
optim['fmod'] = 'model/inv/inv_tomo.vs'             #[OUT] current model [x] for tomography inversion
optim['fmodref'] = 'model/inv/ref_tomo.vs'          #[OUT] reference model mref, just define the file name
optim['fcost'] = 'list/fcost.dat'                   #[OUT] misfit function
optim['fdiff'] = 'list/diff.dat'                    #[OUT] differnece between observed and synthetic data for each pairs
optim['fsensmat'] = 'list/sensmat.npz'              #[OUT] sensitivity matrix [G]
# inverse grids
optim['n3'] = len(np.arange(optim['nxyzfd'][0])[::optim['incxyz'][0]])
optim['n2'] = len(np.arange(optim['nxyzfd'][1])[::optim['incxyz'][1]])
optim['n1'] = len(np.arange(optim['nxyzfd'][2])[::optim['incxyz'][2]])
optim['d3'] = optim['dxyzfd'][0] * optim['incxyz'][0]
optim['d2'] = optim['dxyzfd'][1] * optim['incxyz'][1]
optim['d1'] = optim['dxyzfd'][2] * optim['incxyz'][2]
if optim['modeltype'] != 2:
    optim['numreg'] = 1
else:
  if optim['nxyzfd'][0] == 1:
    optim['numreg'] = 2 # dimensions
  else:
    optim['numreg'] = 3 # dimensions

# write work infomation
foptim = os.path.join(fdir,'optim.json')
fout = open(foptim,'w')
fout.write(json.dumps(optim,indent=4))
fout.close()
print('*'*60+'\n*Output file: '+foptim)

print('Finished...', __file__)
#####################