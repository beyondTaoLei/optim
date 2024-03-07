#!/usr/bin/env python3
"""
generate initial model and input json file for Rosenbrock problem
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
    'submit'                :1,                         #[IN] 1: run inverison; 0: only generate inversion commands
    'optimroot'             :'/home/tao/Nutstore_Files/works/optim', #[IN] software path
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
    'n3g'                   :1,                         #[IN] X- grid number of input model
    'n2g'                   :1300,                      #[IN] Y- grid number of input model
    'n1g'                   :400,                       #[IN] Z- grid number of input model
    'o3g'                   :0.0,                       #[IN] X- original coordinate of input model in meter
    'o2g'                   :0.0,                       #[IN] Y- original coordinate of input model in meter
    'o1g'                   :0.0,                       #[IN] Z- original coordinate of input model in meter
    'd3g'                   :500.0,                     #[IN] X- spacing of input model in meter
    'd2g'                   :500.0,                     #[IN] Y- spacing of input model in meter
    'd1g'                   :500.0,                     #[IN] Z- spacing of input model in meter
    'Diff sensmat part'     :'comment',
    'odir'                  :'obs',                     #[IN] work directory for observed data
    'inc3'                  :1,                         #[IN] X- grid sampling for sensitivity matrix, kernels are mapped into grids [::inc3, ::inc2,::inc1]
    'inc2'                  :10,                        #[IN] Y- grid sampling for sensitivity matrix
    'inc1'                  :1,                         #[IN] Z- grid sampling for sensitivity matrix
    'Optimization'          :'comment',
    'fcost'                 :'list/fcost.dat',          #[OUT] the misfit at point fmod 
    'fdiff'                 :'list/diff.dat',           #[OUT] differnece between observed and synthetic data for each pairs
    'fsensmat'              :'list/sensmat.npz',        #[OUT] sensitivity matrix [G]
    'Descent lsqr or plcg'  :'comment',
    'modeltype'             :3,                         #[IN] 1:smoothest, 2:flattest, 3:smallest perturbation
    'descent_type'          :1,                         #[IN] desecent with format of absolute value(1), percent vaule(2)
    'fmodfdref'             :'model/inv/ref.vs',        #[IN] reference model in input model grid
    'fmod'                  :'model/inv/inv_tomo.vs',   #[OUT] current model [x] for tomography inversion
    'fmodref'               :'model/inv/ref_tomo.vs',   #[OUT] reference model mref, just define the file name
    'fwgtd'                 :'list/wgtd.npz',           #[IN] data weighting
    'fwgtm'                 :'list/wgtm.npz',           #[IN] model weighting
    'percond'               :0,                         #[IN] preconditional operator for conjugate gradient.
    'niter_inner_max'       :100,                       #[IN] maximum number of lsqr
    'lambda0'               :100.0,                     #[IN] initial lambda for L-curve test, suggest which is not more than 1.0E18
    'step'                  :4.0,                       #[IN] the ratio between the neighbor damping factor
    'testmax'               :20,                        #[IN] the maximum test times
    'steplength'            :'comment',
    'eps_scale'             :0.03,                      #[IN] step length in proportion
    'smooth_size'           :'0.0:30000.0:9000.0',      #[IN] smooth size for descent direction, in meter (such as, 6*d3g:6*d2g:6*d1g)
    'mod_limits'            :1,                         #[IN] 1: check the bound of model value;0 no limit for model value
    'mod_lb'                :2500,                      #[IN] lower limit of model value in m/s
    'mod_ub'                :5000,                      #[IN] upper limit of model value in m/s
    'Inversion loop'        :'comment',
    'iter0'                 :1,                         #[IN] the initial iteration at current inversion stage/phase
    'niter_max'             :30,                        #[IN] maximum iterations, for future design
}
optim['prjroot'] = os.getcwd() # project path
optim['fdir'] = fdir
optim['n3'] = len(np.arange(optim['n3g'])[::optim['inc3']])
optim['n2'] = len(np.arange(optim['n2g'])[::optim['inc2']])
optim['n1'] = len(np.arange(optim['n1g'])[::optim['inc1']])
optim['d3'] = optim['d3g'] * optim['inc3']
optim['d2'] = optim['d2g'] * optim['inc2']
optim['d1'] = optim['d1g'] * optim['inc1']
if optim['modeltype'] != 2:
    optim['numreg'] = 1
else:
    optim['numreg'] = 2 # dimensions

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