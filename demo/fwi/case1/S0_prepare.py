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
fdir    =sys.argv[1] # optim_PSTD, optim_PNLCG, optim_LBFGS, optim_TRN 
# make directory
os.makedirs(fdir, exist_ok=True)
os.makedirs('fig', exist_ok=True)
os.makedirs('job', exist_ok=True)
os.makedirs('log', exist_ok=True)
os.makedirs('list', exist_ok=True)
os.makedirs('fdio', exist_ok=True)
os.makedirs('data', exist_ok=True)
os.makedirs('source', exist_ok=True)
try:
    os.remove(os.path.join(fdir, 'step.json'))
except OSError:
    pass

#'/share/intel/2018u4/compilers_and_libraries_2018.5.274/linux/mpi/intel64/bin/mpiexec'
#'/opt/openmpi/bin/mpirun'
#/home/tao/seis_software/seisplatform/FD/SOFI2DAC/bin/sofi2Dac
#/home/tao/seis_software/seisplatform/FD/SOFI2DAC/bin/sofi2Dac_grad
optim={
    'submit'                :1,                     #[IN] 1: run inverison; 0: only generate inversion commands
    'platform'              :'lsf',                 #[IN] remote batch system: so far with [lsf,] 
    'ncores'                :40,                    #[IN] number of CPU cores in one node
    'Path definition'       :'comment',
    'optimroot'             :'/work/ess-leit/optim',            #[IN] software path
    'mpiexec'               :'/share/intel/2018u4/compilers_and_libraries_2018.5.274/linux/mpi/intel64/bin/mpiexec',
    'fdexe'                 :'/work/ess-leit/seisplatform/bin/sofi2Dac',        #[IN] forward modeling operator
    'adexe'                 :'/work/ess-leit/seisplatform/bin/sofi2Dac_grad',   #[IN] adjoint modeling operator
    'Wave simulation'       :'comment',
    'json_test'             :'list/FW.json',        #[IN] json file only for misfit calculation
    'json_grad'             :'list/GRAD.json',      #[IN] json file for gradient calculation
    'fmodfd'                :'model/inv/inv.vp',    #[IN] the current model in work memory
    'nprocx'                :4,                     #[IN] number of MPI processors in X-direction
    'nprocy'                :1,                     #[IN] number of MPI processors in Y-direction
    'n3g'                   :1,                     #[IN] X- grid number of forward modeling
    'n2g'                   :1700,                  #[IN] Y- grid number of forward modeling
    'n1g'                   :350,                   #[IN] Z- grid number of forward modeling
    'd3g'                   :10.0,                  #[IN] X- spacing in meter
    'd2g'                   :10.0,                  #[IN] Y- spacing in meter
    'd1g'                   :10.0,                  #[IN] Z- spacing in meter
    'fevn_info'             :'list/eventinfo.csv',  #[IN] event list file
    'fdio'                  :'fdio',                #[IN] work directory for forward modeling 
    'Frequency filtering'   :'comment',
    'order'                 :6,                     #[IN] the order of the Butterworth lowpass filter
    'zero_phase'            :0,                     #[IN] filter trace in: 0, only forward direction; 1, both forward and reverse direction 
    "Misfit definition"     :'comment',
    'lnorm'                 :1,                     #[IN] misfit norm. 1: L2 norm
    'Gradient'              :'comment',
    'inc3'                  :1,                     #[IN] X- grid sampling for gradient
    'inc2'                  :1,                     #[IN] Y- grid sampling for gradient
    'inc1'                  :1,                     #[IN] Z- grid sampling for gradient
    'fcost'                 :'fdio/comm/misfit.dat',#[OUT] the misfit at point fmod
    'fgrad'                 :'fdio/comm/grad.bin',  #[OUT] the gradient at point fmod
    'fhess'                 :'fdio/comm/hess.bin',  #[OUT] the Hessian vector product at point [fmod, fvctr]
    'fvctr'                 :'fdio/comm/vctr.bin',  #[OUT] the vector pk in LCG algorithm
    'taper'                 :1,                     #[IN] gradient taper, 1: with, 0: without
    'ftaper'                :'taper/taper_all.bin', #[IN] gradient taper file
    'smooth_size'           :5,                     #[IN] median filter window for gradient direction, in odd number (5, 7, ...)
    'Optimization'          :'comment',
    'PNLCG'                 :'comment',
    'powell'                :1,                     #[IN] wheter to check Powell restart condition(1-yes, 0-no)
    'LBFGS'                 :'comment',
    'l'                     :5,                     #[IN] the maximum storing number of the pairs {sk,yk}(~3-7)
    'TRN'                   :'comment',
    'conv_CG'               :0,                     #[IN] initial status
    'niter_max_CG'          :12,                    #[IN] maximum number of inner conjugate gradient loop
    'flag_singularity'      :1,                     #[IN] singularity test, please check the code
    'flag_negative'         :2,                     #[IN] negative curvature test, please check the code
    'flag_truncation'       :1,                     #[IN] truncation test, please check the code
    'flag_min_res'          :1,                     #[IN] record minimum residual when applying the inner CG agrithm
    'pro_sing'              :1.0E-10,               #[IN] for singularity test
    'pro_curv'              :1.0E-10,               #[IN] for negative curvature test (if flag_negative == 1)
    'pro_trun'              :0.5,                   #[IN] for truncation test
    'eta'                   :0.9,                   #[IN] check the code
    'Steplength'            :'comment',
    'try_old_sl'            :1,                     #[IN] whether to use steplength from last iter.(1-yes, 0-no)
    'eps_scale'             :0.002,                 #[IN] max_m0*eps_scale as one unit to test the steplength
    'max_m0'                :4700.0,                #[IN] the maximum of initial model
    'alpha_lb'              :0.1,                   #[IN] lower limit of steplength
    'alpha_ub'              :10,                    #[IN] upper limit of steplength
    'mod_limits'            :1,                     #[IN] 1: check the bound of model value;0 no limit for model value
    'mod_lb'                :1000,                  #[IN] lower limit of model value in m/s
    'mod_ub'                :5000,                  #[IN] upper limit of model value in m/s
    'Inversion loop'        :'comment',
    'iter0'                 :1,                     #[IN] the initial iteration at current inversion stage/phase
    'niter_max'             :30,                    #[IN] maximum iterations, for future design
}
optim['prjroot'] = os.getcwd() # project path
optim['fdir'] = fdir
optim['fwav_list'] = 'list/proc_wav.csv'   #[OUT] source wavelet list
optim['fdata_list'] = 'list/proc_data.csv'  #[OUT] seismograms list
optim['fmod'] = optim['fmodfd'] #[OUT] current model [x] for inversion
optim['n3'] = len(np.arange(optim['n3g'])[::optim['inc3']])
optim['n2'] = len(np.arange(optim['n2g'])[::optim['inc2']])
optim['n1'] = len(np.arange(optim['n1g'])[::optim['inc1']])
optim['d1'] = optim['d1g'] * optim['inc1']
optim['d2'] = optim['d2g'] * optim['inc2']
optim['n'] = optim['n3g']*optim['n2g']*optim['n1g']       #only for TRN now, optimization works on FD grid after gradient interpolation
optim['np'] = optim['nprocx']*optim['nprocy'] # nodes for one simulation
# write work infomation
foptim = os.path.join(fdir,'optim.json')
fout = open(foptim,'w')
fout.write(json.dumps(optim,indent=4))
fout.close()
print("*"*60+"\n*Output file: "+foptim)

# data processing list
fevn_info = optim['fevn_info']
events = pd.read_csv(fevn_info)['station']
## wavelet
fin = [os.path.join('source', '%s.su'%evn) for evn in events]
fout = [os.path.join(optim['fdio'], evn, 'wavelet.su') for evn in events]
df = pd.DataFrame()
df['fin'] = fin
df['fout'] = fout
df.to_csv(optim['fwav_list'],index=False)

## observed data
fin = [os.path.join('data', '%s.su'%evn) for evn in events]
fout = [os.path.join(optim['fdio'], evn, 'obs_final.su') for evn in events]
df = pd.DataFrame()
df['fin'] = fin
df['fout'] = fout
df.to_csv(optim['fdata_list'],index=False)

os.makedirs(os.path.join(optim['fdio'], 'comm'), exist_ok=True)
for evn in events:
    os.makedirs(os.path.join(optim['fdio'], evn), exist_ok=True)

print("Finished...", __file__)
#####################