#!/usr/bin/env python3
"""
Usage:
    python disp_processing/S0_prepare.py
"""
import os
import sys
import json
import numpy as np
import pandas as pd

#Input paras
fdir    =sys.argv[1] # optim
# make directory
os.makedirs(fdir, exist_ok=True)
os.makedirs('log', exist_ok=True)
os.makedirs('job', exist_ok=True)
os.makedirs('list', exist_ok=True)
os.makedirs('fdio', exist_ok=True)
os.makedirs('data', exist_ok=True)
os.makedirs('model', exist_ok=True)
os.makedirs('source', exist_ok=True)
os.makedirs(os.path.join('model', 'inv'), exist_ok=True)
os.makedirs(os.path.join('model', 'init'), exist_ok=True)
os.makedirs(os.path.join('model', 'true'), exist_ok=True)
os.makedirs(os.path.join('fdio', 'comm'), exist_ok=True)

optim ={
    'Program path'              :'comment',
    'optimroot'                 :'/data/ess-leit/software/optim',             #[IN] path to optim package
    'mpiexec'                   :'/share/intel/2020u4/compilers_and_libraries_2020.4.304/linux/mpi/intel64/bin/mpiexec', #[IN] path to C MPI parallel launcher
    'mpiexep'                   :'mpirun',                #[IN] path to Python MPI parallel launcher
    'Forward part'              :'comment',
    'ncores'                    :40,                      #[IN] number of CPU cores in one node
    'nproc'                     :[4, 1],                  #[IN] MPI processors in X- and Z-dir
    'fevn_info'                 :'list/eventinfo.csv',    #[IN] event list file
    'fmodfd'                    :'model/inv/inv.vp',      #[IN] current model in work memory
    'nxzfd'                     :[1700, 350],             #[IN] FD grid number in X- and Z-dir
    'oxzfd'                     :[0.0, 0.0],              #[IN] FD grid origin coordinates
    'dxzfd'                     :[10.0, 10.0],            #[IN] FD grid spacing in X- and Z-dir
    'nt'                        :6500,                    #[IN] time steps
    'delta'                     :1.0e-3,                  #[IN] time spacing in second
    'tshift_adj'                :0.2,                     #[IN] deadline for backward modeling, in second
    'freesurface'               :1,                       #[IN] 0: without; 1: with the free surface
    'pml'                       :15,                      #[IN] PML layers
    'pmlvel'                    :3000.0,                  #[IN] PML velocity in m/s
    'comp'                      :['P'],                   #[IN] measurement components, ['P', 'Vx', 'Vy']
    'Data processing'           :'comment',                 
    'velocity'                  :1,                       #[IN] 0: apply integral, 1: not apply integral for the measurement components
    'dirmute'                   :0,                       #[IN] direct wave muting, 0: no muting, mute (1) by speed, (2) by f-k, (3) by file: mute.h5
    'dirmute_taperlen'          :0.4,                     #[IN] taper length if dirmute=1
    'dirmute_t0'                :0.20,                    #[IN] time offset caused by the delayed excitation of the source (s) if dirmute=1
    'dirmute_lower'             :0.0,                     #[IN] propagating speed of direct wave if dirmute=1; or stop band low corner speed if dirmute=2
    'dirmute_upper'             :1700.0,                  #[IN] stop band high corner speed if dirmute=2  
    'timefilt'                  :1,                       #[IN] whether to apply lowpass filtering, 1: yes, 0: no
    'fc'                        :5.0,                     #[IN] low cut-off frequency in Hz
    'order'                     :6,                       #[IN] filter corners/order, 4 or 6
    'zero_phase'                :0,                       #[IN] zero phase shift when filtering, 1: True, 0: False 
    'trkill'                    :0,                       #[IN] whether to apply trace killing, 1: yes, 0: no
    'trkill_lower'              :20.0,                    #[IN] minimum offset if trkill=1, kill the traces outside [trkill_lower, trkill_upper] 
    'trkill_upper'              :10000.0,                 #[IN] maximum offset if trkill=1
    'timewin'                   :1,                       #[IN] whether to apply time windowing, 1: yes, 0: no
    'timewin_lower'             :0.1,                     #[IN] starting time if timewin=1
    'timewin_upper'             :10.0,                    #[IN] end time if timewin=1
    'Gradient part'             :'comment',                   
    'Lnorm'                     :2,                       #[IN] misfit function, 2: L2 norm, 3: cross-correlation
    'epsilon'                   :0.03,                    #[IN] energy preconditioning, 0: without preconditioning, >0: water level
    'grad_taper_shot'           :1,                       #[IN] whether to apply gradient taper for each shot, 1: yes, 0: no
    'srtradius1'                :1,                       #[IN] starting point of taper far away the source or receiver if grad_taper_shot=1
    'srtradius2'                :10,                      #[IN] end point of taper far away the source or receiver if grad_taper_shot=1
    'grad_filt'                 :5,                       #[IN] >0: gradient smoothing size, 0: no smoothing
    'grad_taper'                :1,                       #[IN] whether to apply taper for final gradient, 1: yes, 0: no
    'grad_taper_file'           :'list/taper.bin',        #[IN] the taper file if grad_taper=1
    'Optimization'              :'comment',                   
    'LBFGS'                     :'comment',                     
    'l'                         :5,                       #[IN] the maximum storing number of the pairs {sk,yk}(~3-7)
    'Steplength'                :'comment',                   
    'try_old_sl'                :1,                       #[IN] whether to use steplength from last iter.(1-yes, 0-no)
    'max_m0'                    :4700.0,                  #[IN] the maximum of initial model
    'eps_scale'                 :0.002,                   #[IN] max_m0*eps_scale as one unit to test the steplength
    'alpha_lb'                  :0.1,                     #[IN] lower limit of steplength
    'alpha_ub'                  :10.0,                    #[IN] upper limit of steplength
    'mod_limits'                :1,                       #[IN] whether to check the bound of model value, 1: yes, 0: no
    'mod_lb'                    :1000.0,                  #[IN] lower limit of model value if mod_limits=1
    'mod_ub'                    :5000.0,                  #[IN] upper limit of model value if mod_limits=1
    'Inversion loop'            :'comment',                     
    'abort'                     :0,                       #[IN] intial stage phase
    'iter0'                     :31,                      #[IN] the initial iteration at current inversion stage/phase
    'niter_min'                 :30,                      #[IN] minimum iterations
    'niter_max'                 :60,                      #[IN] maximum iterations, for future design
    'pro'                       :0.001,                   #[IN] (Ei2 - Ei) / Ei2
    'pro2'                      :0.001,                   #[IN] k2 / k1 
    'pro3'                      :0.001,                   #[IN] average2 / average1 
    }

#intermediate files 
optim['fwiroot']    =os.path.dirname(os.path.dirname(__file__))
optim['fconf']      =os.path.join(optim['fwiroot'], 'conf')
optim['fdroot']     =os.path.join(optim['fwiroot'], 'SOFI2DAC') # FD software path
optim['fdir']       =fdir
optim['prjroot']    =os.getcwd() #'/home/tao/works/software/tomo3d/prj'
optim['prjname']    =optim['prjroot'].split('/')[-1]
optim['fdio']       ='fdio'
optim['fdata']      ='data'
optim['fmodel']     ='model'
optim['fsource']    ='source'
optim['fmod']       =optim['fmodfd']                    # grid of FD = grid of inversion
optim['fcost']      ='list/fcost.dat'                   #[OUT] misfit function
optim['fgrad']      ='list/grad.bin'                    #[OUT] the gradient at point fmod
optim['fmisfit']    ='list/misfit.csv'                  #[OUT] the output from function S3_adjoint_src.py

# optim package
optim['n2'] = optim['nxzfd'][0]
optim['n1'] = optim['nxzfd'][1]
optim['d2'] = optim['dxzfd'][0]
optim['d1'] = optim['dxzfd'][1]
optim['n'] = optim['n2']*optim['n1'] #only for TRN now

foptim = os.path.join(fdir,'optim.json')
fp = open(foptim,'w')
fp.write(json.dumps(optim,indent=4))
fp.close()

print("Finished...", __file__)