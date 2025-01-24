#!/usr/bin/env python3
import os
import sys
import json
import numpy as np
from angle_restart import angle_restart

#Input paras
fdir        =sys.argv[1] #'optim'
foptim      =os.path.join(fdir,'optim.json')
paras       =eval(open(foptim).read())
iterc       =paras['iterc']
iter0       =paras['iter0']
try_old_sl  =paras['try_old_sl']
eps_scale   =paras['eps_scale']
max_m0      =paras['max_m0']
alpha_lb    =paras['alpha_lb']
alpha_ub    =paras['alpha_ub']
fmisfit     =paras['fcost']
fmodo       =paras['fmod']
fmodi       =os.path.join(fdir, 'iter'+str(iterc), 'mod.bin')
fgrad       =os.path.join(fdir, 'iter'+str(iterc), 'grad.bin')
fdesc       =os.path.join(fdir, 'iter'+str(iterc), 'desc.bin')
fstep_estim =os.path.join(fdir, 'step.json')
#read m0 and delta
mod     =np.fromfile(fmodi, np.float32)
grad    =np.fromfile(fgrad, np.float32)
desc    =np.fromfile(fdesc, np.float32)
misfit0 =np.loadtxt(fmisfit).item()

flag =angle_restart(grad, desc)
if flag==1:
    desc.tofile(fdesc[:-4]+'_old.bin')
    desc =-grad
    desc.tofile(fdesc)
    # #reset iter0, special for l-BFGS
    # print("reset iter0 at %d"%iterc)
    # iter0 = iterc # reset step length
    # paras['iter0'] = iterc
    # fout =open(foptim,'w')
    # fout.write(json.dumps(paras,indent=4))
    # fout.close()
    
if try_old_sl==1:
    if iterc>iter0:
        print("use previous step length")
        sl_paras_old = eval(open(fstep_estim).read())
        alphap      =sl_paras_old['optim']
        #alpha_factor=sl_paras_old['alpha_factor']
    else:
        print("use initial step length")
        alphap=1.0
        #alpha_factor=eps_scale*max_m0/np.amax(np.fabs(desc))
else:
    print("use fixed step length")
    alphap=1.0

if eps_scale > 0.0:
    alpha_factor=eps_scale*max_m0/np.amax(np.fabs(desc))
else: # make sure the gradient has already scaled 
    alpha_factor=1.0

sl_paras={
    'max_m0'        :max_m0,    #useless
    'eps_scale'     :eps_scale,
    'beta1'         :2.0,
    'beta2'         :0.5,
    'beta3'         :0.05,
    'alpha_lb'      :alpha_lb,
    'alpha_ub'      :alpha_ub,
    'alphap'        :alphap,
    'alpha0'        :0.0,
    'alpha1'        :alphap,
    'alpha2'        :0.0,       #would be beta1*alpha1 or beta2*alpha1
    'misfit0'       :misfit0,
    'misfit1'       :0.0,       #will be updated later
    'misfit2'       :0.0,       #will be updated later
    'optim'         :alphap,    #will be updated finally
    'alpha_factor'  :alpha_factor
    }
#make new test model
alpha=sl_paras['alpha1']*sl_paras['alpha_factor']
mod =mod+desc*alpha
mod.astype(np.float32).tofile(fmodo)

fout = open(fstep_estim,'w')
fout.write(json.dumps(sl_paras,indent=4))
fout.close()

#print("Finished...", __file__)