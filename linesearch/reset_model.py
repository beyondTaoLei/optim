#!/usr/bin/env python3
"""
suppose that there are some steps in one iteration, if the inversion is 
running at current iteration (see, iterc) and stopped at any step 
by unknown factors, we would restart the update from iteration [iterc]. 
In this case, the current model should be reset from previous iteration.
"""
import os
import sys
import numpy as np

#Input paras
fdir        =sys.argv[1] #'optim'
foptim      =os.path.join(fdir,'optim.json')
paras       =eval(open(foptim).read())
iterc       =paras['iterc']
mod_limits  =paras['mod_limits']
mod_lb      =paras['mod_lb']
mod_ub      =paras['mod_ub']
fmodo       =paras['fmod']
fmodi       =os.path.join(fdir, 'iter'+str(iterc-1), 'mod.bin')
fdesc       =os.path.join(fdir, 'iter'+str(iterc-1), 'desc.bin')
fstep_estim =os.path.join(fdir, 'iter'+str(iterc-1),'step.json')
sl_paras = eval(open(fstep_estim).read())

#read m0 and delta
mod     =np.fromfile(fmodi, np.float32)
desc    =np.fromfile(fdesc, np.float32)

#make new test model
alpha=sl_paras['optim']*sl_paras['alpha_factor']
mod =mod+desc*alpha
#check the value limits
if mod_limits==1:
    print("check the model limits...")
    mod=np.clip(mod, mod_lb, mod_ub)
mod.astype(np.float32).tofile(fmodo)
print("you need to continue the update from iteration ", iterc)

#print("Finished...", __file__)
