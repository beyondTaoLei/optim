#!/usr/bin/env python3
"""
calculates descent direction
"""
import os
import sys
import numpy as np
from scipy.ndimage import gaussian_filter

#Input paras
fdir        =sys.argv[1] #optim_LSQR
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
n2g         =optim['n2g']
n1g         =optim['n1g']
d2g         =optim['d2g']
d1g         =optim['d1g']
iterc       =optim['iterc']
eps_scale   =optim['eps_scale']
smooth_size =optim['smooth_size'] #in meter
mod_limits  =optim['mod_limits']
mod_lb      =optim['mod_lb']
mod_ub      =optim['mod_ub']
fmodo       =optim['fmodfd']

smooth_size =np.array([float(i) for i in smooth_size.split(':')])
sigmas=smooth_size/np.array([d2g, d1g])/np.sqrt(8)

# descent direction
fdesc = os.path.join(fdir, 'iter'+str(iterc), 'desc.bin')
descg = np.fromfile(fdesc, np.float32).reshape([n2g, n1g])

# upate model
descg = gaussian_filter(descg, sigmas, mode='nearest')#sigma
descg = descg.flatten()
max_desc = np.amax(np.fabs(descg))
#print(eps_scale*max_m0/max_desc)
fmodi =os.path.join(fdir, 'iter'+str(iterc), 'modg.bin') 
mod = np.fromfile(fmodi, np.float32)
mod = mod + descg / max_desc * eps_scale * mod.max()
# check the value limits
if mod_limits==1:
    print("check the model limits...")
    print(mod_lb, mod_ub)
    mod=np.clip(mod, mod_lb, mod_ub)
mod.astype(np.float32).tofile(fmodo)

print('Finished...', __file__)