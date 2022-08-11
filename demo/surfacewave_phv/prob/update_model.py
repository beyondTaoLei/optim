#!/usr/bin/env python3
"""
calculates descent direction
"""
import os
import sys
import numpy as np
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter

#Input paras
fdir        =sys.argv[1] #optim_LSQR
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
n3g         =optim['n3g']
n2g         =optim['n2g']
n1g         =optim['n1g']
n3          =optim['n3']
n2          =optim['n2']
n1          =optim['n1']
inc3        =optim['inc3']
inc2        =optim['inc2']
inc1        =optim['inc1']
d3g         =optim['d3g']
d2g         =optim['d2g']
d1g         =optim['d1g']
iterc       =optim['iterc']
eps_scale   =optim['eps_scale']
smooth_size =optim['smooth_size'] #in meter
mod_limits  =optim['mod_limits']
mod_lb      =optim['mod_lb']
mod_ub      =optim['mod_ub']
fmodo       =optim['fmodfd']

# interpretation
smooth_size =np.array([float(i) for i in smooth_size.split(':')])
sigmas=smooth_size/np.array([d3g, d2g, d1g])/np.sqrt(8.0)

# descent direction: from inverison grid to FD grid 
fdesc = os.path.join(fdir, 'iter'+str(iterc), 'desc.bin')
desc = np.fromfile(fdesc, np.float32)
grid3f, grid2f, grid1f = np.mgrid[0:(n3g-1):n3g*1j, 0:(n2g-1):n2g*1j, 0:(n1g-1):n1g*1j]
grid3i, grid2i, grid1i = np.mgrid[0:(n3-1):n3*1j, 0:(n2-1):n2*1j, 0:(n1-1):n1*1j]
grid3i, grid2i, grid1i = inc3*grid3i, inc2*grid2i, inc1*grid1i
points  = (grid3i.flatten(), grid2i.flatten(), grid1i.flatten())
pointsN = (grid3f.flatten(), grid2f.flatten(), grid1f.flatten())
descg = griddata(points, desc, pointsN, method='nearest')

# upate model
descg = descg.reshape([n3g, n2g, n1g])
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