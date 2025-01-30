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
nx, ny      =optim['nxyfd']
dx, dy      =optim['dxyfd']
incx, incy  =optim['incxy']
nxi         =optim['n2']
nyi         =optim['n1']
iterc       =optim['iterc']
eps_scale   =optim['eps_scale']
smooth_size =optim['smooth_size'] #in meter
mod_limits  =optim['mod_limits']
mod_lb      =optim['mod_lb']
mod_ub      =optim['mod_ub']
fmodo       =optim['fmodfd']

sigmas=smooth_size/np.array([dx, dy])/np.sqrt(8)

# descent direction: from inverison grid to FD grid 
grid_xg = np.tile(np.arange(nx),(ny,1)).T
grid_zg = np.tile(np.arange(ny),(nx,1))
grid_x = np.tile(np.arange(0,nx,incx),(nyi,1)).T
grid_z = np.tile(np.arange(0,ny,incy),(nxi,1))
points = np.array([grid_x.flatten(), grid_z.flatten()]).T
fdesc = os.path.join(fdir, 'iter'+str(iterc), 'desc.bin')
desc = np.fromfile(fdesc, np.float32)
descg = griddata(points, desc, (grid_xg, grid_zg), method='nearest') #slowness increment

# upate model
descg = gaussian_filter(descg, sigmas, mode='nearest')#sigma
descg = descg.flatten()
max_desc = np.amax(np.fabs(descg))
#print(eps_scale*max_m0/max_desc)
fmodi =os.path.join(fdir, 'iter'+str(iterc), 'modg.bin') 
mod = np.fromfile(fmodi, np.float32) #vp
# s = s0 + delta s
slowness = 1.0 / mod
meanv = slowness[slowness>0.0].mean()
slowness = slowness + descg / max_desc * eps_scale * meanv
mod = 1.0 / slowness
# check the value limits
if mod_limits==1:
    print("check the model limits...")
    print(mod_lb, mod_ub)
    mod=np.clip(mod, mod_lb, mod_ub)
mod.astype(np.float32).tofile(fmodo)

print('Finished...', __file__)