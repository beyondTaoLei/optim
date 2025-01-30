#!/usr/bin/env python3
"""
update model with steps:
1. interpolate descent direction from inv. grids into input model grids;
2. smooth descent direction;
3. update model in input model grid.
"""
import os
import sys
import numpy as np
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter

#Input paras
fdir              =sys.argv[1] #optim_LSQR
foptim            =os.path.join(fdir,'optim.json')
optim             =eval(open(foptim).read())
nx, ny, nz        =optim['nxyzfd']
dx, dy, dz        =optim['dxyzfd']
incx, incy, incz  =optim['incxyz']
n3                =optim['n3']
n2                =optim['n2']
n1                =optim['n1']
iterc             =optim['iterc']
eps_scale         =optim['eps_scale']
smooth_size       =optim['smooth_size'] #in meter
mod_limits        =optim['mod_limits']
mod_lb            =optim['mod_lb']
mod_ub            =optim['mod_ub']
fmodo             =optim['fmodfd']

# interpretation
sigmas=smooth_size/np.array([dx, dy, dz])/np.sqrt(8.0)

# descent direction: from inverison grid to FD grid 
fdesc = os.path.join(fdir, 'iter'+str(iterc), 'desc.bin')
desc = np.fromfile(fdesc, np.float32)
grid3f, grid2f, grid1f = np.mgrid[0:(nx-1):nx*1j, 0:(ny-1):ny*1j, 0:(nz-1):nz*1j]
grid3i, grid2i, grid1i = np.mgrid[0:(n3-1):n3*1j, 0:(n2-1):n2*1j, 0:(n1-1):n1*1j]
grid3i, grid2i, grid1i = incx*grid3i, incy*grid2i, incz*grid1i
points  = (grid3i.flatten(), grid2i.flatten(), grid1i.flatten())
pointsN = (grid3f.flatten(), grid2f.flatten(), grid1f.flatten())
descg = griddata(points, desc, pointsN, method='nearest')

# upate model
descg = descg.reshape([nx, ny, nz])
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