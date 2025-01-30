#!/usr/bin/env python3
"""
map vp in FD grid to slowness in tomo. grid 
"""
import os
import sys
import numpy as np
from scipy.ndimage import gaussian_filter

#Input paras
fdir        =sys.argv[1] #optim_LSQR
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
nx, ny      =optim['nxyfd']
incx, incy  =optim['incxy']
fmodfd      =optim['fmodfd']
fmod        =optim['fmod']
fmodfdref   =optim['fmodfdref']
fmodref     =optim['fmodref']

# m0 (current model): from FD grid to inverison grid 
modfd = np.fromfile(fmodfd, np.float32).reshape([nx, ny])
modinv = modfd[::incx, ::incy]
modinv = gaussian_filter(modinv, sigma=5/np.sqrt(8.0)) #only radius = 5 bins
modinv[:] = 1.0 / modinv[:]
modinv.tofile(fmod)

# mref (reference model): from FD grid to inverison grid
modfd = np.fromfile(fmodfdref, np.float32).reshape([nx, ny])
modinv = modfd[::incx, ::incy]
modinv = gaussian_filter(modinv, sigma=5/np.sqrt(8.0)) #only radius = 5 bins
modinv[:] = 1.0 / modinv[:]
modinv.tofile(fmodref)