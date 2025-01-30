#!/usr/bin/env python3
"""
maps the input model grids into the inversion model grids, including 
input model and reference model, as
array [nx, ny, nz] ---> array [::incx, ::incy, ::incz] = [n3, n2, n1]
"""
import os
import sys
import numpy as np
from scipy.ndimage import gaussian_filter

#Input paras
fdir              =sys.argv[1] #optim_LSQR
foptim            =os.path.join(fdir,'optim.json')
optim             =eval(open(foptim).read())
nx, ny, nz        =optim['nxyzfd']
incx, incy, incz  =optim['incxyz']
fmodfd            =optim['fmodfd']
fmod              =optim['fmod']
fmodfdref         =optim['fmodfdref']
fmodref           =optim['fmodref']

sigma3 = 2.0/np.sqrt(8.0) if incx >1 else 0.0
sigma2 = 2.0/np.sqrt(8.0) if incy >1 else 0.0
sigma1 = 2.0/np.sqrt(8.0) if incz >1 else 0.0
sigmas = [sigma3, sigma2, sigma1]

# m0 (current model): from FD grid to inverison grid 
modfd = np.fromfile(fmodfd, np.float32).reshape([nx, ny, nz])
modinv = modfd[::incx, ::incy, ::incz]
modinv = gaussian_filter(modinv, sigma=sigmas, mode='nearest') #only radius = 2 bins
modinv.tofile(fmod)

# mref (reference model): from FD grid to inverison grid
modfd = np.fromfile(fmodfdref, np.float32).reshape([nx, ny, nz])
modinv = modfd[::incx, ::incy, ::incz]
modinv = gaussian_filter(modinv, sigma=sigmas, mode='nearest') #only radius = 2 bins
modinv.tofile(fmodref)