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
nxg         =optim['n2g']
nyg         =optim['n1g']
incx        =optim['inc2']
incy        =optim['inc1']
fmodfd      =optim['fmodfd']
fmod        =optim['fmod']
fmodfdref   =optim['fmodfdref']
fmodref     =optim['fmodref']

# m0 (current model): from FD grid to inverison grid 
modfd = np.fromfile(fmodfd, np.float32).reshape([nxg, nyg])
modinv = modfd[::incx, ::incy]
modinv = gaussian_filter(modinv, sigma=5/np.sqrt(8.0)) #only radius = 5 bins
modinv[:] = 1.0 / modinv[:]
modinv.tofile(fmod)

# mref (reference model): from FD grid to inverison grid
modfd = np.fromfile(fmodfdref, np.float32).reshape([nxg, nyg])
modinv = modfd[::incx, ::incy]
modinv = gaussian_filter(modinv, sigma=5/np.sqrt(8.0)) #only radius = 5 bins
modinv[:] = 1.0 / modinv[:]
modinv.tofile(fmodref)