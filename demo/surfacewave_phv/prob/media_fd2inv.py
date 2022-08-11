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
n3g         =optim['n3g']
n2g         =optim['n2g']
n1g         =optim['n1g']
inc3        =optim['inc3']
inc2        =optim['inc2']
inc1        =optim['inc1']
fmodfd      =optim['fmodfd']
fmod        =optim['fmod']
fmodfdref   =optim['fmodfdref']
fmodref     =optim['fmodref']

sigma3 = 2.0/np.sqrt(8.0) if inc3 >1 else 0.0
sigma2 = 2.0/np.sqrt(8.0) if inc2 >1 else 0.0
sigma1 = 2.0/np.sqrt(8.0) if inc1 >1 else 0.0
sigmas = [sigma3, sigma2, sigma1]

# m0 (current model): from FD grid to inverison grid 
modfd = np.fromfile(fmodfd, np.float32).reshape([n3g, n2g, n1g])
modinv = modfd[::inc3, ::inc2, ::inc1]
modinv = gaussian_filter(modinv, sigma=sigmas, mode='nearest') #only radius = 2 bins
modinv.tofile(fmod)

# mref (reference model): from FD grid to inverison grid
modfd = np.fromfile(fmodfdref, np.float32).reshape([n3g, n2g, n1g])
modinv = modfd[::inc3, ::inc2, ::inc1]
modinv = gaussian_filter(modinv, sigma=sigmas, mode='nearest') #only radius = 2 bins
modinv.tofile(fmodref)