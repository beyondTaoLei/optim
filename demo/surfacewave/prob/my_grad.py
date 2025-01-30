#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd
from scipy import sparse
from numpy import linalg as LA
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter

#Input paras
fdir        =sys.argv[1]
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
wdir        =optim['wdir']
odir        =optim['odir']
fperiods    =optim['fperiods']
fstn        =optim['fstn']
fgrad       =optim['fgrad']
nx, ny, nz  =optim['nxyzfd']
dx, dy, dz  =optim['dxyzfd']
smooth_size =optim['smooth_size'] #in meter
precond     =optim['precond']
eps         =optim['EPSILON']
maxmode     =optim['maxmode']
modes       =maxmode + 1

# interpretation
sigmas=smooth_size/np.array([dx, dy, dz])/np.sqrt(8)
df = pd.read_csv(fperiods)
nper = df.shape[0]
periods = np.array(df['T0'])

# read indexes of stations
df = pd.read_csv(fstn)
ntraces = df.shape[0]
# offsets = df['offset']
ixs = df['ix']
iys = df['iy']
    
b = np.zeros([modes, ntraces, nper], np.float32)
for mode in range(modes):
    # phase velocity
    fnm = os.path.join(wdir, 'vmap_M%1d.bin'%mode)
    phvsyn = np.fromfile(fnm, np.float32).reshape([ntraces, nper])
    fnm = os.path.join(odir, 'vmap_M%1d.bin'%mode)
    phvobs = np.fromfile(fnm, np.float32).reshape([ntraces, nper])
    diff = np.where((phvobs>0.0) &(phvsyn>0.0), phvobs-phvsyn, -1000.0)
    b[mode, :, :] = diff

# generate the sensitivity matrix
b[b<-900.0] = 0.0
grad = np.zeros([ntraces, nz], np.float32)
if precond == 1:
    diag1 = np.zeros([ntraces, nz], np.float32)
for mode in range(modes):
    # read kernels
    for iper in range(nper):
        fname='ker_M%1d_T%.2f.bin'%(mode, periods[iper])
        fin = os.path.join(wdir, fname)
        dcdm = np.fromfile(fin, np.float32).reshape([ntraces, nz])
        grad[:]+=np.einsum('i,ij->ij',-b[mode,:,iper], dcdm)
        if precond == 1:
            diag1[:]+= dcdm**2


######### preprocessing #########
# eliminate the outliers
where_are_NaNs=np.isnan(grad)
grad[where_are_NaNs] = 0.0 #unknown error

# preconditioning
if precond == 1:
    mean1 = diag1.mean(axis=0)
    grad[:] = grad/(diag1 +mean1.max()*eps)

# interpolation (ntraces --->(nx, ny) )
gradfinal = np.zeros([nx, ny, nz], np.float32)
grid_3, grid_2 = np.mgrid[0:(nx-1):nx*1j, 0:(ny-1):ny*1j]
# points=np.array([offsets//ny, offsets%ny]).T
points=np.array([ixs, iys]).T
for i in range(nz):
    gradfinal[:, :, i] = griddata(points, grad[:, i], (grid_3, grid_2), method='nearest')

# smooth
gradfinal[:] = gaussian_filter(gradfinal, sigmas, mode='nearest')#sigma
##################################
gradfinal.tofile(fgrad)

print("Finished...", __file__)