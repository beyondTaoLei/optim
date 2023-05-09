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
fcmp        =optim['fcmp']
fgrad       =optim['fgrad']
nper        =optim['nper']
ntraces     =optim['ntraces']
n3          =optim['n3g']
n2          =optim['n2g']
n1          =optim['n1g']
d3          =optim['d3g']
d2          =optim['d2g']
d1          =optim['d1g']
smooth_size =optim['smooth_size'] #in meter
percond     =optim['percond']
eps         =optim['EPSILON']
maxmode     =optim['maxmode']
modes       =maxmode + 1

# interpretation
smooth_size = np.array([float(i) for i in smooth_size.split(':')])
sigmas=smooth_size/np.array([d3, d2, d1])/np.sqrt(8)
df = pd.read_csv(fperiods)
periods = np.array(df['T0'])

# read indexes of CMP
df = pd.read_csv(fcmp)
offsets = df['offset']
    
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
grad = np.zeros([ntraces, n1], np.float32)
if percond == 1:
    diag1 = np.zeros([ntraces, n1], np.float32)
for mode in range(modes):
    # read kernels
    for iper in range(nper):
        fname='ker_M%1d_T%.2f.bin'%(mode, periods[iper])
        fin = os.path.join(wdir, fname)
        dcdm = np.fromfile(fin, np.float32).reshape([ntraces, n1])
        grad[:]+=np.einsum('i,ij->ij',-b[mode,:,iper], dcdm)
        if percond == 1:
            diag1[:]+= dcdm**2


######### preprocessing #########
# eliminate the outliers
where_are_NaNs=np.isnan(grad)
grad[where_are_NaNs] = 0.0 #unknown error

# preconditioning
if percond == 1:
    mean1 = diag1.mean(axis=0)
    grad[:] = grad/(diag1 +mean1.max()*eps)

# interpolation (ntraces --->(n3, n2) )
gradfinal = np.zeros([n3, n2, n1], np.float32)
grid_3, grid_2 = np.mgrid[0:(n3-1):n3*1j, 0:(n2-1):n2*1j]
points=np.array([offsets//n2, offsets%n2]).T
for i in range(n1):
    gradfinal[:, :, i] = griddata(points, grad[:, i], (grid_3, grid_2), method='nearest')

# smooth
gradfinal[:] = gaussian_filter(gradfinal, sigmas, mode='nearest')#sigma
##################################
gradfinal.tofile(fgrad)

print("Finished...", __file__)