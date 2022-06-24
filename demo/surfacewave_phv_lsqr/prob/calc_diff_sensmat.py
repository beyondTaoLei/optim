#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd
from scipy import sparse
from netCDF4 import Dataset
from numpy import linalg as LA

def fnm_comm(mode, per):
    return 'M%1d_T%.2f'%(mode, per)

#Input paras
fdir        =sys.argv[1]
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
wdir        =optim['wdir']
odir        =optim['odir']
nx          =optim['n2g']
nz          =optim['n1g']
dx          =optim['d2g']
fcost       =optim['fcost']
fdiff       =optim['fdiff']
fsensmat    =optim['fsensmat']
maxmode     =optim['maxmode']
modes       =maxmode + 1

# read phase velocity
fnm = os.path.join(wdir, 'vmap_M%1d.nc'%0)
ncin = Dataset(fnm)
nper = ncin.dimensions['K'].size
period = ncin.variables['period'][:,:].data
ncin.close()
periods = period[0, :]
    
b = np.zeros([modes, nx, nper], np.float32)
for mode in range(modes):
    # phase velocity
    fnm = os.path.join(wdir, 'vmap_M%1d.nc'%mode)
    ncin = Dataset(fnm)
    phvsyn = ncin.variables['vel'][:,:].data
    ncin.close()
    fnm = os.path.join(odir, 'vmap_M%1d.nc'%mode)
    ncin = Dataset(fnm)
    phvobs = ncin.variables['vel'][:,:].data
    ncin.close()
    diff = np.where((phvobs>0.0) &(phvsyn>0.0), phvobs-phvsyn, -1000.0)
    b[mode, :, :] = diff

# generate the sensitivity matrix
sensmat = np.zeros([modes, nper, nx, nz], np.float32)
for mode in range(modes):
    # read kernels
    for iper in range(nper):
        fname='ker_'+fnm_comm(mode, periods[iper])+'.nc'
        fin = os.path.join(wdir, fname)
        ncin = Dataset(fin)
        dcdm = ncin.variables['vs'][:,:].data
        ncin.close()
        sensmat[mode, iper, :, :] = dcdm

# save the difference and sensmat
idxs = np.where(b>-900.0) #[modes, nx, nper]
b = b[idxs] # vector 
res = []
rows = np.zeros(nz, np.int32)
offset0 = np.arange(nz)
num = nx*nz
for mode, ix, iper in zip(idxs[0], idxs[1], idxs[2]):
    cols = ix*nz + offset0
    profile = sensmat[mode, iper, ix, :]
    cscmat = sparse.csr_matrix((profile, (rows, cols)), shape=(1, num))
    res.append(cscmat)
    
# save sensitivity matrix with format of csc_matrix
sen_mat = sparse.vstack(res).tocsc()
sparse.save_npz(fsensmat, sen_mat)
np.savetxt(fdiff, b, fmt='%.4e')

misfit=0.5*LA.norm(b)**2
fp = open(fcost,'w')
fp.write('%.6e\n'% misfit)
fp.close()

print("Finished...", __file__)