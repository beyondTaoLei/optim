#!/usr/bin/env python3
"""
calculate phase velocity errors and map the sensitivity kernel 
of each measurement into the inversion grids, as
1D array [nz] ---> 1D array [n1] ---> array([n3, n2, n1]).flatten()


"""
import os
import sys
import numpy as np
import pandas as pd
from scipy import sparse
from numpy import linalg as LA
from scipy.ndimage import gaussian_filter

#Input paras
fdir              =sys.argv[1]
foptim            =os.path.join(fdir,'optim.json')
optim             =eval(open(foptim).read())
wdir              =optim['wdir']
odir              =optim['odir']
fperiods          =optim['fperiods']
fstn              =optim['fstn']
fcost             =optim['fcost']
fdiff             =optim['fdiff']
fsensmat          =optim['fsensmat']
nx, ny, nz        =optim['nxyzfd']
incx, incy, incz  =optim['incxyz']
n3                =optim['n3']
n2                =optim['n2']
n1                =optim['n1']
maxmode           =optim['maxmode']
modes             =maxmode + 1

# interpretation
df = pd.read_csv(fperiods)
nper = df.shape[0]
periods = np.array(df['T0'])
sigmas = np.array([2.0, 2.0, 2.0])/np.sqrt(8.0)

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
sensmat = np.zeros([modes, nper, ntraces, nz], np.float32)
for mode in range(modes):
    # read kernels
    for iper in range(nper):
        fname='ker_M%1d_T%.2f.bin'%(mode, periods[iper])
        fin = os.path.join(wdir, fname)
        dcdm = np.fromfile(fin, np.float32).reshape([ntraces, nz])
        sensmat[mode, iper, :, :] = dcdm

# save the difference and sensmat
idxs = np.where(b>-900.0) #[modes, n2, nper]
b = b[idxs] # vector
res = []
rows = np.zeros(n1, np.int32)
offset0 = np.arange(n1)
num = n3*n2*n1
nz2 = n1*incz
profile = np.zeros(nz2, np.float32)
for mode, itr, iper in zip(idxs[0], idxs[1], idxs[2]):
    # i3, i2 = offsets[itr]//ny, offsets[itr]%ny
    # I3, I2 = i3//incx, i2//incy
    I3, I2 = ixs[itr]//incx, iys[itr]//incy
    cols = (I2 + I3*n2)*n1 + offset0
    profile[:nz] = sensmat[mode, iper, itr, :]
    profile2 = profile.reshape(n1, incz).sum(axis=1)
    cscmat = sparse.csr_matrix((profile2, (rows, cols)), shape=(1, num))
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