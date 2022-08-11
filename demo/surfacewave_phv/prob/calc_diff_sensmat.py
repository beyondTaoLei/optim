#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd
from scipy import sparse
from numpy import linalg as LA
from scipy.ndimage import gaussian_filter

#Input paras
fdir        =sys.argv[1]
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
wdir        =optim['wdir']
odir        =optim['odir']
fperiods    =optim['fperiods']
fcmp        =optim['fcmp']
fcost       =optim['fcost']
fdiff       =optim['fdiff']
fsensmat    =optim['fsensmat']
nper        =optim['nper']
ntraces     =optim['ntraces']
n3g         =optim['n3g']
n2g         =optim['n2g']
n1g         =optim['n1g']
inc3        =optim['inc3']
inc2        =optim['inc2']
inc1        =optim['inc1']
n3          =optim['n3']
n2          =optim['n2']
n1          =optim['n1']
maxmode     =optim['maxmode']
modes       =maxmode + 1

# interpretation
df = pd.read_csv(fperiods)
periods = np.array(df['T0'])
sigmas = np.array([2.0, 2.0, 2.0])/np.sqrt(8.0)

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
sensmat = np.zeros([modes, nper, ntraces, n1g], np.float32)
for mode in range(modes):
    # read kernels
    for iper in range(nper):
        fname='ker_M%1d_T%.2f.bin'%(mode, periods[iper])
        fin = os.path.join(wdir, fname)
        dcdm = np.fromfile(fin, np.float32).reshape([ntraces, n1g])
        sensmat[mode, iper, :, :] = dcdm

# save the difference and sensmat
idxs = np.where(b>-900.0) #[modes, n2, nper]
b = b[idxs] # vector
res = []
rows = np.zeros(n1, np.int32)
offset0 = np.arange(n1)
num = n3*n2*n1
n1g2 = n1*inc1
profile = np.zeros(n1g2, np.float32)
for mode, itr, iper in zip(idxs[0], idxs[1], idxs[2]):
    i3, i2 = offsets[itr]//n2g, offsets[itr]%n2g
    I3, I2 = i3//inc3, i2//inc2
    cols = (I2 + I3*n2)*n1 + offset0
    profile[:n1g] = sensmat[mode, iper, itr, :]
    profile2 = profile.reshape(n1, inc1).sum(axis=1)
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