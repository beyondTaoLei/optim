#!/usr/bin/env python3
"""
calculates descent direction
"""
import os
import sys
import csv
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from numpy import linalg as LA
from scipy.sparse.linalg import lsqr

#########################
#Input paras
fdir        =sys.argv[1] #optim_LSQR
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
fsensmat    =optim['fsensmat']
fwgtd       =optim['fwgtd']
fwgtm       =optim['fwgtm']
iterc       =optim['iterc']
niter_inner =optim['niter_inner_max']
lambda0     =optim['lambda0'] #suggest which is not more than 1.0E18
des_type    =optim['descent_type']
flag        =optim['modeltype']#1:smoothest, 2:flattest, 3:smallest
numreg      =optim['numreg']

# sensitivity matrix
idx = fsensmat.rfind('.')
suf = fsensmat[idx:]
if suf == '.npz':
    G = sparse.load_npz(fsensmat)
    N, M = G.shape
    # data residual: b = dobs - dsyn
    fdiff = os.path.join(fdir, 'iter'+str(iterc), 'diff.dat')
    b = np.loadtxt(fdiff, dtype=np.float32)
elif suf == '.h5':
    hf = h5py.File(fsensmat, 'r')
    G = hf['mat'][()]
    b = hf['b'][()]
    hf.close()
    G = sparse.csc_matrix(G)
    N, M = G.shape
elif suf == '.list':
    df = pd.read_csv(fsensmat)
    fsensmats   =list(df['list'])
    nums        =np.array(df['num'], np.int32)
    nlist = len(fsensmats)
    N = nums.sum()
    idxs = np.hstack([np.array([0], np.int32), nums.cumsum()])
    hf = h5py.File(fsensmats[0], 'r')
    N1, M = hf['mat'].shape
    hf.close()
    G = np.zeros([N, M], np.float32)
    b = np.zeros(N, np.float32)
    for i in range(nlist):
        idx1, idx2 = idxs[i], idxs[i+1]
        hf = h5py.File(fsensmats[i], 'r')
        G[idx1:idx2,:] = hf['mat'][()]
        b[idx1:idx2] = hf['b'][()]
        hf.close()
    G = sparse.csc_matrix(G)

print('finished matrix loading ...')

# data weighting
Wd = sparse.load_npz(fwgtd)

# model weighting
Wm = sparse.load_npz(fwgtm) # [nr, nr]

# current and reference model
fmod = os.path.join(fdir, 'iter'+str(iterc), 'mod.bin')
m0 = np.fromfile(fmod, np.float32)
if flag==1 or flag==2:
    mref = np.zeros_like(m0)
elif flag==3:
    fmodref = os.path.join(fdir, 'iter'+str(iterc), 'ref.bin')
    mref = np.fromfile(fmodref, np.float32)
else:
    raise ValueError("\nSo far we haven't implemented for multiple damping types\n")

# regularization term
yl, yr = [], []
for i in range(numreg):
    fnm = 'regularL'+str(i+1)+'.npz'
    fLm = os.path.join(fdir, 'regular', fnm)
    Lm = sparse.load_npz(fLm)
    print(Lm.shape, m0.shape)
    if des_type == 1:
        damphess = Lm
    else:
        damphess = Lm.dot(sparse.diags(m0))
    yl.append(Wm.dot(damphess))
    yr.append(Wm.dot(Lm.dot(mref - m0)))
# merge
ylhs = sparse.vstack(yl)
yrhs = np.hstack(yr)

# solve with lsqr
print("\ntest the regularization parameter lambda\n")
print("lambda0:",lambda0)
row_list=[["idx", "lambda", "phid", "phim", "bnorm"]]
xlhs = Wd.dot(G)
xrhs = Wd.dot(b)
Gext = sparse.vstack([xlhs, lambda0 * ylhs]).tocsc() #[N+M, M]
bext = np.hstack([xrhs, lambda0 * yrhs]) #[N+M]
Result = lsqr(Gext, bext, 0.0, 1e-08, 1e-08, 10000.0, niter_inner)
desc = Result[0]
phid = LA.norm(Gext[:N,:].dot(desc)-bext[:N])
phim = LA.norm(Gext[N:,:].dot(desc)-bext[N:]) / lambda0
bnorm = LA.norm(bext[:N]) #||Wd b||
row_list.append(list([0, format(lambda0,'12.4e'), format(phid,'12.4e'), 
                        format(phim,'12.4e'), format(bnorm,'12.4e')]))
print('%d %12.4e %12.4e %12.4e %12.4e'%(0, lambda0, phid, phim, bnorm))
fout=os.path.join(fdir, 'iter'+str(iterc), 'L_curve.dat')
with open(fout, 'w', newline='') as fp1:
    writer = csv.writer(fp1)
    writer.writerows(row_list)
fp1.close()
fdesc=os.path.join(fdir, 'iter'+str(iterc), 'desc.bin')
desc.astype(np.float32).tofile(fdesc)

#print("Finished...", __file__)