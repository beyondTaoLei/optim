#!/usr/bin/env python3
"""
calculates descent direction
"""
import os
import sys
import csv
import numpy as np
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
G = sparse.load_npz(fsensmat)
N, M = G.shape

# data weighting
if os.path.exists(fwgtd):
    Wd = sparse.load_npz(fwgtd)
else:
    Wd = sparse.eye(N)

# model weighting
if os.path.exists(fwgtm):
    Wm = sparse.load_npz(fwgtm)
else:
    Wm = sparse.eye(M)

# data residual: b = dobs - dsyn
fdiff = os.path.join(fdir, 'iter'+str(iterc), 'diff.dat')
b = np.loadtxt(fdiff, dtype=np.float32)

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
Gext = sparse.vstack([xlhs, lambda0 * ylhs]) #[N+M, M]
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