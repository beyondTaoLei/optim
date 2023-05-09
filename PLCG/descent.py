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

def BasicHv(W, G, v):
    """
    return: vector G' W' W G v
    """
    temp = G.dot(v)
    temp = W.dot(temp)
    temp = W.transpose().dot(temp)
    a = G.transpose().dot(temp)
    
    return a

def func_phid(Wd, G, b, desc):
    """
    phid = ||Wd (Gx - b) ||, here b = dobs - dsyn.
    """
    temp = G.dot(desc)-b
    Gx_b = Wd.dot(temp)
    
    return LA.norm(Gx_b)

def func_phim(Wm, Lms, m0, desc, mref, des_type):
    """
    phim = ||Wm Lm ( m0 + x - mref)||
    """
    if des_type == 1:
        mm = m0 + desc - mref
    else:
        mm = m0 * (1.0 + desc) - mref
    
    out = 0.0
    for Lm in Lms:
        temp = Lm.dot(mm)
        mterm = Wm.dot(temp)
        out += LA.norm(mterm)
    
    return out

def diag_Hess(G, Wd, Wm, Lms, m0, des_type):
    """
    preconditional term: the diagonal of the Hessian matrix
    """
    N, M = G.shape
    Wd2 = Wd.diagonal()**2
    Wm2 = Wm.diagonal()**2
    diag1 = np.zeros(M, np.float32)
    diag2 = np.zeros(M, np.float32)
    inc=100 #memory safety
    for i in range(0, N, inc):
        W = Wd2[i:i+inc]
        H = G[i:i+inc, :]
        diag1[:] += sparse.csc_matrix(W).dot(H.power(2))#[1, inc] *[inc, M]
    
    for Lm in Lms:
        if des_type == 1:
            Lmm = Lm
        else:
            Lmm = Lm.dot(sparse.diags(m0))
        for i in range(0, M, inc):
            W = Wm2[i:i+inc]
            H = Lmm[i:i+inc, :]
            diag2[:] += sparse.csc_matrix(W).dot(H.power(2))
    
    return diag1, diag2
    
def HessV(G, Wd, Wm, Lms, mu, m0, des_type, v):
    """
    problem Gx = b, here b = dobs - dsyn.
    if searching for the absolute increment
        f = ||Wd (Gx - b) || + mu ||Wm Lm ( m0 + x - mref)||
        Hess = (G'Wd'WdG+ mu Lm'Wm'WmLm)
        rhs = G' Wd' Wd b + mu Lm' Wm' Wm Lm (mref-m0)
    if searching for the rate of increase
        f = ||Wd (Gx - b) || + mu ||Wm Lm ( [m0] (1 + x) - mref)||
        Hess = (G'Wd'WdG+ mu [m0]' Lm' Wm' Wm Lm [m0])
        rhs = G' Wd' Wd b + mu [m0]' Lm' Wm' Wm Lm (mref-m0)
    [m0] is diagonal matrix of m0
    return Hess*v
    """
    a1 = BasicHv(Wd, G, v)
    a2 = np.zeros_like(a1)
    if des_type == 1:
        for Lm in Lms:
            a2 += BasicHv(Wm, Lm, v)
        out = a1 + mu * a2
    else:
        mm = m0 * v
        for Lm in Lms:
            a2 += BasicHv(Wm, Lm, mm)
        out = a1 + mu * m0 * a2
    
    return out
    
def RHS(G, Wd, Wm, Lms, mu, m0, mref, des_type, b):
    """
    update the right hand side of equation
    return rhs, which is defined in HessV function
    """
    temp = Wd.dot(b)
    temp = Wd.transpose().dot(temp)
    b1 = G.transpose().dot(temp)
    b2 = np.zeros_like(b1)
    mm = mref-m0
    for Lm in Lms:
        b2 += BasicHv(Wm, Lm, mm)
    if des_type == 1:
        out = b1 + mu * b2
    else:
        out = b1 + mu * m0 * b2
    
    return out
    
def PLCG(G, Wd, Wm, Lms, mu, P, b, m0, mref, des_type, kmax=100, tol=1e-5):
    """
    preconditional CG
    P(G'Wd'WdG + muWm'Wm)x = P[G'Wd'Wd*(dobs-dsyn)+muWm'Wm(mref-mc)]
    G:  [N x M] sensitivity matrix
    Wd: [N x N] data weighting
    Wm: [M x M] model weighting (regularization term)
    P:  [M] preconditional vector
    b:  [N] misfit vector (dobs - dsyn)
    m0: [M] initial model
    mref: [M] reference model
    """
    #update the right hand side of equation
    B = RHS(G, Wd, Wm, Lms, mu, m0, mref, des_type, b)
    #initialization
    #xk = m0.copy()
    #rk = HessV(G, Wd, Wm, Lms, mu, m0, des_type, xk) - B
    xk = np.zeros_like(m0)
    rk = -B
    yk = rk / P
    pk = -yk
    ratio = LA.norm(rk)/LA.norm(B)
    
    #iteration loop
    k = 0
    #print(k, xk, ratio)
    while k < kmax and ratio > tol:
        apk = HessV(G, Wd, Wm, Lms, mu, m0, des_type, pk)
        rkyk = np.inner(rk, yk)
        alpha = rkyk / np.inner(pk, apk)
        
        xk += alpha * pk
        rk += alpha * apk
        
        yk = rk / P
        beta = np.inner(rk, yk) / rkyk
        pk[:] = -yk + beta * pk
        ratio = LA.norm(rk)/LA.norm(B)
        
        k += 1
        #print(k, xk, ratio)
    
    return xk

#########################
#Input paras
fdir        =sys.argv[1] #optim_PLCG
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
percond     =optim['percond']
fsensmat    =optim['fsensmat']
fwgtd       =optim['fwgtd']
fwgtm       =optim['fwgtm']
iterc       =optim['iterc']
niter_inner =optim['niter_inner_max']
lambda0     =optim['lambda0'] #suggest which is not more than 1.0E18
des_type    =optim['descent_type'] #1: absolute 0: rate of increase
flag        =optim['modeltype']#1:smoothest, 2:flattest, 3:smallest
numreg      =optim['numreg']
tol         =0.01 #optim['tolerance'] tolerance for ||Ax-B||/||B||

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

#data residual: b = dobs - dsyn
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
Lms = []
for i in range(numreg):
    fnm = 'regularL'+str(i+1)+'.npz'
    fLm = os.path.join(fdir, 'regular', fnm)
    Lms.append(sparse.load_npz(fLm))

# preconditional opt
P = np.ones(M)
diag1, diag2 = diag_Hess(G, Wd, Wm, Lms, m0, des_type)

# solve with plcg
print("\ntest the regularization parameter lambda\n")
print("lambda0:",lambda0)
row_list=[["idx", "lambda", "phid", "phim", "bnorm"]]

mu = lambda0 * lambda0
# preconditional term
if percond == 1: #update P
    P[:] = diag1 + mu * diag2
    #P2 = np.sort(P[P>0.0], kind='quicksort') #should be careful for quicksort method
    #P = P + P2[int(M*0.03)]
    
desc = PLCG(G, Wd, Wm, Lms, mu, P, b, m0, mref, des_type, niter_inner, tol)
phid = func_phid(Wd, G, b, desc)
phim = func_phim(Wm, Lms, m0, desc, mref, des_type)
bnorm = LA.norm(b)
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
    
    