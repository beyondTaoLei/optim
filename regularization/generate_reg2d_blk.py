#!/usr/bin/env python3
import os
import sys
import json
import numpy as np
from scipy import sparse
from util import opt_identity, opt_deriv, opt_laplace

#Input paras
fdir        =sys.argv[1] #optim
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
n1          =optim['n1']
n2          =optim['n2']
d1          =optim['d1']
d2          =optim['d2']
flag        =optim['modeltype']#1:smoothest, 2:flattest, 3:smallest
fmapr2b     =optim['fmapr2b']
fmapb2r     =optim['fmapb2r']

mapr2b = sparse.load_npz(fmapr2b)
mapb2r = sparse.load_npz(fmapb2r)
os.makedirs(os.path.join(fdir, 'regular'), exist_ok=True)

# choose one type of damping
if flag==1:
    print("generate the smoothest model")
    Lm = opt_laplace(n1, n2, d1, d2)
elif flag==2:
    print("generate the flattest model")
    Lm, Lm2, Lm3 = opt_deriv(n1, n2, d1, d2)
elif flag==3:
    print("generate the smallest model")
    Lm = opt_identity(n1, n2, d1, d2)
else:
    raise ValueError("\nSo far we haven't implemented for multiple damping types\n")
    
# output
foutL = os.path.join(fdir, 'regular', 'regularL1.npz')
sparse.save_npz(foutL, mapr2b * Lm * mapb2r)
if flag ==2:
    foutL = os.path.join(fdir, 'regular', 'regularL2.npz')
    sparse.save_npz(foutL, mapr2b * Lm2 * mapb2r)
    foutL = os.path.join(fdir, 'regular', 'regularL3.npz')
    sparse.save_npz(foutL, mapr2b * Lm3 * mapb2r)

print("Finished...", __file__)