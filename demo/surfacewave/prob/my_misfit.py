#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd
from scipy import sparse
from numpy import linalg as LA

#Input paras
fdir        =sys.argv[1]
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
wdir        =optim['wdir']
odir        =optim['odir']
fperiods    =optim['fperiods']
fstn        =optim['fstn']
fcost       =optim['fcost']
maxmode     =optim['maxmode']
modes       =maxmode + 1

# interpretation
df = pd.read_csv(fperiods)
nper = df.shape[0]

# read indexes of stations
df = pd.read_csv(fstn)
ntraces = df.shape[0]

b = np.zeros([modes, ntraces, nper], np.float32)
for mode in range(modes):
    # phase velocity
    fnm = os.path.join(wdir, 'vmap_M%1d.bin'%mode)
    phvsyn = np.fromfile(fnm, np.float32).reshape([ntraces, nper])
    fnm = os.path.join(odir, 'vmap_M%1d.bin'%mode)
    phvobs = np.fromfile(fnm, np.float32).reshape([ntraces, nper])
    diff = np.where((phvobs>0.0) &(phvsyn>0.0), phvobs-phvsyn, -1000.0)
    b[mode, :, :] = diff


# save the difference and sensmat
idxs = np.where(b>-900.0) #[modes, n2, nper]
b = b[idxs] # vector 
    
misfit=0.5*LA.norm(b)**2
fp = open(fcost,'w')
fp.write('%.6e\n'% misfit)
fp.close()

print("Finished...", __file__)