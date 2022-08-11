#!/usr/bin/env python3
import os
import sys
import numpy as np
from scipy import sparse
from scipy.ndimage import gaussian_filter

#Input paras
fdir        =sys.argv[1] #optim
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
nx          =optim['n2']
ny          =optim['n1']
fwgtm       =optim['fwgtm']

istart = 17
taper = np.zeros([nx ,ny], np.float32)
taper[:, :istart] = 100.0
taper[:, istart:] = 1.0
taper = gaussian_filter(taper, 3/np.sqrt(8.0), mode='nearest')

Wm = sparse.diags(taper.flatten())
sparse.save_npz(fwgtm, Wm)

#print("Finished...", __file__)