#!/usr/bin/env python3
"""
calculates descent direction
"""
import os
import sys
import numpy as np

#Input paras
fdir    =sys.argv[1] #'optim'
foptim  =os.path.join(fdir,'optim.json')
optim   =eval(open(foptim).read())
iterc   =optim['iterc']

#read gradient 
grad =np.fromfile(os.path.join(fdir, 'iter'+str(iterc), 'grad.bin'), np.float32)
#generate the descent direction
descent =-1.*grad
#save descent
descent.tofile(os.path.join(fdir, 'iter'+str(iterc), 'desc.bin'))

#print("Finished...", __file__)