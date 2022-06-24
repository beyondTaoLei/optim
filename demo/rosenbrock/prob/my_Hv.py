#!/usr/bin/env python3
import os
import sys
import numpy as np
from rosenbrock import hess, hess2

#Input paras
fdir    =sys.argv[1]
foptim  =os.path.join(fdir,'optim.json')
optim   =eval(open(foptim).read())
fmod    =optim['fmod']
fvctr   =optim['fvctr']
fhess   =optim['fhess']

mod=np.fromfile(fmod, np.float32)
vctr=np.fromfile(fvctr, np.float32)
hes=hess2(mod, vctr)
#print(mod, vctr, hes)
hes.astype(np.float32).tofile(fhess)

print("Finished...", __file__)