#!/usr/bin/env python3
import os
import sys
import numpy as np
from rosenbrock import grad

#Input paras
fdir    =sys.argv[1]
foptim  =os.path.join(fdir,'optim.json')
optim   =eval(open(foptim).read())
fmod    =optim['fmod']
fgrad   =optim['fgrad']

mod=np.fromfile(fmod, np.float32)
grad=grad(mod)
#print(mod, grad)
grad.astype(np.float32).tofile(fgrad)

#print("Finished...", __file__)