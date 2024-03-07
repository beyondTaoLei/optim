#!/usr/bin/env python3
import os
import sys
import numpy as np
from rosenbrock import func

#Input paras
fdir    =sys.argv[1]
foptim  =os.path.join(fdir,'optim.json')
optim   =eval(open(foptim).read())
fmod    =optim['fmod']
fcost   =optim['fcost']

mod=np.fromfile(fmod, np.float32)
misfit=func(mod)
#print(mod, misfit)
fp = open(fcost,'w')
fp.write('%.6e\n'% misfit)
fp.close()

#print("Finished...", __file__)