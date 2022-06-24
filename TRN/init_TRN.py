#!/usr/bin/env python3
"""
copy extern files to optim file system at each iteration
"""
import os, sys, json
import numpy as np
from shutil import copyfile

#Input paras
fdir    =sys.argv[1] #'optim'
iterc   =int(sys.argv[2])
foptim  =os.path.join(fdir,'optim.json')
optim   =eval(open(foptim).read())
n       =optim['n']
iter0   =optim['iter0']
fmod    =optim['fmod']
fcost   =optim['fcost']
fgrad   =optim['fgrad']
fiter   ='iter'+str(iterc)
nwork   =10

#save arrays to optim file system
os.makedirs(os.path.join(fdir, fiter),exist_ok=True)
copyfile(fmod,  os.path.join(fdir, fiter, 'mod.bin'))
copyfile(fcost, os.path.join(fdir, fiter, 'fcost.dat'))
copyfile(fgrad, os.path.join(fdir, fiter, 'grad.bin'))


#allocate space for TN algorithm
"""
mod
grad
descent
descent_prev
residual
d
Hd
eisenvect
descent_rec
residual_rec
"""
worktn=np.zeros([nwork, n],np.float32)
mod =np.fromfile(os.path.join(fdir, fiter, 'mod.bin'), np.float32)
grad =np.fromfile(os.path.join(fdir, fiter, 'grad.bin'), np.float32)
worktn[0,:]=mod
worktn[1,:]=grad
worktn.tofile(os.path.join(fdir, 'worktn.bin'))

#update work infomation
optim['CG_phase']='INIT'
optim['nwork'] =nwork
optim['iterc'] =iterc
if iterc==iter0:
    optim['f0']=np.loadtxt(optim['fcost']).item()
optim['fk']=np.loadtxt(optim['fcost']).item()
optim['norm_grad']=np.linalg.norm(grad).item()
#print(optim)
#print(type(optim['fk']),optim['norm_grad'],type(optim['norm_grad']))
fout =open(foptim,'w')
fout.write(json.dumps(optim,indent=4))
fout.close()

#print("Finished...", __file__)