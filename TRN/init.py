#!/usr/bin/env python3
"""
copy extern files to optim file system at each iteration
"""
import os
import sys
import json
import numpy as np
import shutil
from shutil import copyfile
from print_info import print_info
from forcing_term_TRN import forcing_term_TRN

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
nwork   =9

#save arrays to optim file system
fcur = os.path.join(fdir, fiter)
if os.path.exists(fcur):
    shutil.rmtree(fcur, ignore_errors=True)
os.makedirs(fcur,exist_ok=True)
copyfile(fmod,  os.path.join(fdir, fiter, 'mod.bin'))
copyfile(fcost, os.path.join(fdir, fiter, 'fcost.dat'))
copyfile(fgrad, os.path.join(fdir, fiter, 'grad.bin'))

# read model and gradient
mod =np.fromfile(os.path.join(fdir, fiter, 'mod.bin'), np.float32)
grad =np.fromfile(os.path.join(fdir, fiter, 'grad.bin'), np.float32)
#update optim
optim['fdir']=fdir
optim['CG_phase']='INIT'
optim['nwork'] =nwork
optim['iterc'] =iterc
optim['fdesclog']= os.path.join(fdir, fiter, 'desc.log')
fwork = os.path.join(fdir, 'worktn.bin')
# update eta
if iterc > iter0 and optim['eta_update'] == 1:
    worktn_old = np.fromfile(fwork, np.float32).reshape([nwork, n])
    residual_m1 = worktn_old[4,:]
    norm_grad_m1 = optim['norm_grad']
    eta0 = optim['eta']
    optim['eta'] = forcing_term_TRN(grad, residual_m1, norm_grad_m1, eta0)
#print information
str1 = '\n'+'*'*30+" iter = "+str(iterc)+' '+'*'*30
print_info(optim['fdesclog'], str1)
str1 = "iter: %3d, current eta: %.4f"%(iterc, optim['eta'])
print_info(optim['fdesclog'], str1)

#allocate space for TN algorithm
"""
mod
grad
descent
descent_prev
residual
d
Hd
descent_rec
residual_rec
"""
worktn=np.zeros([nwork, n],np.float32)
worktn[0,:]=mod
worktn[1,:]=grad
worktn.tofile(fwork)

#update optim infomation
optim['norm_grad']=np.linalg.norm(grad).item()
fout =open(foptim,'w')
fout.write(json.dumps(optim,indent=4))
fout.close()

#print("Finished...", __file__)