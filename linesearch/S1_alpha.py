#!/usr/bin/env python3
import os
import sys
import json
import numpy as np

#Input paras
fdir        =sys.argv[1] #'optim'
foptim      =os.path.join(fdir,'optim.json')
paras       =eval(open(foptim).read())
iterc       =paras['iterc']
fmisfit     =paras['fcost']
fmodo       =paras['fmod']
fmodi       =os.path.join(fdir, 'iter'+str(iterc), 'mod.bin')
fdesc       =os.path.join(fdir, 'iter'+str(iterc), 'desc.bin')
fstep_estim =os.path.join(fdir, 'step.json')
#read m0 and delta
mod     =np.fromfile(fmodi, np.float32)
desc    =np.fromfile(fdesc, np.float32)
misfit1 =np.loadtxt(fmisfit).item()

sl_paras = eval(open(fstep_estim).read())
sl_paras['misfit1']=misfit1
if misfit1<sl_paras['misfit0']:
    sl_paras['alpha2']=sl_paras['beta1']*sl_paras['alphap']
else:
    sl_paras['alpha2']=sl_paras['beta2']*sl_paras['alphap']
#make new test model
alpha=sl_paras['alpha2']*sl_paras['alpha_factor']
mod =mod+desc*alpha
mod.astype(np.float32).tofile(fmodo)
#print('S1',mod)

fout = open(fstep_estim,'w')
fout.write(json.dumps(sl_paras,indent=4))
fout.close()

#print("Finished...", __file__)