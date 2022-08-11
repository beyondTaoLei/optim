#!/usr/bin/env python3
import os
import sys
import json
import numpy as np

#Input paras
fdirA       =sys.argv[1] #'optimvp'
fdirB       =sys.argv[2] #'optimvs'
foptim  =os.path.join(fdirA,'optim.json')
paras       =eval(open(foptim).read())
iterc       =paras['iterc']
fmodA       =os.path.join(fdirA, 'iter'+str(iterc), 'mod.bin')
fdescA      =os.path.join(fdirA, 'iter'+str(iterc), 'desc.bin')
fmodB       =os.path.join(fdirB, 'iter'+str(iterc), 'mod.bin')
fdescB      =os.path.join(fdirB, 'iter'+str(iterc), 'desc.bin')
fstep_estim =os.path.join(fdirA, 'test/step.json')
fmisfit     =os.path.join(fdirA, 'test/misfit1.dat')
fmod2A      =os.path.join(fdirA, 'test/m2A.bin')
fmod2B      =os.path.join(fdirA, 'test/m2B.bin')
#read m0 and delta
modA     =np.fromfile(fmodA, np.float32)
descA    =np.fromfile(fdescA, np.float32)
modB     =np.fromfile(fmodB, np.float32)
descB    =np.fromfile(fdescB, np.float32)
misfit1 =np.loadtxt(fmisfit).item()

sl_paras = eval(open(fstep_estim).read())
sl_paras['misfit1']=misfit1
if misfit1<sl_paras['misfit0']:
    sl_paras['alpha2']=sl_paras['beta1']*sl_paras['alphap']
else:
    sl_paras['alpha2']=sl_paras['beta2']*sl_paras['alphap']
#make new test model
alpha=sl_paras['alpha2']*sl_paras['alpha_factorA']
modA =modA+descA*alpha
modA.astype(np.float32).tofile(fmod2A)
alpha=sl_paras['alpha2']*sl_paras['alpha_factorB']
modB =modB+descB*alpha
modB.astype(np.float32).tofile(fmod2B)
#print('S1',modA)

fout = open(fstep_estim,'w')
fout.write(json.dumps(sl_paras,indent=4))
fout.close()

#print("Finished...", __file__)