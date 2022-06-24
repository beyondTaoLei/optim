#!/usr/bin/env python3
import os, sys, json
import numpy as np
from ployfit2 import ployfit2
from shutil import copyfile

#Input paras
fdirA       =sys.argv[1] #'optimvp'
fdirB       =sys.argv[2] #'optimvs'
foptim      =os.path.join(fdirA,'optim.json')
paras       =eval(open(foptim).read())
iterc       =paras['iterc']
mod_limits  =paras['mod_limits']
mod_lbA     =paras['mod_lb']
mod_ubA     =paras['mod_ub']
mod_lbB     =paras['mod_lbB']
mod_ubB     =paras['mod_ubB']
fmodA       =os.path.join(fdirA, 'iter'+str(iterc), 'mod.bin')
fdescA      =os.path.join(fdirA, 'iter'+str(iterc), 'desc.bin')
fmodB       =os.path.join(fdirB, 'iter'+str(iterc), 'mod.bin')
fdescB      =os.path.join(fdirB, 'iter'+str(iterc), 'desc.bin')
fstep_estim =os.path.join(fdirA, 'test/step.json')
fmisfit     =os.path.join(fdirA, 'test/misfit2.dat')
fmod2A      =os.path.join(fdirA, 'test/mod_optimA.bin')
fmod2B      =os.path.join(fdirA, 'test/mod_optimB.bin')
#read m0 and delta
modA    =np.fromfile(fmodA, np.float32)
descA   =np.fromfile(fdescA, np.float32)
modB    =np.fromfile(fmodB, np.float32)
descB   =np.fromfile(fdescB, np.float32)
misfit2 =np.loadtxt(fmisfit).item()

sl_paras = eval(open(fstep_estim).read())
sl_paras['misfit2']=misfit2
epst=np.array([sl_paras['alpha0'], sl_paras['alpha1'], sl_paras['alpha2']])
misf=np.array([sl_paras['misfit0'], sl_paras['misfit1'], sl_paras['misfit2']])
sl_paras['optim']=ployfit2(epst,misf,sl_paras['alpha_lb'],sl_paras['alpha_ub'])
#make new test model
alpha=sl_paras['optim']*sl_paras['alpha_factorA']
modA =modA+descA*alpha
alpha=sl_paras['optim']*sl_paras['alpha_factorB']
modB =modB+descB*alpha
#check the value limits
if mod_limits==1:
    print("check the model limits...")
    modA=np.clip(modA, mod_lbA, mod_ubA)
    modB=np.clip(modB, mod_lbB, mod_ubB)
modA.astype(np.float32).tofile(fmod2A)
modB.astype(np.float32).tofile(fmod2B)
#print('S2',modA)

fout = open(fstep_estim,'w')
fout.write(json.dumps(sl_paras,indent=4))
fout.close()
copyfile(fstep_estim,  os.path.join(fdirA, 'iter'+str(iterc),'step.json'))
#debug
print('[%7.4f %7.4f %7.4f]'%(misf[0],misf[1],misf[2]), '[%7.4f %7.4f %7.4f] %7.4f'%(epst[0],epst[1],epst[2], sl_paras['optim']))
#print("eps_scale,max_m0,max_desc,alpha_factor alpha*eps_scale",sl_paras['eps_scale'],\
    #sl_paras['max_m0'], np.amax(np.fabs(descA)),sl_paras['alpha_factor'],sl_paras['optim']*sl_paras['eps_scale'])

#print("Finished...", __file__)
