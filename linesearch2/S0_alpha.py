#!/usr/bin/env python3
import os, sys, json
import numpy as np
from angle_restart import angle_restart

#Input paras
fdirA       =sys.argv[1] #'optimvp'
fdirB       =sys.argv[2] #'optimvs'
foptim      =os.path.join(fdirA,'optim.json')
paras       =eval(open(foptim).read())
iterc       =paras['iterc']
iter0       =paras['iter0']
try_old_sl  =paras['try_old_sl']
eps_scaleA  =paras['eps_scale']
max_m0A     =paras['max_m0']
eps_scaleB  =paras['eps_scaleB']
max_m0B     =paras['max_m0B']
alpha_lb    =paras['alpha_lb']
alpha_ub    =paras['alpha_ub']
fmodA       =os.path.join(fdirA, 'iter'+str(iterc), 'mod.bin')
fgradA      =os.path.join(fdirA, 'iter'+str(iterc), 'grad.bin')
fdescA      =os.path.join(fdirA, 'iter'+str(iterc), 'desc.bin')
fmodB       =os.path.join(fdirB, 'iter'+str(iterc), 'mod.bin')
fgradB      =os.path.join(fdirB, 'iter'+str(iterc), 'grad.bin')
fdescB      =os.path.join(fdirB, 'iter'+str(iterc), 'desc.bin')
fstep_estim =os.path.join(fdirA, 'test/step.json')
fmisfit     =os.path.join(fdirA, 'test/misfit0.dat')
fmod2A      =os.path.join(fdirA, 'test/m1A.bin')
fmod2B      =os.path.join(fdirA, 'test/m1B.bin')
#read m0 and delta
modA     =np.fromfile(fmodA, np.float32)
gradA    =np.fromfile(fgradA, np.float32)
descA    =np.fromfile(fdescA, np.float32)
modB     =np.fromfile(fmodB, np.float32)
gradB    =np.fromfile(fgradB, np.float32)
descB    =np.fromfile(fdescB, np.float32)
misfit0 =np.loadtxt(fmisfit).item()

flag =angle_restart(gradA, descA)
if flag==1:
    descA =-gradA
    descA.tofile(fdescA)
flag =angle_restart(gradB, descB)
if flag==1:
    descB =-gradB
    descB.tofile(fdescB)

if try_old_sl==1:
    if iterc>iter0:
        print("use previous step length")
        sl_paras_old = eval(open(fstep_estim).read())
        alphap      =sl_paras_old['optim']
        #alpha_factorA=sl_paras_old['alpha_factorA']
    else:
        print("use initial step length")
        alphap=1.0
        #alpha_factorA=eps_scaleA*max_m0A/np.amax(np.fabs(descA))
else:
    print("use fixed step length")
    alphap=1.0
alpha_factorA=eps_scaleA*max_m0A/np.amax(np.fabs(descA))
alpha_factorB=eps_scaleB*max_m0B/np.amax(np.fabs(descB))
sl_paras={
    'max_m0A'       :max_m0A,    #useless
    'eps_scaleA'    :eps_scaleA,
    'max_m0B'       :max_m0B,    #useless
    'eps_scaleB'    :eps_scaleB,
    'beta1'         :2.0,
    'beta2'         :0.5,
    'beta3'         :0.05,
    'alpha_lb'      :alpha_lb,
    'alpha_ub'      :alpha_ub,
    'alphap'        :alphap,
    'alpha0'        :0.0,
    'alpha1'        :alphap,
    'alpha2'        :0.0,       #would be beta1*alpha1 or beta2*alpha1
    'misfit0'       :misfit0,
    'misfit1'       :0.0,       #will be updated later
    'misfit2'       :0.0,       #will be updated later
    'optim'         :alphap,    #will be updated finally
    'alpha_factorA'  :alpha_factorA,
    'alpha_factorB'  :alpha_factorB
    }
#make new test model
alpha=sl_paras['alpha1']*sl_paras['alpha_factorA']
modA =modA+descA*alpha
modA.astype(np.float32).tofile(fmod2A)
alpha=sl_paras['alpha1']*sl_paras['alpha_factorB']
modB =modB+descB*alpha
modB.astype(np.float32).tofile(fmod2B)
#print('S0',modA)

fout = open(fstep_estim,'w')
fout.write(json.dumps(sl_paras,indent=4))
fout.close()

#print("Finished...", __file__)