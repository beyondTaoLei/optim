#!/usr/bin/env python3
import os, sys, json
import numpy as np
from ployfit2 import ployfit2
from shutil import copyfile

#Input paras
fdir        =sys.argv[1] #'optim'
foptim      =os.path.join(fdir,'optim.json')
paras       =eval(open(foptim).read())
iterc       =paras['iterc']
mod_limits  =paras['mod_limits']
mod_lb      =paras['mod_lb']
mod_ub      =paras['mod_ub']
fmisfit     =paras['fcost']
fmodo       =paras['fmod']
fmodi       =os.path.join(fdir, 'iter'+str(iterc), 'mod.bin')
fdesc       =os.path.join(fdir, 'iter'+str(iterc), 'desc.bin')
fstep_estim =os.path.join(fdir, 'step.json')
#read m0 and delta
mod     =np.fromfile(fmodi, np.float32)
desc    =np.fromfile(fdesc, np.float32)
misfit2 =np.loadtxt(fmisfit).item()

sl_paras = eval(open(fstep_estim).read())
sl_paras['misfit2']=misfit2
epst=np.array([sl_paras['alpha0'], sl_paras['alpha1'], sl_paras['alpha2']])
misf=np.array([sl_paras['misfit0'], sl_paras['misfit1'], sl_paras['misfit2']])
sl_paras['optim']=ployfit2(epst,misf,sl_paras['alpha_lb'],sl_paras['alpha_ub'])
#make new test model
alpha=sl_paras['optim']*sl_paras['alpha_factor']
mod =mod+desc*alpha
#check the value limits
if mod_limits==1:
    print("check the model limits...")
    mod=np.clip(mod, mod_lb, mod_ub)
mod.astype(np.float32).tofile(fmodo)
#print('S2',mod)

fout = open(fstep_estim,'w')
fout.write(json.dumps(sl_paras,indent=4))
fout.close()
copyfile(fstep_estim,  os.path.join(fdir, 'iter'+str(iterc),'step.json'))
#debug
print('[%7.4f %7.4f %7.4f]'%(misf[0],misf[1],misf[2]), '[%7.4f %7.4f %7.4f] %7.4f'%(epst[0],epst[1],epst[2], sl_paras['optim']))
#print("eps_scale,max_m0,max_desc,alpha_factor alpha*eps_scale",sl_paras['eps_scale'],\
    #sl_paras['max_m0'], np.amax(np.fabs(desc)),sl_paras['alpha_factor'],sl_paras['optim']*sl_paras['eps_scale'])

#print("Finished...", __file__)
