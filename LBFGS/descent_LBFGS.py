#!/usr/bin/env python3
"""
calculates the descent direction based on L-BFGS algorithm      

This routines is used by the l-BFGS routine to      
compute l-BFGS descent direction following the two  
recursion loop algorithm                             
Nocedal, 2nd edition, 2006 Algorithm 7.5 p. 179     
"""
    
import os, sys, json
import numpy as np
#from funcs import save_LBFGS
#from funcs import update_LBFGS
from funcs import descent_LBFGS_kernel

#Input paras
fdir    =sys.argv[1] #'optim'
foptim  =os.path.join(fdir,'optim.json')
optim = eval(open(foptim).read())
fmod        =optim['fmod']
iter0       =optim['iter0']
iterc       =optim['iterc']
cpt_lbfgs   =optim['cpt_lbfgs']
l           =optim['l']

#read gradient and model
x =np.fromfile(os.path.join(fdir, 'iter'+str(iterc), 'mod.bin'), np.float32)
grad =np.fromfile(os.path.join(fdir, 'iter'+str(iterc), 'grad.bin'), np.float32)
n=x.shape[0]

if iterc==iter0: #regard current point as initial one
    cpt_lbfgs=0
    descent =-1.*grad
    print("run L-BFGS and regard current point as initial one")
else: #can employe the previous info.;
    #LBFGS update
    sk = np.zeros([l,n], np.float32)
    yk = np.zeros([l,n], np.float32)
    cpt_lbfgs = min(iterc-iter0, l)
    x_prev = np.fromfile(os.path.join(fdir, 'iter'+str(iterc-1), 'mod.bin'), np.float32)
    grad_prev = np.fromfile(os.path.join(fdir, 'iter'+str(iterc-1), 'grad.bin'), np.float32)
    sk[cpt_lbfgs-1,:]=x-x_prev
    yk[cpt_lbfgs-1,:]=grad-grad_prev
    sk[cpt_lbfgs-1,:].tofile(os.path.join(fdir, 'iter'+str(iterc-1), 's.bin'))
    yk[cpt_lbfgs-1,:].tofile(os.path.join(fdir, 'iter'+str(iterc-1), 'y.bin'))
    if cpt_lbfgs > 1:
        for i in range(cpt_lbfgs-2, -1, -1):
            sk[i,:] = np.fromfile(os.path.join(fdir, 'iter'+str(iterc-cpt_lbfgs+i), 's.bin'), np.float32)
            yk[i,:] = np.fromfile(os.path.join(fdir, 'iter'+str(iterc-cpt_lbfgs+i), 'y.bin'), np.float32)
    #Computation of the new descent direction
    descent =descent_LBFGS_kernel(grad, sk, yk, cpt_lbfgs, l)
    print("run L-BFGS and employe the previous info.")

#save descent
descent.tofile(os.path.join(fdir, 'iter'+str(iterc), 'desc.bin'))

#update work infomation
optim['cpt_lbfgs']=cpt_lbfgs
fout = open(foptim,'w')
fout.write(json.dumps(optim,indent=4))
fout.close()

#print("Finished...", __file__)