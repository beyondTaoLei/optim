#!/usr/bin/env python3
"""
calculates descent direction
"""
import os, sys
import numpy as np
from beta_PNLCG import beta_PNLCG

#Input paras
fdir    =sys.argv[1] #'optim'
foptim  =os.path.join(fdir,'optim.json')
optim   =eval(open(foptim).read())
iter0   =optim['iter0']
iterc   =optim['iterc']
powell  =optim['powell']

#read gradient 
grad =np.fromfile(os.path.join(fdir, 'iter'+str(iterc), 'grad.bin'), np.float32)
grad_preco =grad.copy()
#calculate the beta for conjugate gradient method
if iterc>iter0: #can employe the previous info.;
    fprev=os.path.join(fdir, 'iter'+str(iterc-1))
    if os.path.exists(os.path.join(fprev,'desc.bin')):
        #print("I will generate the conjugate gradient direction")
        grad_prev =np.fromfile(os.path.join(fprev,'grad.bin'), np.float32)
        descent_prev =np.fromfile(os.path.join(fprev,'desc.bin'), np.float32)
        beta =beta_PNLCG(grad_prev, descent_prev, grad, grad_preco, powell)
        descent =-1.*grad_preco + beta*descent_prev
    else: #regard current point as initial one
        descent =-1.*grad_preco
        raise ValueError("I do not have previous gradient and descent direction")
else: #regard current point as initial one
    print("I will generate the steepest descent direction")
    descent =-1.*grad_preco
#save descent
descent.tofile(os.path.join(fdir, 'iter'+str(iterc), 'desc.bin'))

#print("Finished...", __file__)
    
    