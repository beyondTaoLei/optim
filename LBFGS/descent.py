#!/usr/bin/env python3
"""
calculates the descent direction based on L-BFGS algorithm      

This routines is used by the l-BFGS routine to      
compute l-BFGS descent direction following the two  
recursion loop algorithm                             
Nocedal, 2nd edition, 2006 Algorithm 7.5 p. 179     
"""
    
import os
import sys
import json
import numpy as np

def descent_LBFGS_kernel(grad, sk, yk, cpt_lbfgs, l):
    
    descent=np.zeros(len(grad),np.float32)
                     
    borne_i=cpt_lbfgs-1

    #------------------------------------------#
    # SAFEGUARD                                #
    #------------------------------------------#
    norml2_sk = np.linalg.norm(sk[borne_i,:])
    norml2_yk = np.linalg.norm(yk[borne_i,:])

    if norml2_sk==0. or norml2_yk==0.:
        descent[:]=-1.*grad    
    else:
        #------------------------------------------#
        # First phase of the recursion loop        #
        #------------------------------------------#
        alpha = np.zeros(cpt_lbfgs)
        rho = np.zeros(cpt_lbfgs)
        q = grad.copy()
        
        #should be careful for range() function
        for i in range(borne_i+1): 
            ik=borne_i-i
            rho[ik] = np.inner(yk[ik,:],sk[ik,:])
            rho[ik]=1./rho[ik]
            alpha[ik] = np.inner(sk[ik,:],q)
            alpha[ik]=rho[ik]*alpha[ik]
            q -= alpha[ik]*yk[ik,:]
        
        gamma_num = np.inner(sk[borne_i,:],yk[borne_i,:])
        gamma_den=np.linalg.norm(yk[borne_i,:])
        #------------------------------------------#
        # Scaling by gamma                         #
        #------------------------------------------#
        gamma=gamma_num/(gamma_den**2)     
        descent[:]=gamma*q
        #------------------------------------------#
        # Second phase of the recursion loop       #
        #------------------------------------------#
        #should be careful for range() function
        for i in range(borne_i+1): 
            beta = np.inner(yk[i,:],descent)
            beta=rho[i]*beta
            descent += (alpha[i]-beta)*sk[i,:]
            
        descent[:] = -1.*descent
        #print('grad'+str(grad))
        #print('sk12'+str(sk[:3,:]))
        #print('yk'+str(yk[:3,:]))
        #print('desc'+str(descent))
    return descent

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