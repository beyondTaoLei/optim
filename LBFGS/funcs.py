#!/usr/bin/env python3
"""
calculates the descent direction based on L-BFGS algorithm      
*****************************************************
This routines is used by the l-BFGS routine to      
compute l-BFGS descent direction following the two  
recursion loop algorithm                             
Nocedal, 2nd edition, 2006 Algorithm 7.5 p. 179     
-----------------------------------------------------
"""
    
import os, sys
import numpy as np

def update_LBFGS(x, grad, sk, yk, cpt_lbfgs, l):
    """
    This routines is used by the l-BFGS routine to      
    compute the new pairs of models and gradient for    
    the inverse Hessian approximation                   
    Nocedal, 2nd edition, 2006 Algorithm 7.5 p. 179     
    """
    
    #print(str(cpt_lbfgs)+' sk1\n'+str(sk[:3,:]))
    if cpt_lbfgs <= l-1:
        #---------------------------------------------------#
        # if the number of stored pairs does not exceed the #
        # maximum value, then compute a new pair sk yk and  #
        # update the counter cpt_lbfgs                      #
        #---------------------------------------------------#
        #print('x'+str(x))
        sk[cpt_lbfgs,:]=x-sk[cpt_lbfgs,:]
        yk[cpt_lbfgs,:]=grad-yk[cpt_lbfgs,:]
        cpt_lbfgs=cpt_lbfgs+1
    else:
        #---------------------------------------------------#
        # otherwise, simply update the lth pair             #
        #---------------------------------------------------#
        sk[l-1,:]=x-sk[l-1,:]
        yk[l-1,:]=grad-yk[l-1,:]
    #print('x:'+str(x))    
    #print(str(cpt_lbfgs)+' sk2\n'+str(sk[:3,:])+'\n')
    return cpt_lbfgs

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

def save_LBFGS(x, grad, sk, yk, cpt_lbfgs, l):
    if cpt_lbfgs < l:
        #---------------------------------------------------#
        # if the number of stored pairs does not exceed the #
        # maximum value, then save x and grad               #
        #---------------------------------------------------#
        #print(cpt_lbfgs,sk.shape,x.shape)
        sk[cpt_lbfgs,:]=x.copy()
        yk[cpt_lbfgs,:]=grad.copy()   
    else:
        #---------------------------------------------------#
        # otherwise, erase the oldest pair and save the     #
        # new one (shift)                                   #
        #---------------------------------------------------#
        #should be careful for range() function
        for i in range(l-1):
            sk[i,:]=sk[i+1,:]
            yk[i,:]=yk[i+1,:]
            
        sk[l-1,:]=x.copy()
        yk[l-1,:]=grad.copy()

            