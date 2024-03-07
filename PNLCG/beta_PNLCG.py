#!/usr/bin/env python3

import math
import numpy as np

def beta_PNLCG(grad_prev, descent_prev, grad, grad_preco, powell):
    """
    calculates beta in nonlinear conjugate gradient     
    """
    tau1 = 0.2
    
    sk = grad-grad_prev
    
    gkpgk_prev  = np.inner(grad_prev,grad_prev)     #(gk'gk)        scalar product of previous gradient                   
    gkpgk       = np.inner(grad,grad_preco)         #(gk+1'gk+1)    scalar product of current gradient numerator of DY
    gkpgk2      = np.inner(grad_preco,grad_prev)    #(gk+1'gk)        scalar product of previous gradient

    gkpsk       = np.inner(grad,sk)                 # (gk+1'yk)     numerator of HS                                       
    skpk        = np.inner(sk,descent_prev)         # (dk'yk)       denominator of HS,DY                                  
    
    #beta=gkpgk/skpk
    beta_FR = gkpgk/gkpgk_prev
    beta_PRP = gkpsk/gkpgk_prev
    #beta_HS = gkpsk/skpk
    #beta_DY = gkpgk/skpk
    if math.isclose(np.linalg.norm(sk),0.0):
        beta_HS = 0.0
        beta_DY = 0.0
    else:
        beta_HS = gkpsk/skpk
        beta_DY = gkpgk/skpk

    #beta=beta_HS 
    #beta=beta_DY
    beta = beta_PRP
    if beta < 0.0: beta=0.0

    #------------------------------------------------#
    #         Powell restart condition               # 
    #------------------------------------------------#
    if powell==1 and gkpgk2/gkpgk_prev > tau1:
        beta = 0.0
        print('<gk+1 gk>/<gk gk>,tau1:%f,%f\n'%(gkpgk2/gkpgk_prev,tau1))
        print('RESTART CG BASED ON POWELL RESTART CONDITION\n')
    
    #------------------------------------------------------------#
    # Safeguard (may be useful in some cases)                    #
    #------------------------------------------------------------# 
    if(beta >= 1e5) or (beta <= -1e5):  
        beta=0.
    print('beta_FR,beta_PRP,beta_DY,beta_HS,beta=%f,%f,%f,%f,%f\n'%(beta_FR,beta_PRP,beta_DY,beta_HS,beta))
    
    return beta