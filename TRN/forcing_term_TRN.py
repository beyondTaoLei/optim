#!/usr/bin/env python3
import numpy as np
from math import sqrt

def forcing_term_TRN(grad, residual_m1, norm_grad_m1, eta):
    """
    This file contains the routine for the computation  
    of forcing term following the first formula of      
    Eisenstat and Walker, see                           
    Choosing the forcing term in an Inexact Newton      
    method, Eisenstat and Walker, 1994,                 
    SIAM Journal on Scientific Computing 17 (1), 16-32  
    """
    
    #-----------------------------------------------------#
    # Computation of the forcing term eta following #
    # the formula                                         #
    #-----------------------------------------------------#
    eta_save=eta  
    eisenvect=grad - residual_m1 #current gradient - previous residual
    norm_eisenvect = np.linalg.norm(eisenvect)
    eta=norm_eisenvect/norm_grad_m1

    #-----------------------------------------------------#
    # Additional safeguard if eta is too large      #       
    #-----------------------------------------------------#
    eta_save_power=eta_save**((1.+sqrt(5.))/2.)
    if eta_save_power>0.1:
        eta=max(eta,eta_save_power)
    if eta>1.:
        eta=0.9
    
    return eta