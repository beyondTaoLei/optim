#!/usr/bin/env python3
import numpy as np
import argparse, os, sys
from math import sqrt

#grad_preco,optim

def forcing_term_TRN(optim):
    """
    This file contains the routine for the computation  
    of forcing term following the first formula of      
    Eisenstat and Walker, see                           
    Choosing the forcing term in an Inexact Newton      
    method, Eisenstat and Walker, 1994,                 
    SIAM Journal on Scientific Computing 17 (1), 16-32  
    """
    
    #-----------------------------------------------------#
    # Computation of the forcing term optim.eta following #
    # the formula                                         #
    #-----------------------------------------------------#
    eta_save=optim.eta  
    optim.eisenvect[:]=optim.grad-optim.residual
    norm_eisenvect = np.linalg.norm(optim.eisenvect)
    optim.eta=norm_eisenvect/optim.norm_grad

    #-----------------------------------------------------#
    # Additional safeguard if optim.eta is too large      #       
    #-----------------------------------------------------#
    eta_save_power=eta_save**((1.+sqrt(5.))/2.)
    if eta_save_power>0.1:
        optim.eta=max(optim.eta,eta_save_power)
    if optim.eta>1.:
        optim.eta=0.9     