#!/usr/bin/env python3
"""
the inversion example for Rosenbrock problem
Usage:
    python S1_optim_LBFGS.py optim_LBFGS 1 50
"""

import os
import sys
import json
import shutil

def print_info(p, key, iterc, level=1):
    """
    print information returned by os.system 
    """
    
    pre='*' * level + ' '
    info_fail = '%5s fails: %30s at iter. %4d'%(pre, key, iterc)
    info_succeed = '%5s succeeds: %30s at iter. %4d'%(pre, key, iterc)
    if p!=0:
        raise ValueError(info_fail)
    else:
        print(info_succeed)

#Input paras
fdir = sys.argv[1] # optim, as example
iter_start = int(sys.argv[2]) # inversion starts from iter_start for this submission
iter_end = int(sys.argv[3]) # inversion ends to iter_end for this submission

"""
Definition of functions
"""
cmd_misfit = 'python3 prob/my_misfit.py ' + fdir
cmd_gradient = 'python3 prob/my_grad.py ' + fdir
cmd_descent = ['python3 ../../LBFGS/init_LBFGS.py ' + fdir +' iterc',
               'python3 ../../LBFGS/descent_LBFGS.py ' + fdir
              ]
cmd_misfit_current = 'python3 ../../linesearch/S0_alpha.py ' + fdir
cmd_misfit_test1 = ['python3 prob/my_misfit.py ' + fdir,
                    'python3 ../../linesearch/S1_alpha.py ' + fdir
                   ]
cmd_misfit_test2 = ['python3 prob/my_misfit.py ' + fdir,
                    'python3 ../../linesearch/S2_alpha.py ' + fdir
                   ]
cmd_modelupdate = "echo 'nothing done'"

"""
Initialization for current stage
"""
# make directory
for iterc in range(iter_start, iter_end+1, 1):
    os.makedirs(os.path.join(fdir, 'iter'+str(iterc)), exist_ok=True)

"""
Inversion loop
"""
for iterc in range(iter_start, iter_end+1, 1):
    # print information
    print('*'*80)
    print('starts to run optimization at iter. %4d of [%4d %4d]'%(iterc, iter_start, iter_end))    
    
    # calculate the misfit
    p = os.system(cmd_misfit)
    print_info(p, '[misfit]', iterc, 1)
    
    # calculate the gradient
    p = os.system(cmd_gradient)
    print_info(p, '[gradient]', iterc, 1)
    
    # calculate the descent direction
    p = os.system(cmd_descent[0].replace('iterc', str(iterc)))
    print_info(p, '[descent: copy files]', iterc, 1)
    
    p = os.system(cmd_descent[1])
    print_info(p, '[descent: main]', iterc, 1)
    
    # calculate the step length
    ## collect current misfit
    p = os.system(cmd_misfit_current)
    print_info(p, '[0th step length]', iterc, 2)
    
    ## calculate misfit at 1st test
    p = os.system(cmd_misfit_test1[0])
    print_info(p, '[1st step lengt: misfit]', iterc, 2)
    
    p = os.system(cmd_misfit_test1[1])
    print_info(p, '[1st step length: alpha]', iterc, 2)
    
    ## calculate misfit at 2nd test
    p = os.system(cmd_misfit_test2[0])
    print_info(p, '[2nd step length: misfit]', iterc, 2)
    
    p = os.system(cmd_misfit_test2[1])
    print_info(p, '[2nd step length: alpha]', iterc, 2)
    
    # update the model
    p = os.system(cmd_modelupdate)
    print_info(p, '[model update]', iterc, 1)
    
    # save current information
    print('finish optimization at iter. %4d of [%4d %4d] \n'%(iterc, iter_start, iter_end))
    