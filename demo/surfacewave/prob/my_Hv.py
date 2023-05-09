#!/usr/bin/env python3
import os
import sys
import numpy as np

def run_current_step(jobc):
    """
    submit the current task.
    """
    for job in jobc:
        cmd = 'bash ' + job
        print(cmd)
        p = os.system(cmd)
        if p!=0:
            raise ValueError('fails to submit job: '+job)

def model_perturbation_deta(x, v, flag):
    #/* master node will call this function */
    epsilon=1.0E-4; #1.D-5 1.D-6
    
    if flag == 1:
        MCHPR = 1.0
    else:
        norm_x = np.linalg.norm(x)
        norm_v = np.linalg.norm(v)
        if norm_v<=1.0E-10: # norm_v is close to zero
            print("norm_v=%e ",norm_v)
            msg = "!!!1)For Newton method you should check the codes!!!"
            raise ValueError(msg)
        else:
            MCHPR=2*epsilon*(1.0+norm_x)/norm_v
            print("norm_x=%e,norm_v=%e,MCHPR=%e"%(norm_x,norm_v,MCHPR))
            if MCHPR <= 1.0E-12:
                print("MCHPR=%e "%MCHPR)
                msg = "!!!2)For Newton method you should check the codes!!!"
                raise ValueError(msg)
    return MCHPR

def Hv_FD(fmod, fvctr, fgrad, job_list_grad):
    """
    2nd order central finite difference
    """
    x = np.fromfile(fmod, np.float32)
    v = np.fromfile(fvctr, np.float32)
    MCHPR = model_perturbation_deta(x, v, 0)
    grads = []
    for i in [-1, 1]:
        x2 = x + i * MCHPR * v
        x2.astype(np.float32).tofile(fmod)
        # command for misfit
        run_current_step(job_list_grad)
        grad = np.fromfile(fgrad, np.float32)
        grads.append(grad)
    
    Hv=(grads[1]-grads[0])/(2.0*MCHPR)
    return Hv
    

#Input paras
fdir    =sys.argv[1]
foptim  =os.path.join(fdir,'optim.json')
optim   =eval(open(foptim).read())
fmod    =optim['fmod']
fgrad   =optim['fgrad']
fvctr   =optim['fvctr']
fhess   =optim['fhess']
dirjob  = 'job'
jobnm1 = os.path.join(dirjob, 'fd.bash')
jobnm2 = os.path.join(dirjob, 'grad.bash')
job_list_grad=[jobnm1, jobnm2]
Hv = Hv_FD(fmod, fvctr, fgrad, job_list_grad)
Hv.astype(np.float32).tofile(fhess)

print("Finished...", __file__)