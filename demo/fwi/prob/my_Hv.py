#!/usr/bin/env python3
import os
import sys
import glob
import time
import subprocess
import numpy as np

def run_current_step(jobc, iterc):
    """
    submit the current task.
    """
    outs, errs = [], []
    for job in jobc:
        name = job.split('/')[-1][:-4]
        out = 'log/out.%s.it%d'%(name, iterc)
        err = 'log/err.%s.it%d'%(name, iterc)
        outs.append(out)
        errs.append(err)
        cmd = 'bsub -J %s <%s -o %s -e %s >>log/lsf.log'%(name, job, out, err)
        print(cmd)
        p = os.system(cmd)
        if p != 0:
            raise ValueError('fails to submit job: '+job)
    # whether to finish
    num = len(outs)
    flag = 1
    while flag == 1:
        numo, nume = 0, 0
        for i in range(num):
            if os.path.exists(outs[i]):
                numo += 1
            if os.path.exists(errs[i]):
                nume += 1
        #check the error files
        if numo == nume == num:
            for ferr in errs:
                if os.path.getsize(ferr) != 0:
                    raise ValueError('!!! errors: ' + ferr + ' !!!')
            flag = 0 # finished
        else:
            print("*** waiting %s ***"%name)
            time.sleep(30)
            flag = 1 # not finished

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

#Input paras
fdir    =sys.argv[1]
foptim  =os.path.join(fdir,'optim.json')
optim   =eval(open(foptim).read())
fmod    =optim['fmod']
fgrad   =optim['fgrad']
fvctr   =optim['fvctr']
fhess   =optim['fhess']
iterc   =optim['iterc']
dirjob  = 'job'
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
    # rm the log
    for fnm in glob.glob('log/*.it%d'%(iterc+101+i)):
        os.remove(fnm)
    # command for misfit
    jobc = glob.glob(os.path.join(dirjob, 'misfit_grad_part*.lsf'))
    run_current_step(jobc, iterc+101+i)
    jobnm = os.path.join(dirjob, 'sum_misfit_grad.lsf')
    run_current_step([jobnm, ], iterc+101+i)
    grad = np.fromfile(fgrad, np.float32)
    grads.append(grad)

Hv=(grads[1]-grads[0])/(2.0*MCHPR)
Hv.astype(np.float32).tofile(fhess)

print("Finished...", __file__)