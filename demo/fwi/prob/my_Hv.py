#!/usr/bin/env python3
import os
import sys
import glob
import time
import numpy as np

def run_current_step(jobc, iterc):
    """
    submit the current task.
    """
    num = len(jobc)
    names, outs, errs = [], [], []
    for i in range(num):
        name = jobc[i].split('/')[-1][:-4]
        names.append(name)
        outs.append('log/out.%s.it%d'%(name, iterc))
        errs.append('log/err.%s.it%d'%(name, iterc))
    # delete the old files
    for fnm in outs + errs:
        try:
            os.remove(fnm)
        except OSError:
            pass
    # submit jobs
    for i in range(num):
        cmd = 'bsub -J %s <%s -o %s -e %s >>log/lsf.log'%(names[i], jobc[i], outs[i], errs[i])
        print(cmd)
        p = os.system(cmd)
        if p != 0:
            raise ValueError('fails to submit job: '+jobc)
    # whether to finish
    flag = 1
    while flag == 1:
        numo = sum([os.path.exists(fnm) for fnm in outs])
        nume = sum([os.path.exists(fnm) for fnm in errs])
        #check the error files
        if numo == nume == num: # files .out and .err exist
            errsize = sum([os.path.getsize(fnm) for fnm in errs])
            if errsize == 0:
                flag = 0 # finished
            else:
                flag = -1 # error happens
                raise ValueError('Error happens: ' + errs[0])
        else:
            print("*** waiting %s ***"%names[0])
            time.sleep(30)
            flag = 1 # not finished
    #return flag

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
    tag = 500*i+1500
    x2 = x + i * MCHPR * v
    x2.astype(np.float32).tofile(fmod)
    # command for misfit
    jobc = glob.glob(os.path.join(dirjob, 'misfit_grad_part*.lsf'))
    run_current_step(jobc, tag)
    jobnm = os.path.join(dirjob, 'sum_misfit_grad.lsf')
    run_current_step([jobnm, ], tag)
    grad = np.fromfile(fgrad, np.float32)
    grads.append(grad)

Hv=(grads[1]-grads[0])/(2.0*MCHPR)
Hv.astype(np.float32).tofile(fhess)

print("Finished...", __file__)