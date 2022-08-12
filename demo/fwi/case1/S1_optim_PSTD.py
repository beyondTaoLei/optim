#!/usr/bin/env python3
"""
the inversion example for Rosenbrock problem
Usage:

"""
import os
import sys
import glob
import time
import numpy  as np
import pandas as pd

def split_jobs(jobs, cell, ncores):
    """
    split jobs into subjobs for one lsf file
    """
    nsrc = len(jobs)
    num_in_group = ncores//cell
    subjobs = [jobs[i:i+num_in_group] for i in range(0, nsrc, num_in_group)]
    return subjobs

def generate_cmdfile(jobnm, joblist, queuenm, ncores, wtime, headers, parallel=1):
    """
    generate jobs with lsf system
    if necessary, change corresponding parameters.
    """
    hdconf ='''#!/bin/bash
    #BSUB -q queuenm
    #BSUB -n ncores
    #BSUB -R "span[ptile=ncores]"
    #BSUB -W wtime
    #BSUB -e %J.err
    #BSUB -o %J.out

    module load mpi/intel/2018.4
    module load intel/2018.4
    module load fftw/3.3.8
    module load matlab/2017b
    
    '''
    hdconf = hdconf.replace('queuenm',queuenm)
    hdconf = hdconf.replace('ncores','%d'%ncores)
    hdconf = hdconf.replace('wtime',wtime)
    hdconf = hdconf.replace("    ","")
    
    fp = open(jobnm,'w')
    a = fp.write(hdconf)
    for s in headers:
        a = fp.write("%s\n"%s)
    a = fp.write("\n")
    
    if parallel == 1:
        for s in joblist[:-1]:
            a = fp.write('%s&\n'% s) 
        a = fp.write('%s\n'% joblist[-1])
    else:
        for s in joblist:
            a = fp.write('%s\n'% s)
    fp.close()
    
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
            raise ValueError('fails to submit job: '+job)
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

#INPut paras
fdir        =sys.argv[1] # optim
iter_start  =int(sys.argv[2]) # inversion starts from iter_start for this submission
iter_end    =int(sys.argv[3]) # inversion ends to iter_end for this submission
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
submit      =optim['submit']
ncores      =optim['ncores']
mpiexec     =optim['mpiexec']
prjroot     =optim['prjroot']
optimroot   =optim['optimroot']
lnorm       =optim['lnorm']
fdexe       =optim['fdexe']
adexe       =optim['adexe']
NP          =optim['np']
fevn_info   =optim['fevn_info']
platform    =optim['platform']
json_test   =optim['json_test']
json_grad   =optim['json_grad']
py3         ='python3'
json_test = os.path.join('${prjroot}', optim['json_test'])
json_grad = os.path.join('${prjroot}', optim['json_grad'])

if submit == 0:
    iter_end = iter_start

fwiroot = os.path.join(optimroot, 'demo', 'fwi', 'prob')
headers = [
    "prjroot=%s"%prjroot,
    "fwiroot=%s"%fwiroot,
    "optimroot=%s"%optimroot,
    "fdexe=%s"%fdexe,
    "adexe=%s"%adexe,
    "mpiexec=%s"%mpiexec,
    "py3=%s"%py3,
    "NP=%d"%NP
    ]
if platform == 'lsf':
    queue_we = 'short'
    queue_single = 'smp'
    wtime = '2:00'
    dirjob = 'job'
else:
    raise ValueError('implement later')

#source list
evn_info = pd.read_csv(fevn_info)
test = np.array(evn_info['test'])
shot_list = list(evn_info['station'])
shot_test_list = list(evn_info['station'][test==1])
switch_event_dir = ['cd '+os.path.join('${prjroot}', 'fdio', evt) for evt in shot_list]
switch_test_event_dir = ['cd '+os.path.join('${prjroot}', 'fdio', evt) for evt in shot_test_list]

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
    for fnm in glob.glob('log/*.it%d'%iterc):
        os.remove(fnm)
    for fnm in glob.glob('job/*'):
        os.remove(fnm)
    
    # calculate the misfit and gradient
    jobs = []
    cmd_forward = '${mpiexec}' + ' -np ${NP} ${fdexe}'+ ' '+ json_grad
    cmd_data_process = 'mv seis_p_shot1.su syn_final.su' # replace with your data process
    cmd_misfit_adjsrc = '${py3} ' + os.path.join('${fwiroot}', 'calc_misfit_res.py 1 '+str(lnorm))
    cmd_adjoint = '${mpiexec}' + ' -np ${NP} ${adexe}'+ ' '+ json_grad
    for isrc in range(len(switch_event_dir)):
        cmd = [
            switch_event_dir[isrc],
            cmd_forward,
            cmd_data_process,
            cmd_misfit_adjsrc,
            cmd_adjoint
            ]
        jobs.append('&& '.join(cmd))
    subjobs = split_jobs(jobs, NP, ncores)
    for i in range(len(subjobs)):
        jobnm = os.path.join(dirjob, 'misfit_grad_part%d.lsf'%(i+1))
        generate_cmdfile(jobnm, subjobs[i], queue_we, ncores, wtime, headers, 1)
    if submit == 1: #run the tasks
        jobc = glob.glob(os.path.join(dirjob, 'misfit_grad_part*.lsf'))
        run_current_step(jobc, iterc)
    
    # sum up the misfit and gradient
    subjob = [
        '${py3} ' + os.path.join('${fwiroot}', 'sum_up_misfit.py')+' '+fdir+' all',
        '${py3} ' + os.path.join('${fwiroot}', 'sum_up_gradient.py')+' '+fdir,
        '${py3} ' + os.path.join('${fwiroot}', 'process_gradient.py')+' '+fdir,
        ]
    jobnm = os.path.join(dirjob, 'sum_misfit_grad.lsf')
    generate_cmdfile(jobnm, subjob, queue_single, 1, wtime, headers, 0)
    if submit == 1: #run the tasks
        run_current_step([jobnm, ], iterc)
    
    # calculate the descent direction
    subjob = [
        '${py3} ' + os.path.join('${optimroot}', 'PSTD', 'init.py')+' '+fdir+' '+str(iterc), 
        '${py3} ' + os.path.join('${optimroot}', 'PSTD', 'descent.py')+' '+fdir
        ]
    jobnm = os.path.join(dirjob, 'descent.lsf')
    generate_cmdfile(jobnm, subjob, queue_single, 1, wtime, headers, 0)
    if submit == 1: #run the tasks
        run_current_step([jobnm, ], iterc)
    
    # calculate the step length
    ## collect current misfit
    subjob = [
        '${py3} '+ os.path.join('${fwiroot}', 'sum_up_misfit.py')+' '+fdir+' test',
        '${py3} '+ os.path.join('${optimroot}', 'linesearch', 'S0_alpha.py')+' '+fdir
        ]
    jobnm = os.path.join(dirjob, 'misfit0.lsf')
    generate_cmdfile(jobnm, subjob, queue_single, 1, wtime, headers, 0)
    if submit == 1: #run the tasks
        run_current_step([jobnm, ], iterc)
    
    ## calculate misfit at 1st test
    ### calculate misfit
    jobs = []
    cmd_forward = '${mpiexec}' + ' -np ${NP} ${fdexe}'+ ' '+ json_test
    cmd_data_process = 'mv seis_p_shot1.su syn_final.su' # replace with your data process
    cmd_misfit_only = '${py3} ' + os.path.join('${fwiroot}', 'calc_misfit_res.py 0 '+str(lnorm))
    for isrc in range(len(switch_test_event_dir)):
        cmd = [
            switch_test_event_dir[isrc],
            cmd_forward,
            cmd_data_process,
            cmd_misfit_only
            ]
        jobs.append('&& '.join(cmd))
    subjobs = split_jobs(jobs, NP, ncores)
    for i in range(len(subjobs)):
        jobnm = os.path.join(dirjob, 'misfit1_part%d.lsf'%(i+1))
        generate_cmdfile(jobnm, subjobs[i], queue_we, ncores, wtime, headers, 1)
    if submit == 1: #run the tasks
        jobc = glob.glob(os.path.join(dirjob, 'misfit1_part*.lsf'))
        run_current_step(jobc, iterc)
    
    ### collect current misfit
    subjob = [
        '${py3} '+ os.path.join('${fwiroot}', 'sum_up_misfit.py')+' '+fdir+' test',
        '${py3} '+ os.path.join('${optimroot}', 'linesearch', 'S1_alpha.py')+' '+fdir
        ]
    jobnm = os.path.join(dirjob, 'misfit1.lsf')
    generate_cmdfile(jobnm, subjob, queue_single, 1, wtime, headers, 0)
    if submit == 1: #run the tasks
        run_current_step([jobnm, ], iterc)
    
    ## calculate misfit at 2nd test
    ### calculate misfit
    jobs = []
    cmd_forward = '${mpiexec}' + ' -np ${NP} ${fdexe}'+ ' '+ json_test
    cmd_data_process = 'mv seis_p_shot1.su syn_final.su' # replace with your data process
    cmd_misfit_only = '${py3} ' + os.path.join('${fwiroot}', 'calc_misfit_res.py 0 '+str(lnorm))
    for isrc in range(len(switch_test_event_dir)):
        cmd = [
            switch_test_event_dir[isrc],
            cmd_forward,
            cmd_data_process,
            cmd_misfit_only
            ]
        jobs.append('&& '.join(cmd))
    subjobs = split_jobs(jobs, NP, ncores)
    for i in range(len(subjobs)):
        jobnm = os.path.join(dirjob, 'misfit2_part%d.lsf'%(i+1))
        generate_cmdfile(jobnm, subjobs[i], queue_we, ncores, wtime, headers, 1)
    if submit == 1: #run the tasks
        jobc = glob.glob(os.path.join(dirjob, 'misfit2_part*.lsf'))
        run_current_step(jobc, iterc)
    
    ### collect current misfit
    subjob = [
        '${py3} '+ os.path.join('${fwiroot}', 'sum_up_misfit.py')+' '+fdir+' test',
        '${py3} '+ os.path.join('${optimroot}', 'linesearch', 'S2_alpha.py')+' '+fdir
        ]
    jobnm = os.path.join(dirjob, 'misfit2.lsf')
    generate_cmdfile(jobnm, subjob, queue_single, 1, wtime, headers, 0)
    if submit == 1: #run the tasks
        run_current_step([jobnm, ], iterc)
        
    # save current information
    print('finish optimization at iter. %4d of [%4d %4d] \n'%(iterc, iter_start, iter_end))
    