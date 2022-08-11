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
    for job in jobc:
        name = job.split('/')[-1][:-4]
        out = 'log/out.%s.it%d'%(name, iterc)
        err = 'log/err.%s.it%d'%(name, iterc)
        cmd = 'bsub -J %s <%s -K -o %s -e %s&'%(name, job, out, err)
        print(cmd)
        p = os.system(cmd)
        if p != 0:
            raise ValueError('fails to submit job: '+job)
    # whether to finish
    while check_finished(jobc, iterc) == 1:
        print("*** waiting %s ***"%jobc[0][4:7])
        time.sleep(30)

def check_finished(jobc, iterc):
    """
    check if current step is finished
    """
    numo, nume = 0, 0
    outs, errs = [], []
    for job in jobc:
        name = job.split('/')[-1][:-4]
        outs.append('log/out.%s.it%d'%(name, iterc))
        errs.append('log/err.%s.it%d'%(name, iterc))
        if os.path.exists(outs[-1]):
            numo += 1
        if os.path.exists(errs[-1]):
            nume += 1
    
    #check the error files
    if numo == nume == len(jobc):
        err_size = 0
        for ferr in errs:
            err_size += os.path.getsize(ferr)
        if err_size == 0:
            return 0
        else:
            info_fail = '!!! fails: %30s at iter. %4d !!!'%(jobc[0], iterc)
            raise ValueError(info_fail)
    else:
        return 1

#Input paras
fdir        =sys.argv[1] # optim
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
submit      =optim['submit']
ncores      =optim['ncores']
mpiexec     =optim['mpiexec']
prjroot     =optim['prjroot']
optimroot   =optim['optimroot']
fdexe       =optim['fdexe']
adexe       =optim['adexe']
NP          =optim['np']
fevn_info   =optim['fevn_info']
platform    =optim['platform']
json_test   =optim['json_test']
py3         ='python3'
json_test = os.path.join('${prjroot}', optim['json_test'])
iterc = 0

fwiroot = os.path.join(optimroot, 'demo', 'FWI_acoustic', 'prob')
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
shot_list = list(evn_info['station'])
switch_event_dir = ['cd '+os.path.join('${prjroot}', 'fdio', evt) for evt in shot_list]

"""
Inversion loop
"""
# print information
print('*'*80)
for fnm in glob.glob('log/*.it%d'%iterc):
    os.remove(fnm)
for fnm in glob.glob('job/*'):
    os.remove(fnm)

# calculate the misfit and gradient
jobs = []
for isrc in range(len(switch_event_dir)):
    cmd = [
        switch_event_dir[isrc],
        '${mpiexec}' + ' -np ${NP} ${fdexe}'+ ' '+ json_test,
        'mv seis_p_shot1.su '+ os.path.join('${prjroot}', 'data', '%s.su'%shot_list[isrc]) # replace with your data process,
        ]
    jobs.append('&& '.join(cmd))
subjobs = split_jobs(jobs, NP, ncores)
for i in range(len(subjobs)):
    jobnm = os.path.join(dirjob, 'fd%d.lsf'%(i+1))
    generate_cmdfile(jobnm, subjobs[i], queue_we, ncores, wtime, headers, 1)
if submit == 1: #run the tasks
    jobc = glob.glob(os.path.join(dirjob, 'fd*.lsf'))
    run_current_step(jobc, iterc)
    
print("Finished...", __file__)