#!/usr/bin/env python3
"""
the inversion example for Rosenbrock problem
Usage:
    cp model/true/true.vp model/inv/inv.vp
    cp model/true/true.vs model/inv/inv.vs
    cp model/true/true.rho model/inv/inv.rho
    
    python /home/tao/Nutstore_Files/works/seisplatform/mod_OPTIM/demo/firstarrival/S1_generate_data.py optim
"""
import os
import sys
import json
import glob
import time
import shutil
import numpy  as np
import pandas as pd

def split_jobs(jobs, cell, ncores):
    """
    split jobs into subjobs for one lsf file
    """
    nsrc=len(jobs)
    num_in_group = ncores//cell
    subjobs=[jobs[i:i+num_in_group] for i in range(0, nsrc, num_in_group)]
    return subjobs

def generate_cmdfile(jobnm, joblist, headers, parallel=1):
    """
    generate FD commands into file
    """
    fp = open(jobnm,'w')
    for s in headers:
        a=fp.write("%s\n"%s)
    a=fp.write("\n")
    
    if parallel == 1:
        for s in joblist[:-1]:
            a=fp.write('%s&\n'% s) 
        a=fp.write('%s\n'% joblist[-1])
    else:
        for s in joblist:
            a=fp.write('%s\n'% s)
    fp.close()
    
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

#Input paras
fdir        =sys.argv[1] # optim
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
ncores      =optim['ncores']
prjroot     =optim['prjroot']
optimroot   =optim['optimroot']
fevn_info   =optim['fevn_info']
py3         ='python3'
fdir_loc   =os.path.join('../..', fdir)
dirjob= 'job'

fdroot = os.path.join(optimroot, 'demo', 'firstarrival', 'prob')
headers =[
    "prjroot='%s'"%prjroot,
    "fdroot='%s'"%fdroot,
    "optimroot='%s'"%optimroot,
    "py3='%s'"%py3
    ]

#source list
evn_info=pd.read_csv(fevn_info)
shot_list=list(evn_info['station'])
switch_event_dir=['cd '+os.path.join('${prjroot}', 'fdio', evt) for evt in shot_list]

# print information
print('*'*80)
for fnm in glob.glob('job/*'):
    os.remove(fnm)

# calculate the misfit and gradient
jobs=[]
cmd_forward = '${py3}'+' '+ os.path.join('${fdroot}', 'forward2d.py '+ fdir_loc)
rename_file = 'mv time.csv obs.csv'
for isrc in range(len(switch_event_dir)):
    cmd = [switch_event_dir[isrc], cmd_forward, rename_file]
    jobs.append('&& '.join(cmd))

subjobs = split_jobs(jobs, 1, ncores)
for i in range(len(subjobs)):
    jobnm = os.path.join(dirjob, 'fd%02d.bash'%(i+1))
    generate_cmdfile(jobnm, subjobs[i], headers, 1)

jobc=sorted(glob.glob(os.path.join(dirjob, 'fd*.bash')))
run_current_step(jobc)
    