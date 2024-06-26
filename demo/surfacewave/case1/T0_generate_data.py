#!/usr/bin/env python3
"""
the inversion example for Rosenbrock problem
Usage:
    
"""
import os
import sys
import glob
import numpy  as np

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
prjroot     =optim['prjroot']
optimroot   =optim['optimroot']
mpiexec     =optim['mpiexec']
NP          =optim['ncores']
wdir        =optim['wdir']
odir        =optim['odir']
py3         ='python3'
dirjob      = 'job'

fdroot = os.path.join(optimroot, 'demo', 'surfacewave', 'prob')
headers =[
    "prjroot='%s'"%prjroot,
    "fdroot='%s'"%fdroot,
    "py3='%s'"%py3,
    "mpiexec='%s'"%mpiexec,
    "NP=%d"%NP
    ]

os.makedirs(odir, exist_ok=True)
"""
Inversion loop
"""
# print information
print('*'*80)
for fnm in glob.glob('job/*'):
    os.remove(fnm)

# calculate the difference and sensitivity matrix
subjob=['${mpiexec}' + ' -np ${NP} ${py3} ' + os.path.join('${fdroot}', 'forward.py')+' '+fdir,
        'cp -r ' + wdir + '/* ' + odir +'/'
       ]
jobnm = os.path.join(dirjob, 'fd0.bash')
generate_cmdfile(jobnm, subjob, headers, 0)
jobc=sorted(glob.glob(os.path.join(dirjob, 'fd0.bash')))
run_current_step(jobc)
    
    ## save current information
    #print('finish optimization at iter. %4d of [%4d %4d] \n'%(iterc, iter_start, iter_end))
    