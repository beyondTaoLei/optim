#!/usr/bin/env python3
"""
the inversion example for Rosenbrock problem
Usage:
    cp model/init/init.vp model/inv/inv.vp
    cp model/init/init.vs model/inv/inv.vs
    cp model/init/init.rho model/inv/inv.rho
    
    python /home/tao/Nutstore_Files/works/seisplatform/mod_OPTIM/demo/surfacewave_phv_lsqr/S1_optim_LSQR.py optim 1 5
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
iter_start  =int(sys.argv[2]) # inversion starts from iter_start for this submission
iter_end    =int(sys.argv[3]) # inversion ends to iter_end for this submission
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
submit      =optim['submit']
prjroot     =optim['prjroot']
optimroot   =optim['optimroot']
py3         ='python3'
dirjob      = 'job'
lambda_idx  = 12

if submit == 0:
    iter_end = iter_start

fdroot = os.path.join(optimroot, 'demo', 'surfacewave_phv_lsqr', 'prob')
headers =[
    "prjroot='%s'"%prjroot,
    "fdroot='%s'"%fdroot,
    "optimroot='%s'"%optimroot,
    "py3='%s'"%py3
    ]

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
    for fnm in glob.glob('job/*'):
        os.remove(fnm)
    
    # calculate the difference and sensitivity matrix
    subjob= [
            '${py3} ' + os.path.join('${fdroot}', 'forward.py')+' '+fdir,
            '${py3} ' + os.path.join('${fdroot}', 'kernel.py')+' '+fdir +' 0',
            '${py3} ' + os.path.join('${fdroot}', 'calc_diff_sensmat.py')+' '+fdir
            ]
    jobnm = os.path.join(dirjob, 'fd.bash')
    generate_cmdfile(jobnm, subjob, headers, 0)
    if submit == 1: #run the tasks
        run_current_step([jobnm, ])
    
    # calculate the descent direction with LSQR
    subjob= [
            '${py3} ' + os.path.join('${optimroot}', 'regularization', 'generate_reg2d.py')+' '+fdir, 
            '${py3} ' + os.path.join('${optimroot}', 'LSQR', 'init_LSQR.py')+' '+fdir+' '+str(iterc),
            '${py3} ' + os.path.join('${optimroot}', 'LSQR', 'descent_LSQR.py')+' '+fdir
            ]
    jobnm = os.path.join(dirjob, 'descent.bash')
    generate_cmdfile(jobnm, subjob, headers, 0)
    if submit == 1: #run the tasks
        run_current_step([jobnm, ])
    
    # update the model
    subjob =[
            '${py3} ' + os.path.join('${fdroot}', 'select_Lcurve_result.py')+' '+fdir + ' '+str(lambda_idx),
            '${py3} ' + os.path.join('${fdroot}', 'update_model.py'+ ' '+ fdir)
            ]
    jobnm = os.path.join(dirjob, 'update_model.bash')
    generate_cmdfile(jobnm, subjob, headers, 0)
    if submit == 1: #run the tasks
        run_current_step([jobnm, ])
    
    ## save current information
    #print('finish optimization at iter. %4d of [%4d %4d] \n'%(iterc, iter_start, iter_end))
    