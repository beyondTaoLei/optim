#!/usr/bin/env python3
"""
the inversion example for surface wave inversion problem
Usage:
    
"""
import os
import sys
import glob
import numpy  as np
import pandas as pd

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
    
#Input paras
fdir        =sys.argv[1] # optim
iter_start  =int(sys.argv[2]) # inversion starts from iter_start for this submission
iter_end    =int(sys.argv[3]) # inversion ends to iter_end for this submission
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
prjroot     =optim['prjroot']
optimroot   =optim['optimroot']
mpiexec     =optim['mpiexec']
NP          =optim['ncores']
niter_max_CG=optim['niter_max_CG']
py3         ='python3'
dirjob      = 'job'

fdroot = os.path.join(optimroot, 'demo', 'surfacewave', 'prob')
headers =[
        "prjroot='%s'"%prjroot,
        "fdroot='%s'"%fdroot,
        "optimroot='%s'"%optimroot,
        "py3='%s'"%py3,
        "mpiexec='%s'"%mpiexec,
        "NP=%d"%NP
        ]

"""
Initialization for current stage
"""
# make directory
for fnm in glob.glob('job/*'):
  os.remove(fnm)
for iterc in range(iter_start, iter_end+1, 1):
  os.makedirs(os.path.join(fdir, 'iter'+str(iterc)), exist_ok=True)


"""
Inversion loop
"""
df = pd.DataFrame(columns=['iter', 'step', 'status', 'flag', 'cmd'])
idx = -1
for iterc in range(iter_start, iter_end+1, 1):
  # print information
  print('*'*80)
  print('generate the commands at iter. %4d of [%4d %4d]'%(iterc, iter_start, iter_end))    
  
  # calculate the misfit
  subjob= [
          '${mpiexec}' + ' -np ${NP} ${py3} ' + os.path.join('${fdroot}', 'forward.py')+' '+fdir,
          '${py3} ' + os.path.join('${fdroot}', 'my_misfit.py')+' '+fdir
          ]
  jobnm = os.path.join(dirjob, 'fd.bash')
  generate_cmdfile(jobnm, subjob, headers, 0)
  idx += 1
  df.loc[idx] = [iterc, 1, 0, 0, 'bash ' + jobnm]
  
  # calculate the gradient
  subjob= [
          '${mpiexec}' + ' -np ${NP} ${py3} ' + os.path.join('${fdroot}', 'kernel.py')+' '+fdir +' 0',
          '${py3} ' + os.path.join('${fdroot}', 'my_grad.py')+' '+fdir
          ]
  jobnm = os.path.join(dirjob, 'grad.bash')
  generate_cmdfile(jobnm, subjob, headers, 0)
  idx += 1
  df.loc[idx] = [iterc, 2, 0, 0, 'bash ' + jobnm]
  
  # calculate the descent direction
  ##########################################
  subjob= [
          '${py3} ' + os.path.join('${optimroot}', 'TRN', 'init.py')+' '+fdir+' '+str(iterc),
          '${py3} ' + os.path.join('${optimroot}', 'TRN', 'descent.py')+' '+fdir,
          '${py3} ' + os.path.join('${fdroot}', 'my_Hv.py')+' '+fdir
          ]
  jobnm = os.path.join(dirjob, 'descent'+str(iterc)+'.bash')
  generate_cmdfile(jobnm, subjob[:2], headers, 0)
  idx += 1
  df.loc[idx] = [iterc, 3, 0, 0, 'bash ' + jobnm]
  for icg in range(niter_max_CG):
    jobnm = os.path.join(dirjob, 'descenticg.bash')
    generate_cmdfile(jobnm, [subjob[2], subjob[1]], headers, 0)
    idx += 1
    df.loc[idx] = [iterc, 4, 0, 'innerCG', 'bash ' + jobnm]
  
  # calculate the step length
  ## collect current misfit and generate 1st test model
  subjob =['${py3} ' + os.path.join('${optimroot}', 'linesearch', 'S0_alpha.py')+' '+fdir, ]
  jobnm = os.path.join(dirjob, 'misfit0.bash')
  generate_cmdfile(jobnm, subjob, headers, 0)
  idx += 1
  df.loc[idx] = [iterc, 5, 0, 0, 'bash ' + jobnm]
  
  ## calculate misfit at 1st test
  subjob= [
          '${mpiexec}' + ' -np ${NP} ${py3} ' + os.path.join('${fdroot}', 'forward.py')+' '+fdir,
          '${py3} ' + os.path.join('${fdroot}', 'my_misfit.py')+' '+fdir
          ]
  jobnm = os.path.join(dirjob, 'misfit1.bash')
  generate_cmdfile(jobnm, subjob, headers, 0)
  idx += 1
  df.loc[idx] = [iterc, 6, 0, 0, 'bash ' + jobnm]

  ## collect current misfit and generate 2nd test model
  subjob= [
          '${py3} ' + os.path.join('${optimroot}', 'linesearch', 'S1_alpha.py')+' '+fdir
          ]
  jobnm = os.path.join(dirjob, 'misfit2.bash')
  generate_cmdfile(jobnm, subjob, headers, 0)
  idx += 1
  df.loc[idx] = [iterc, 7, 0, 0, 'bash ' + jobnm]
  
  ## calculate misfit at 2nd test
  subjob= [
          '${mpiexec}' + ' -np ${NP} ${py3} ' + os.path.join('${fdroot}', 'forward.py')+' '+fdir,
          '${py3} ' + os.path.join('${fdroot}', 'my_misfit.py')+' '+fdir
          ]
  jobnm = os.path.join(dirjob, 'misfit3.bash')
  generate_cmdfile(jobnm, subjob, headers, 0)
  idx += 1
  df.loc[idx] = [iterc, 8, 0, 0, 'bash ' + jobnm]

  ## collect current misfit and generate optimal model
  subjob= [
          '${py3} ' + os.path.join('${optimroot}', 'linesearch', 'S2_alpha.py')+' '+fdir
          ]
  jobnm = os.path.join(dirjob, 'misfit4.bash')
  generate_cmdfile(jobnm, subjob, headers, 0)
  idx += 1
  df.loc[idx] = [iterc, 9, 0, 0, 'bash ' + jobnm]

#write the tasks
jobNM = os.path.join(dirjob, 'job.csv')
df.to_csv(jobNM, index=False)