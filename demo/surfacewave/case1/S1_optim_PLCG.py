#!/usr/bin/env python3
"""
solve surface wave inversion problem with Preconditioning Linear Conjugate Gradient optimization
Usage:
    
"""
import os
import sys
import glob
import numpy  as np
import pandas as pd

def lcv_lambda(a0,step, num):
    a=np.zeros(num,np.float32)
    for i in range(num):
        a[i]=a0/(step**i)
    return a

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
py3         ='python3'
dirjob      = 'job'
lambda0 = 6.0e-3
lambdas = lcv_lambda(lambda0, 1.5, iter_end-iter_start+1)

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
  
  # calculate the difference and sensitivity matrix
  subjob= [
          '${mpiexec}' + ' -np ${NP} ${py3} ' + os.path.join('${fdroot}', 'forward.py')+' '+fdir,
          '${mpiexec}' + ' -np ${NP} ${py3} ' + os.path.join('${fdroot}', 'kernel.py')+' '+fdir +' 0',
          '${py3} ' + os.path.join('${fdroot}', 'calc_diff_sensmat.py')+' '+fdir
          ]
  jobnm = os.path.join(dirjob, 'fd.bash')
  generate_cmdfile(jobnm, subjob, headers, 0)
  idx += 1
  df.loc[idx] = [iterc, 1, 0, 0, 'bash ' + jobnm]
  
  # map the input grids to inversion grids
  subjob = ['${py3} ' + os.path.join('${fdroot}', 'media_fd2inv.py'+ ' '+ fdir), ]
  jobnm = os.path.join(dirjob, 'mediamap.bash')
  generate_cmdfile(jobnm, subjob, headers, 0)
  idx += 1
  df.loc[idx] = [iterc, 2, 0, 0, 'bash ' + jobnm]
  
  # calculate the descent direction with PLCG
  lam = lambdas[iterc-iter_start]
  subjob= [
          '${py3} ' + os.path.join('${fdroot}', 'update_optim.py'+ ' '+ fdir+ ' '+ 'lambda0'+ ' '+ str(lam)),
          '${py3} ' + os.path.join('${optimroot}', 'regularization', 'generate_reg2d.py')+' '+fdir, 
          '${py3} ' + os.path.join('${optimroot}', 'PLCG', 'init.py')+' '+fdir+' '+str(iterc),
          '${py3} ' + os.path.join('${optimroot}', 'PLCG', 'descent.py')+' '+fdir
          ]
  jobnm = os.path.join(dirjob, 'descent'+str(iterc)+'.bash')
  generate_cmdfile(jobnm, subjob, headers, 0)
  idx += 1
  df.loc[idx] = [iterc, 3, 0, 0, 'bash ' + jobnm]
  
  # update the model
  subjob =['${py3} ' + os.path.join('${fdroot}', 'update_model.py'+ ' '+ fdir), ]
  jobnm = os.path.join(dirjob, 'update_model.bash')
  generate_cmdfile(jobnm, subjob, headers, 0)
  idx += 1
  df.loc[idx] = [iterc, 4, 0, 0, 'bash ' + jobnm]

#write the tasks
jobNM = os.path.join(dirjob, 'job.csv')
df.to_csv(jobNM, index=False)
    