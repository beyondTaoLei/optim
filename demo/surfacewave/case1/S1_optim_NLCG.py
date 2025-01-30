#!/usr/bin/env python3
"""
solve surface wave inversion problem with NLCG optimization
Usage:
    
"""
import os
import sys
import glob
import numpy  as np
import pandas as pd

def generate_cmdfile(jobnm, cmdlist, headers):
  """
  generate jobs with bash script
  if necessary, change corresponding parameters.
  """
  prjroot = headers['prjroot']
  fp = open(jobnm,'w')
  # bash headers
  a = fp.write('#!/bin/bash \n\n')
  # environment variables
  for key, value in headers.items():
    s = key + '=' + str(value)
    a = fp.write("%s\n"%s)
  a = fp.write("\n")
  # write users' commands
  a = fp.write('commands=( \n')
  num = len(cmdlist)
  for i in range(num):
    a = fp.write('    "'+cmdlist[i] + '"\n')
  a = fp.write(') \n\n')
  # Run and check each Python command
  a = fp.write('for cmd in "${commands[@]}"; do\n')
  a = fp.write('    eval $cmd \n')
  a = fp.write('    if [ $? -ne 0 ]; then\n')
  a = fp.write('        echo "Exiting unsuccessfully." \n')
  a = fp.write('        echo "$cmd failed. Exiting..." >&2 \n')
  a = fp.write('        exit 1\n')
  a = fp.write('    fi\n')
  a = fp.write('done\n\n')
  a = fp.write('echo "Exiting successfully."\n')
  fp.close()

#Input paras
fdir        =sys.argv[1] # optim
iter_start  =int(sys.argv[2]) # inversion starts from iter_start for this submission
iter_end    =int(sys.argv[3]) # inversion ends to iter_end for this submission
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
prjroot     =optim['prjroot']
demoroot    =optim['demoroot']
optimroot   =optim['optimroot']
mpiexep     =optim['mpiexep']
NP          =optim['ncores']
py3         ='python3'
dirjob      ='job'

headers = {
  "prjroot"   :prjroot,
  "demoroot"  :demoroot,
  "optimroot" :optimroot,
  "mpiexep"   :mpiexep,
  "py3"       :py3,
  "NP"        :NP
}

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
df = pd.DataFrame(columns=['iter', 'step', 'inner', 'status', 'cmd', 'outfile', 'errfile'])
idx = -1
for iterc in range(iter_start, iter_end+1, 1):
  # print information
  print('*'*80)
  print('starts to generate the commands at iter. %4d of [%4d %4d]'%(iterc, iter_start, iter_end)) 
  
  # calculate the misfit
  subjob = [
            '${mpiexep}' + ' -np ${NP} ${py3} ' + os.path.join('${demoroot}', 'prob', 'forward.py')+' '+fdir,
            '${py3} ' + os.path.join('${demoroot}', 'prob', 'my_misfit.py')+' '+fdir
            ]
  jobnm = os.path.join(dirjob, 'fd.bash')
  name = os.path.basename(jobnm).split('.')[0]
  out = 'log/out.%s.it%d'%(name, iterc)
  err = 'log/err.%s.it%d'%(name, iterc)
  cmd = 'bash %s >%s 2>%s'%(jobnm, out, err)
  generate_cmdfile(jobnm, subjob, headers)
  idx += 1
  df.loc[idx] = [iterc, 1, 0, 0, cmd, out, err]
  
  # calculate the gradient
  subjob= [
          '${mpiexep}' + ' -np ${NP} ${py3} ' + os.path.join('${demoroot}', 'prob', 'kernel.py')+' '+fdir +' 0',
          '${py3} ' + os.path.join('${demoroot}', 'prob', 'my_grad.py')+' '+fdir
          ]
  jobnm = os.path.join(dirjob, 'grad.bash')
  name = os.path.basename(jobnm).split('.')[0]
  out = 'log/out.%s.it%d'%(name, iterc)
  err = 'log/err.%s.it%d'%(name, iterc)
  cmd = 'bash %s >%s 2>%s'%(jobnm, out, err)
  generate_cmdfile(jobnm, subjob, headers)
  idx += 1
  df.loc[idx] = [iterc, 2, 0, 0, cmd, out, err]

  # calculate the descent direction
  subjob= [
          '${py3} ' + os.path.join('${optimroot}', 'NLCG', 'init.py')+' '+fdir+' '+str(iterc),
          '${py3} ' + os.path.join('${optimroot}', 'NLCG', 'descent.py')+' '+fdir
          ]
  jobnm = os.path.join(dirjob, 'descent'+str(iterc)+'.bash')
  name = os.path.basename(jobnm).split('.')[0]
  out = 'log/out.%s.it%d'%(name, iterc)
  err = 'log/err.%s.it%d'%(name, iterc)
  cmd = 'bash %s >%s 2>%s'%(jobnm, out, err)
  generate_cmdfile(jobnm, subjob, headers)
  idx += 1
  df.loc[idx] = [iterc, 3, 0, 0, cmd, out, err]

  # calculate the step length
  ## collect current misfit and generate 1st test model
  subjob =['${py3} ' + os.path.join('${optimroot}', 'linesearch', 'S0_alpha.py')+' '+fdir, ]
  jobnm = os.path.join(dirjob, 'misfit0.bash')
  name = os.path.basename(jobnm).split('.')[0]
  out = 'log/out.%s.it%d'%(name, iterc)
  err = 'log/err.%s.it%d'%(name, iterc)
  cmd = 'bash %s >%s 2>%s'%(jobnm, out, err)
  generate_cmdfile(jobnm, subjob, headers)
  idx += 1
  df.loc[idx] = [iterc, 4, 0, 0, cmd, out, err]

  ## calculate misfit at 1st test
  subjob= [
          '${mpiexep}' + ' -np ${NP} ${py3} ' + os.path.join('${demoroot}', 'prob', 'forward.py')+' '+fdir,
          '${py3} ' + os.path.join('${demoroot}', 'prob', 'my_misfit.py')+' '+fdir
          ]
  jobnm = os.path.join(dirjob, 'misfit1.bash')
  name = os.path.basename(jobnm).split('.')[0]
  out = 'log/out.%s.it%d'%(name, iterc)
  err = 'log/err.%s.it%d'%(name, iterc)
  cmd = 'bash %s >%s 2>%s'%(jobnm, out, err)
  generate_cmdfile(jobnm, subjob, headers)
  idx += 1
  df.loc[idx] = [iterc, 5, 0, 0, cmd, out, err]

  ## collect current misfit and generate 2nd test model
  subjob= [
          '${py3} ' + os.path.join('${optimroot}', 'linesearch', 'S1_alpha.py')+' '+fdir
          ]
  jobnm = os.path.join(dirjob, 'misfit2.bash')
  name = os.path.basename(jobnm).split('.')[0]
  out = 'log/out.%s.it%d'%(name, iterc)
  err = 'log/err.%s.it%d'%(name, iterc)
  cmd = 'bash %s >%s 2>%s'%(jobnm, out, err)
  generate_cmdfile(jobnm, subjob, headers)
  idx += 1
  df.loc[idx] = [iterc, 6, 0, 0, cmd, out, err]
  
  ## calculate misfit at 2nd test
  subjob= [
          '${mpiexep}' + ' -np ${NP} ${py3} ' + os.path.join('${demoroot}', 'prob', 'forward.py')+' '+fdir,
          '${py3} ' + os.path.join('${demoroot}', 'prob', 'my_misfit.py')+' '+fdir
          ]
  jobnm = os.path.join(dirjob, 'misfit3.bash')
  name = os.path.basename(jobnm).split('.')[0]
  out = 'log/out.%s.it%d'%(name, iterc)
  err = 'log/err.%s.it%d'%(name, iterc)
  cmd = 'bash %s >%s 2>%s'%(jobnm, out, err)
  generate_cmdfile(jobnm, subjob, headers)
  idx += 1
  df.loc[idx] = [iterc, 7, 0, 0, cmd, out, err]

  ## collect current misfit and generate optimal model
  subjob= [
          '${py3} ' + os.path.join('${optimroot}', 'linesearch', 'S2_alpha.py')+' '+fdir,
          '${py3} ' + os.path.join('${demoroot}', 'prob', 'check_stopping.py')+' '+fdir
          ]
  jobnm = os.path.join(dirjob, 'misfit4.bash')
  name = os.path.basename(jobnm).split('.')[0]
  out = 'log/out.%s.it%d'%(name, iterc)
  err = 'log/err.%s.it%d'%(name, iterc)
  cmd = 'bash %s >%s 2>%s'%(jobnm, out, err)
  generate_cmdfile(jobnm, subjob, headers)
  idx += 1
  df.loc[idx] = [iterc, 8, 0, 0, cmd, out, err]

#write the tasks
jobNM = os.path.join(dirjob, 'job.csv')
df.to_csv(jobNM, index=False)