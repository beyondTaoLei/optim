#!/usr/bin/env python3
"""
the inversion example for Rosenbrock problem
Usage:

"""
import os
import sys
import glob
import numpy  as np
import pandas as pd

def split_jobs(lst, n):
  # Calculate the size of each chunk (rounded up)
  chunk_size = -(-len(lst) // n)  # Using double negative for ceiling division
  subjobs = [lst[i:i + chunk_size] for i in range(0, len(lst), chunk_size)]
  return subjobs

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
  
#INPut paras
fdir        =sys.argv[1] # optim
iter_start  =int(sys.argv[2]) # inversion starts from iter_start for this submission
iter_end    =int(sys.argv[3]) # inversion ends to iter_end for this submission
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
prjroot     =optim['prjroot']
demoroot    =optim['demoroot']
optimroot   =optim['optimroot']
fsrc        =optim['fsrc']
ncores      =optim['ncores']
py3         ='python3'
dirjob      ='job'
fdir_loc    =os.path.join(prjroot, fdir)
headers = {
  "prjroot"   :prjroot,
  "demoroot"  :demoroot,
  "optimroot" :optimroot,
  "py3"       :py3
}

#source list
evn_info = pd.read_csv(fsrc)
shot_list = list(evn_info['station'])

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
  
  # run FD simulation
  jobs = []
  cmd_forward = '${py3}'+' '+ os.path.join('${demoroot}', 'prob', 'forward2d.py')+ ' ' + fdir_loc
  for evt in shot_list:
    cmd = 'cd '+os.path.join('${prjroot}', 'fdio', evt) + ' && ' + cmd_forward
    jobs.append(cmd)
  subjobs = split_jobs(jobs, ncores)
  for i in range(len(subjobs)):
    jobnm = os.path.join(dirjob, 'FD_p%02d.bash'%(i+1))
    name = os.path.basename(jobnm).split('.')[0]
    out = 'log/out.%s.it%d'%(name, iterc)
    err = 'log/err.%s.it%d'%(name, iterc)
    cmd = 'bash %s >%s 2>%s'%(jobnm, out, err)
    generate_cmdfile(jobnm, subjobs[i], headers)
    idx += 1
    df.loc[idx] = [iterc, 1, 0, 0, cmd, out, err]
  
  # calculate the residuals and the sensitivity matrix
  subjob = [
      '${py3} ' + os.path.join('${demoroot}', 'prob', 'calc_diff_sensmat.py')+ ' '+ fdir,
      ]
  jobnm = os.path.join(dirjob, 'calc_diff_sensmat.bash')
  name = os.path.basename(jobnm).split('.')[0]
  out = 'log/out.%s.it%d'%(name, iterc)
  err = 'log/err.%s.it%d'%(name, iterc)
  cmd = 'bash %s >%s 2>%s'%(jobnm, out, err)
  generate_cmdfile(jobnm, subjob, headers)
  idx += 1
  df.loc[idx] = [iterc, 2, 0, 0, cmd, out, err]
  
  # map the FD model to inversion model
  subjob = [
      '${py3} ' + os.path.join('${demoroot}', 'prob', 'media_fd2inv.py')+ ' '+ fdir,
      ]
  jobnm = os.path.join(dirjob, 'media_fd2inv.bash')
  name = os.path.basename(jobnm).split('.')[0]
  out = 'log/out.%s.it%d'%(name, iterc)
  err = 'log/err.%s.it%d'%(name, iterc)
  cmd = 'bash %s >%s 2>%s'%(jobnm, out, err)
  generate_cmdfile(jobnm, subjob, headers)
  idx += 1
  df.loc[idx] = [iterc, 3, 0, 0, cmd, out, err]
  
  # calculate the descent direction with LSQR
  subjob=[
      '${py3} ' + os.path.join('${optimroot}', 'regularization', 'generate_reg2d.py')+' '+fdir, 
      '${py3} ' + os.path.join('${optimroot}', 'LSQR', 'init.py')+' '+fdir+' '+str(iterc),
      '${py3} ' + os.path.join('${optimroot}', 'LSQR', 'descent.py')+' '+fdir
      ]
  jobnm = os.path.join(dirjob, 'descent_it'+str(iterc)+'.bash')
  name = os.path.basename(jobnm).split('.')[0]
  out = 'log/out.%s.it%d'%(name, iterc)
  err = 'log/err.%s.it%d'%(name, iterc)
  cmd = 'bash %s >%s 2>%s'%(jobnm, out, err)
  generate_cmdfile(jobnm, subjob, headers)
  idx += 1
  df.loc[idx] = [iterc, 4, 0, 0, cmd, out, err]
  
  # update the model and check the abort criterion
  subjob=[
      '${py3} ' + os.path.join('${demoroot}', 'prob', 'update_model.py')+ ' '+ fdir,
      '${py3} ' + os.path.join('${demoroot}', 'prob', 'check_stopping.py')+' '+fdir
      ]
  jobnm = os.path.join(dirjob, 'update_model.bash')
  name = os.path.basename(jobnm).split('.')[0]
  out = 'log/out.%s.it%d'%(name, iterc)
  err = 'log/err.%s.it%d'%(name, iterc)
  cmd = 'bash %s >%s 2>%s'%(jobnm, out, err)
  generate_cmdfile(jobnm, subjob, headers)
  idx += 1
  df.loc[idx] = [iterc, 5, 0, 0, cmd, out, err]

#write the tasks
jobNM = os.path.join(dirjob, 'job.csv')
df.to_csv(jobNM, index=False)