#!/usr/bin/env python3
"""
solve surface wave inversion problem with LSQR optimization
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

#Input paras
fdir        =sys.argv[1] # optim
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
prjroot     =optim['prjroot']
demoroot    =optim['demoroot']
optimroot   =optim['optimroot']
mpiexep     =optim['mpiexep']
NP          =optim['ncores']
py3         ='python3'
dirjob      ='job'
iterc       =optim['iterc']
lambda0     =100.0
step        =4.0
num_lam     =20
lambdas = lcv_lambda(lambda0, step, num_lam)

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

"""
Inversion loop
"""
df = pd.DataFrame(columns=['iter', 'step', 'inner', 'status', 'cmd', 'outfile', 'errfile'])
idx = -1
for itest in range(num_lam):
  # print information
  print('*'*80)
  print('starts to generate the commands at iter. %4d of %4d'%(itest, itest)) 
  
  # calculate the descent direction with LSQR
  lam = lambdas[itest]
  fdesci = os.path.join(fdir, 'iter'+str(iterc), 'desc.bin')
  fdesco = os.path.join(fdir, 'iter'+str(iterc), 'L_curve'+str(itest)+'.bin')
  flambi = os.path.join(fdir, 'iter'+str(iterc), 'L_curve.dat')
  flambo = os.path.join(fdir, 'iter'+str(iterc), 'L_curve'+str(itest)+'.dat')
  subjob  = [
            '${py3} ' + os.path.join('${demoroot}', 'prob', 'update_optim.py'+ ' '+ fdir+ ' '+ 'lambda0'+ ' '+ str(lam)),
            '${py3} ' + os.path.join('${optimroot}', 'LSQR', 'descent.py')+' '+fdir,
            'mv ' + fdesci + ' ' + fdesco,
            'mv ' + flambi + ' ' + flambo
            ]
  jobnm = os.path.join(dirjob, 'descent_lam'+str(itest)+'.bash')
  name = os.path.basename(jobnm).split('.')[0]
  out = 'log/out.%s.it%d'%(name, iterc)
  err = 'log/err.%s.it%d'%(name, iterc)
  cmd = 'bash %s >%s 2>%s'%(jobnm, out, err)
  generate_cmdfile(jobnm, subjob, headers)
  idx += 1
  df.loc[idx] = [iterc, itest, 0, 0, cmd, out, err]

#write the tasks
jobNM = os.path.join(dirjob, 'lam.csv')
df.to_csv(jobNM, index=False)
    