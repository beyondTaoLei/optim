#!/usr/bin/env python3
"""
the inversion example for surface wave inversion problem
Usage:
    
"""
import os
import sys
import pandas as pd

#Input paras
fdir        =sys.argv[1] # optim
foptim      =os.path.join(fdir,'optim.json')
# optim       =eval(open(foptim).read())

#Input paras
dirjob = 'job'
jobNM = os.path.join(dirjob, 'job.csv')

"""
run the tasks
status:
 0: not start
-1: have failed
 1: be finished
"""
df = pd.read_csv(jobNM)
for idx in range(df.shape[0]):
  iter    = df.loc[idx]['iter']
  step    = df.loc[idx]['step']
  status  = df.loc[idx]['status']
  cmd     = df.loc[idx]['cmd']
  flag    = df.loc[idx]['flag']
  if status == 0 or status == -1:
    if status == 0:
      print('*** stage: ', idx, iter, step, status, cmd, 'will start')
    else:
      print('*** stage: ', idx, iter, step, status, cmd, 'will restart')
    # whether is truncated Newton approach
    if flag == 'innerCG':
      print('*** run inner loop ***')
      optim = eval(open(foptim).read()) # the lastest optim.json
      if optim['conv_CG'] == 1:
        df.loc[idx, 'status'] = 1
        df.to_csv(jobNM, index=False)
      else:
        p = os.system(cmd)
        if p!=0:
          df.loc[idx, 'status'] = -1
          df.to_csv(jobNM, index=False)
          raise ValueError('fails to submit job: '+cmd)
        else:
          df.loc[idx, 'status'] = 1
          df.to_csv(jobNM, index=False)
    else:
      p = os.system(cmd)
      if p!=0:
        df.loc[idx, 'status'] = -1
        df.to_csv(jobNM, index=False)
        raise ValueError('fails to submit job: '+cmd)
      else:
        df.loc[idx, 'status'] = 1
        df.to_csv(jobNM, index=False)
  elif status == 1:
    print('*** stage: ', idx, iter, step, status, cmd, 'has already finished')
  else:
    raise ValueError('we can not recognize the status code: '+status)