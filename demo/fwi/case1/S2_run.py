#!/usr/bin/env python3
"""
the inversion example for surface wave inversion problem
Usage:
    
"""
import os
import sys
import time
import numpy as np
import pandas as pd

def check_status_end(fnm, flag_finished, flag_failed, nL=10):
  a = 0
  if os.path.exists(fnm):
    with open(fnm, 'r', encoding='utf-8') as fp:
      lines = fp.readlines()
      if lines:
        for line in lines[-nL:]:
          if flag_finished in line:
            a = 1
            break
          elif flag_failed in line:
            a = -1
            break
  return a

#Input paras
fdir          =sys.argv[1] # optim
foptim        =os.path.join(fdir,'optim.json')
flag_finished ='Exiting successfully'
flag_failed   ='Exiting unsuccessfully'
time_sleep    =60 # seconds
nL            =10 # number of lines to check in the end of output file

#Input paras
dirjob = 'job'
jobNM = os.path.join(dirjob, 'job.csv')

"""
['iter', 'step', 'inner', 'status', 'cmd', 'outfile', 'errfile']
run the tasks
status:
 0: not start
-1: have failed
 1: be finished
"""

df = pd.read_csv(jobNM)
iter_start  = df['iter'].min()
iter_end    = df['iter'].max()
step_start  = df['step'].min()
step_end    = df['step'].max()
for iterc in range(iter_start, iter_end+1, 1):
  # first check 'abort'
  optim = eval(open(foptim).read())
  if optim['abort'] == 1:
    print('Current stage is over!')
    break
  else:
    print('Current stage is running!')
  # filter the current iteration
  dfit = df[df['iter']==iterc]
  for istep in range(step_start, step_end+1, 1):
    ## filter the current step
    dfst = dfit[dfit['step'] == istep]
    inner_start = dfst['inner'].min()
    inner_end = dfst['inner'].max()
    for icg in range(inner_start, inner_end+1, 1):
      ### filter the current inner loop
      dfin = dfst[dfst['inner'] == icg]
      if inner_end != 0:
        ### second check 'conv_CG'
        optim = eval(open(foptim).read())
        if optim['conv_CG'] == 1:
          print('*** inner loop: ', istep, icg, ' has been finished ***')
          df.loc[dfst.index, 'status'] = 1
          df.to_csv(jobNM, index=False)
          break
      #### filter the current job
      dfjob = dfin[dfin['status'] != 1]
      jobc = dfjob['cmd'].to_list()
      outs = dfjob['outfile'].to_list()
      errs = dfjob['errfile'].to_list()
      num = len(jobc)
      #### submit jobs
      for i in range(num):
        print(jobc[i])
        # first clean
        if os.path.exists(outs[i]):
          os.remove(outs[i])
        if os.path.exists(errs[i]):
          os.remove(errs[i])
        # then submit
        if num > 1:
          p = os.system(jobc[i]+'&')
        else:
          p = os.system(jobc[i])
        if p != 0:
          raise ValueError('fails to submit job: ' + jobc[i])
      # print(iterc, istep, jobc, outs, errs,'\n')
      stat = np.zeros(num, np.int32)
      flag = 1
      while flag == 1:
        # refresh the disk
        for i in range(num):
          stat[i] = check_status_end(outs[i], flag_finished, flag_failed, nL)
        # update the status
        df.loc[dfjob.index, 'status'] = stat
        df.to_csv(jobNM, index=False)
        if -1 in stat:
          j = np.argmin(abs(stat+1))
          raise ValueError('Error happens: ' + errs[j])
        elif 0 in stat:
          j = np.argmin(abs(stat))
          name = outs[0].split('/')[1].split('.')[1]
          print("*** waiting %s ***"%name)
          time.sleep(time_sleep)
          flag = 1
        else:
          flag = 0
          print('****** job at (', iterc, istep, icg, ') is finished ****** \n', jobc)
