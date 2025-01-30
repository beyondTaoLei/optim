#/home/tao/seis_software
import os
import sys
import json
import shutil
import numpy as np
import pandas as pd
from util import linefitting

#Input paras
fdir        =sys.argv[1]
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
niter_min   =optim['niter_min']
pro         =optim['pro']
pro2        =optim['pro2']
pro3        =optim['pro3']
iter0       =optim['iter0']
iterc       =optim['iterc']
fiter       ='iter'+str(iterc)

# read current information
fcost = os.path.join(fdir, fiter, 'fcost.dat')
misfit0 = np.loadtxt(fcost).item()
# write 
fhist = os.path.join(fdir, 'history.csv')
if not os.path.exists(fhist):
  data = [[iterc, misfit0], ]
  df = pd.DataFrame(data, columns=['iter', 'misfit'])
else:
  df = pd.read_csv(fhist)
  df = df[df['iter'] < iterc]
  df.loc[len(df.index)] = [iterc, misfit0]
  df['iter'] = df['iter'].astype(np.int32)

df.to_csv(fhist, index=False)

# check stopping
check_len = 4
flag = 0
if iterc >= iter0 + niter_min - 1:
  # check stopping 1
  misfit2 = df[df['iter']==iterc-2]['misfit'].item()
  diff = abs((misfit2 - misfit0)/misfit2)
  print('misfit2, misfit0, diff:', misfit2, misfit0, diff)
  if diff <= pro:
    flag = 1
  else:
    # check stopping 2
    if pro2 > 0.0 or pro3 > 0.0:
      df1 = df.loc[(df['iter'] >= iter0) & (df['iter'] < iter0+check_len)]
      df2 = df.loc[(df['iter'] > iterc-check_len) & (df['iter'] <= iterc)]
      Y1 = df1['misfit'].to_numpy()
      Y2 = df2['misfit'].to_numpy()
      Xs = np.arange(check_len) + 1
      a1, b1, ave1 = linefitting(Xs, Y1)
      a2, b2, ave2 = linefitting(Xs, Y2)
      print('a1, b1, ave1:', a1, b1, ave1)
      print('a2, b2, ave2:', a2, b2, ave2)
      if(a1 < 0.0 and a2 < 0.0):
        if(a2 / a1 < pro2 or ave2 / ave1 < pro3):
          print('inversion stop code: 1')
          flag = 1
      elif(a1 < 0.0 and a2 > 0.0):
        if(abs(a2 / a1) < pro2 and ave2 / ave1 < pro3):
          print('inversion stop code: 2')
          flag = 1
      elif(a1 < 0.0 and a2 == 0.0):
        print('inversion stop code: 3')
        flag = 1
      elif(a1 > 0.0):
        if(ave2 / ave1 < pro3):
          print('inversion stop code: 4')
          flag = 1

if flag == 1:
  print('Inversion will goes to next stage!')
  optim.update({ 'abort' : 1})
  fout = open(foptim,'w')
  fout.write(json.dumps(optim,indent=4))
  fout.close()
else:
  print('inversion goes on ...')

print("Finished...", __file__)        