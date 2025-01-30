#!/usr/bin/env python3
"""
the inversion example for Rosenbrock problem
Usage:
    
"""
import os
import sys
import numpy as np
import pandas as pd

fdir        =sys.argv[1]
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
nx, ny, nz  =optim['nxyzfd']
ox, oy, oz  =optim['oxyzfd']
dx, dy, dz  =optim['dxyzfd']

#output files
os.makedirs('list', exist_ok=True)
fstn= 'list/stn.csv'
fperiods = 'list/periods.csv'

#Y = np.arange(ny)*dy + oy
Y = np.arange(10, ny, 10)*dy + oy
X = np.zeros_like(Y)
Z = np.zeros_like(Y)

df = pd.DataFrame()
df['station'] = ['sz%03d'%i for i in range(len(X))]
df['x'] = X
df['y'] = Y
df['z'] = Z
# I3s = ((X-ox)/dx).round().astype(np.int32)
# I2s = ((Y-oy)/dy).round().astype(np.int32)
# I1s = ((Z-oz)/dz).round().astype(np.int32)
# df['offset'] = I2s + I3s*ny
df['ix'] = ((X-ox)/dx).round().astype(np.int32)
df['iy'] = ((Y-oy)/dy).round().astype(np.int32)
df['iz'] = ((Z-oz)/dz).round().astype(np.int32)
df['offset'] = df['iy'] + df['ix']*ny
df.to_csv(fstn,float_format="%.4f",index=False)

#periods
T = np.arange(3, 60,2)
T0 = ['%.2f'%i for i in T]
phase = ['ph_%d'%i for i in range(len(T))]
d={'phase':phase}
df=pd.DataFrame(data=d)
df["T0"] =T0
df.to_csv(fperiods,index=False)

    