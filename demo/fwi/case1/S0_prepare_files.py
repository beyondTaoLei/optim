#!/usr/bin/env python3
"""
the inversion example for Rosenbrock problem
Usage:
    
"""
import os
import numpy as np
import pandas as pd

os.makedirs('list', exist_ok=True)
fsrc = 'list/source.dat'
fevn_info = 'list/eventinfo.csv'
frec = 'list/receiver.dat'

# source 
X = np.arange(700, 16500+1, 200)*1.0
num = len(X)
Y = np.ones(num) * 20.0
delay = 0.2
fc = 15.0
amp = 1.00e+10
test_inc = 3
fp = open(fsrc,'w')
tmp=fp.write('%d\n'% num)
for i in range(num):
    tmp=fp.write('%16.2f %5.2f %10.2f %10.2f %10.2f %10.2e\n'%(X[i], 0.00, Y[i], delay, fc, amp))
fp.close()

test = np.zeros(num, np.int32)
test[::test_inc]=1
name = ['sz%03d'%(i+1) for i in range(num)]
d = {'station': name}
df = pd.DataFrame(data=d)
df['x'] = X
df['y'] = Y
df['tshift'] = delay
df['fc'] = fc
df['amp'] = amp
df['test'] = test
df.to_csv(fevn_info,float_format="%.2f",index=False)

# receiver
X = np.arange(500, 16500+1, 20)*1.0
num = len(X)
Y = np.ones(num) * 20.0
fp = open(frec,'w')
for i in range(num):
    tmp=fp.write('%16.2f %10.2f\n'%(X[i], Y[i]))
fp.close()

    