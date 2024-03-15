#!/usr/bin/env python3
"""
the inversion example for Rosenbrock problem
Usage:
    
"""
import os
import numpy as np
import pandas as pd

os.makedirs('list', exist_ok=True)
#Z, Y, X: (1, 2, 3)
n1, n2, n3 = 400, 1300, 1
o1, o2, o3 = 0.0, 0.0, 0.0
d1, d2, d3 = 500.0, 500.0, 500.0
fcmp= 'list/cmp.csv'
fperiods = 'list/periods.csv'

#Y = np.arange(n2)*d2 + o2
Y = np.arange(10, n2, 10)*d2 + o2
X = np.zeros_like(Y)
Z = np.zeros_like(Y)

df = pd.DataFrame()
df['station'] = ['sz%03d'%i for i in range(len(X))]
df['x'] = X
df['y'] = Y
df['z'] = Z
I3s = ((X-o3)/d3).round().astype(np.int32)
I2s = ((Y-o2)/d2).round().astype(np.int32)
I1s = ((Z-o1)/d1).round().astype(np.int32)
df['offset'] = I2s + I3s*n2
df.to_csv(fcmp,float_format="%.4f",index=False)

#periods
T = np.arange(3, 60,2)
T0 = ['%.2f'%i for i in T]
phase = ['ph_%d'%i for i in range(len(T))]
d={'phase':phase}
df=pd.DataFrame(data=d)
df["T0"] =T0
df.to_csv(fperiods,index=False)

    