#!/usr/bin/env python3
"""
the inversion example for Rosenbrock problem
Usage:
    cp model/true/true.vp model/inv/inv.vp
    cp model/true/true.vs model/inv/inv.vs
    cp model/true/true.rho model/inv/inv.rho
    
    python /home/tao/Nutstore_Files/works/seisplatform/mod_OPTIM/demo/surfacewave_phv_lsqr/S1_generate_data.py optim
"""
import os
import math
import numpy as np
import pandas as pd
import matplotlib.pylab as plt

def remove_duplicates(offsets):
    #https://stackoverflow.com/questions/30003068/how-to-get-a-list-of-all-indices-of-repeated-elements-in-a-numpy-array
    # creates an array of indices, sorted by unique element
    idx_sort = np.argsort(offsets)
    # sorts records array so all unique elements are together 
    sorted_records_array = offsets[idx_sort]
    # returns the unique values, the index of the first occurrence of a value, and the count for each element
    vals, idx_start, count = np.unique(sorted_records_array, return_counts=True, return_index=True)
    # splits the indices into separate arrays
    res = np.split(idx_sort, idx_start[1:])
    idxs = [ii[0] for ii in res]
    return idxs

os.makedirs('list', exist_ok=True)
os.makedirs('obs', exist_ok=True)

# read coordinates
sxy0 = np.loadtxt('/home/tao/seis_software/DATA/LongBeach/origin/sxy.txt')
xmin, ymin = sxy0[:,0].min(), sxy0[:,1].min()
sxy = 0.0 * sxy0
sxy[:, 0] = sxy0[:, 0] - xmin 
sxy[:, 1] = sxy0[:, 1] - ymin
ntraces = len(sxy)
print('x(min, max):', sxy[:,0].min(), sxy[:,0].max())
print('y(min, max):', sxy[:,1].min(), sxy[:,1].max())
#Z, Y, X: (1, 2, 3)
n1, n2, n3 = 100, 340, 240
o1, o2, o3 = 0.0, 0.0, 0.0
d1, d2, d3 = 10.0, 100.0, 100.0
fcmp= 'list/cmp.csv'
fperiods = 'list/periods.csv'

offsets = np.zeros(ntraces, np.int32)
for i in range(ntraces):
    xo = math.floor((sxy[i,0] - o3) / d3) * d3 + o3
    yo = math.floor((sxy[i,1] - o2) / d2) * d2 + o2
    xs = np.array([xo, xo + d3, xo + d3, xo])
    ys = np.array([yo, yo, yo + d2, yo + d2])
    dist = (xs-sxy[i,0])**2 + (ys-sxy[i,1])**2
    idx = dist.argmin()
    I3 = int((xs[idx]- o3) / d3)
    I2 = int((ys[idx]- o2) / d2)
    offsets[i] = I2 + I3*n2

#remove duplicates
idx_nm = np.array([i+1 for i in range(ntraces)], np.int32)
idxs = remove_duplicates(offsets)
print('ntraces(old. new):', ntraces, len(idxs))
ntraces = len(idxs)
df = pd.DataFrame()
df['station'] = ['sz%03d'%i for i in range(ntraces)]
df['x0'] = sxy0[:, 0][idxs]
df['y0'] = sxy0[:, 1][idxs]
df['x'] = sxy[:, 0][idxs]
df['y'] = sxy[:, 1][idxs]
df['z'] = 0.0
df['offset'] = offsets[idxs]
idx_nm = idx_nm[idxs]
df['name'] = idx_nm
df.to_csv(fcmp,float_format="%.4f",index=False) #order with offset


#period from short to long
pre = '/home/tao/seis_software/DATA/LongBeach_DC/data'
min1, max1 = 10, -10
for i in range(ntraces):
    fnm = os.path.join(pre, str(idx_nm[i])+'.txt')
    a = np.loadtxt(fnm)
    min1 = min(a[:,0].min(), min1)
    max1 = max(a[:,0].max(), max1)
    #plt.plot(a[:,0],'.')
    
#plt.show()    
print('freq. min, max:', min1, max1) # 0.6 3.8
print(1/min1, 1/max1)
T = np.arange(27, 167, 5) * 0.01
nper = len(T)
print('T:', T)
print('F:', 1/T)
T0 = ['%.2f'%i for i in T]
phase = ['ph_%d'%i for i in range(nper)]
d={'phase':phase}
df=pd.DataFrame(data=d)
df["T0"] =T0
df.to_csv(fperiods,index=False)

#phase map 
print('ntraces, nper:', ntraces, nper)
phvmap0 = np.zeros([ntraces, nper], np.float32)#m/s
phvmap1 = np.zeros([ntraces, nper], np.float32)#m/s
for i in range(ntraces):
    fnm = os.path.join(pre, str(idx_nm[i])+'.txt')
    a = np.loadtxt(fnm)
    phv0 = a[a[:,2]==0.0][::-1,:]
    phv1 = a[a[:,2]==1.0][::-1,:]
    if len(phv0)>0:
        print(fnm)
        phvmap0[i,:] = np.interp(T, 1.0/phv0[:,0], phv0[:,1], -1.0, -1.0)
        if len(phv1)>0:
            phvmap1[i,:] = np.interp(T, 1.0/phv1[:,0], phv1[:,1], -1.0, -1.0)
    else:
        print('No fund.:', fnm)

#save phase velocity map
phvmap0[:] *= 1000.0 #m/s
phvmap1[:] *= 1000.0 #m/s
phvmap0.tofile('obs/vmap_M0.bin')
phvmap1.tofile('obs/vmap_M1.bin')


    

    