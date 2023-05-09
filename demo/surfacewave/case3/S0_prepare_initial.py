#!/usr/bin/env python3
import os
import json
import scipy.io
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter

def vp2rho(vp):
    """
    Empirical Relations between Elastic Wavespeeds and Density in the Earth's Crust
    in equation (1)
    vp: input in m/s and valid for vp between [1.5 8.5]*1000.0 m/s
    """
    vpkm=vp/1000.0 #km/s
    vp2=vpkm*vpkm
    vp3=vpkm*vp2
    vp4=vp2*vp2
    vp5=vp2*vp3
    rho=1.6612*vpkm-0.4721*vp2+0.0671*vp3-0.0043*vp4+0.000106*vp5 #g/cm3
    return rho*1000.0 #

def vs2vp(vs):
    """
    Empirical Relations between Elastic Wavespeeds and Density in the Earth's Crust
    in equation (9)
    vs: input in m/s  and valid for vs between [0 4.5]*1000.0 m/s
    """
    vskm=vs/1000.0
    vs2=vskm*vskm
    vs3=vskm*vs2
    vs4=vs2*vs2
    vp=0.9409 + 2.0947*vskm-0.8206*vs2+ 0.2683*vs3 - 0.0251*vs4
    return vp*1000.0 #m/s

os.makedirs('list', exist_ok=True)
os.makedirs('obs', exist_ok=True)
os.makedirs('model', exist_ok=True)
os.makedirs('model/init', exist_ok=True)
os.makedirs('model/inv', exist_ok=True)

#Z, Y, X: (1, 2, 3)
n1, n2, n3 = 100, 340, 240
o1, o2, o3 = 0.0, 0.0, 0.0
d1, d2, d3 = 10.0, 100.0, 100.0
# output slices
snap1 = [0, 20, 40, 80, 99]
snap2 = [100, 200, 300]
snap3 = [100, 150, 200]
sigmas=np.array([8,8,6])/np.sqrt(8.0) # in grid
# input files
fcmp= 'list/cmp.csv'
fperiods = 'list/periods.csv'
fobs = 'obs/vmap_M0.bin'
# output files
finit = 'model/init/init.vs'
fsnap = 'model/init/snap'


# cmp info.
df = pd.read_csv(fcmp)
offsets = df['offset']
ntr = len(offsets)

# periods
df = pd.read_csv(fperiods)
T = df["T0"]
nper = len(T)

# observed data
phvs = np.fromfile(fobs, np.float32).reshape([ntr, nper])# m/s

#initial models
X = np.arange(n3)*d3 + o3
Y = np.arange(n2)*d2 + o2
Z = np.arange(n1)*d1 + o1
vs3d = np.zeros([ntr, n1], np.float32)
for itr in range(ntr):
    phv = phvs[itr,:]
    idxs = np.where(phv>0)[0]
    if len(idxs) > 0:
        phvtmp = phv[idxs]
        Ttmp = T[idxs]
        depth = 1.0*Ttmp*phvtmp/3.0
        vs = 1.1*phvtmp
        vs3d[itr, :] = np.interp(Z, depth, vs, -1.0, -1.0)
    else:
        print('No data:', itr, ' of ', ntr)

# ntr to n2
vs = np.zeros([n3, n2, n1], np.float32)
grid_3, grid_2 = np.mgrid[0:(n3-1):n3*1j, 0:(n2-1):n2*1j]
recordy, recordn = [], []
for i in range(n1):
    mod = vs3d[:, i]
    idxs = np.where(mod>0.0)[0]
    if len(idxs) > 0:
        print(i, ' of ', n1, '( ', len(idxs), n2*n3, ' )')
        offset_tmp = offsets[idxs]
        mod_tmp = mod[idxs]
        points=np.array([offset_tmp//n2, offset_tmp%n2]).T
        vs[:, :, i] = griddata(points, mod_tmp, (grid_3, grid_2), method='nearest')
        recordy.append(i)
    else:
        print(i, ' of ', n1, '( No data )')
        recordn.append(i)

print('No value in depth index:', recordn)
numno = len(recordn)
if numno > 0:
    Zy = np.array(recordy)*d1 + o1
    Zn = np.array(recordn)*d1 + o1
    for i in range(numno):
        dist = np.abs(Zn[i]-Zy)
        idx = dist.argmin()
        print(recordy[idx], ' ---> ', recordn[i])
        if Zn[i] > Zy[idx]:
            vs[:, :, recordn[i]] = 1.01 * vs[:, :, recordy[idx]]
        else:
            vs[:, :, recordn[i]] = 0.99 * vs[:, :, recordy[idx]]

# smooth
vs = gaussian_filter(vs, sigmas, mode='nearest')#sigma

vp=vs2vp(vs.flatten())
rho=vp2rho(vp)
vs.tofile(finit)
vp.tofile(finit.replace('.vs', '.vp'))
rho.tofile(finit.replace('.vs', '.rho'))

# snap in z dir.
for i in snap1: 
    vs[:,:,i].tofile(fsnap+'_z'+'%.2f'%(o1+i*d1)+'m.bin')

# snap in y dir.
for i in snap2: 
    vs[:,i,:].tofile(fsnap+'_y'+'%.2f'%(o2+i*d2)+'m.bin')

# snap in x dir.    
for i in snap3: 
    vs[i,:,:].tofile(fsnap+'_x'+'%.2f'%(o3+i*d3)+'m.bin')

paras ={
    'n3':n3,
    'n2':n2,
    'n1':n1,
    'd3':d3,
    'd2':d2,
    'd1':d1,
    'o3':o3,
    'o2':o2, 
    'o1':o1,
    'unit':'m,m/s,m/s,kg/m3'
    }
fout = open(finit.replace('.vs', '.json'),'w')
fout.write(json.dumps(paras,indent=4))
fout.close()

print("ntr, nper, n1, n2, n3:", ntr, nper, n1, n2, n3)

