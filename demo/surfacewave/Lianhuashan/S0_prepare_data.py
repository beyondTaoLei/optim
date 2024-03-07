#!/usr/bin/env python3
import os
import json
import scipy.io
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from scipy.ndimage import gaussian_filter, median_filter

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
n1, n2, n3 = 100, 1400, 1
o1, o2, o3 = 0.0, 0.0, 0.0
d1, d2, d3 = 100.0, 100.0, 100.0
fcmp= 'list/cmp.csv'
fperiods = 'list/periods.csv'
fobs = 'obs/vmap_M0.bin'
finit = 'model/init/init.vs'

fstn = 'list/stationinfo.csv'
df = pd.read_csv(fstn)
stns = list(df['station'])


Y = df['X'].to_numpy() # meter
X = np.zeros_like(Y)
Z = np.zeros_like(Y)
d = {'X':list(X)}
df = pd.DataFrame(data=d)
df['Y'] = Y
df['Z'] = Z
I3s = ((X-o3)/d3).round().astype(np.int32)
I2s = ((Y-o2)/d2).round().astype(np.int32)
I1s = ((Z-o1)/d1).round().astype(np.int32)
df['offset'] = I2s + I3s*n2
df.to_csv(fcmp,float_format="%.4f",index=False)

# periods
T = np.arange(0.5, 5.1, 0.2)
#T = np.arange(2, 31, 1)*0.1
T0 = ['%.2f'%i for i in T]
phase = ['ph_%d'%i for i in range(len(T))]
d={'phase':phase}
df=pd.DataFrame(data=d)
df["T0"] =T0
df.to_csv(fperiods,index=False)

# observed data
nper = len(T)
ntr = 112
phvs = np.zeros([ntr, nper], np.float32)
for itr in range(ntr):
    fnm = os.path.join('/home/tao/seis_software/DATA/Lianhuashan/curve', stns[itr]+'curve.txt')
    print(stns[itr]+'curve.txt')
    a = np.loadtxt(fnm)[::-1,:]
    a[:, 0] = 1.0 / a[:, 0]
    phv = np.interp(T, a[:,0], a[:,1])
    phvs[itr,:] = phv # m/s

phvs.tofile(fobs[:-4]+'_ori.bin')
sigmas=np.array([4,2])/np.sqrt(8.0)
phvs_old = phvs.copy()
print(phvs_old[:,-1])
phvs = gaussian_filter(phvs, sigmas, mode='nearest')
phvs=np.where(phvs>0.0, phvs, phvs_old)
phvs.tofile(fobs)
print(phvs[:,-1])

#initial models
dep = np.arange(n1)*d1 + o1
vs2d = np.zeros([ntr, n1], np.float32)
for itr in range(ntr):
    phv = phvs[itr,:]
    depth = 1.0*T[phv>0]*phv[phv>0]/3.0
    vs = 1.1*phv[phv>0]
    vs2d[itr, :] = np.interp(dep, depth, vs)

# ntr to n2
Y2 = np.arange(n2)*d2 + o2
vs2d2 = np.zeros([n2, n1], np.float32)
print(Y)
print(Y2)
#plt.plot(Y, np.zeros(ntr),'.')
#plt.plot(Y2, np.zeros(n2), '.')
#plt.show()
for i in range(n1):
    vs2d2[:, i] = np.interp(Y2, Y[::-1], vs2d[::-1, i])

sigmas=np.array([10,8])/np.sqrt(8.0)
vs2d2 = gaussian_filter(vs2d2, sigmas, mode='nearest')#sigma

vp=vs2vp(vs2d2.flatten())
rho=vp2rho(vp)
vs2d2.tofile(finit)
vp.tofile(finit.replace('.vs', '.vp'))
rho.tofile(finit.replace('.vs', '.rho'))
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

