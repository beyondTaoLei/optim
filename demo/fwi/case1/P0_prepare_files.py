#!/usr/bin/env python3
"""
Usage:
"""
import os
import sys
import json
import h5py
import numpy as np
import pandas as pd
from shutil import copyfile
import matplotlib.pylab as plt
from matplotlib.path import Path
from scipy.signal.windows import tukey
from scipy.ndimage import gaussian_filter

def ttaperL(t, t1, nalpha=20):
  """
  select one phase
  """
  nt = len(t)
  taper = np.zeros_like(t, dtype = np.float32)
  it1 = np.argmin(abs(t-t1))
  a = tukey(2*nalpha, 1.0)
  taper[it1:it1+nalpha] = a[:nalpha]
  taper[it1+nalpha:] = 1.0
  return taper

os.makedirs('list', exist_ok=True)
os.makedirs('source', exist_ok=True)

fdir                  =sys.argv[1]
foptim                =os.path.join(fdir,'optim.json')
optim                 =eval(open(foptim).read())
nxzfd                 =optim['nxzfd']
oxzfd                 =optim['oxzfd']
dxzfd                 =optim['dxzfd']
fwiroot               =optim['fwiroot']
fdio                  =optim['fdio']
delta                 =optim['delta']
nt                    =optim['nt']

sys.path.insert(0, os.path.join(fwiroot, 'prob'))
from util import write_su
from wavelet import wavelet

################################
#           geometry           #
################################
fevn_info = 'list/eventinfo.csv'
frec = 'list/receiver.dat'
# fwavelet = 'source/sz001.su'

# events
delay = 0.2
fc = 15.0
amp = 1.00e+10
test_inc = 3
Sx = np.arange(80)*200.0 + 700.0
Sz = 20.0 * np.ones_like(Sx) #depth
nsrc = len(Sx)
test = np.zeros(nsrc, np.int32)
test[::test_inc]=1
name=["sz%03d"%(i+1) for i in range(nsrc)]
df=pd.DataFrame(data={'station':name})
df["x"]       =Sx
df["z"]       =Sz
df['tshift']  =delay
df['fc']      =fc
df['amp']     =amp
df['test']    =test
df.to_csv(fevn_info,float_format="%.2f",index=False)

# wavelet
shape       =1 #New Ricker Wavelet, equal to SOFI2D
df = pd.read_csv(fevn_info)
name    = list(df['station'])
X       = list(df['x'])
Z       = list(df['z'])
tshift  = list(df['tshift'])
fc      = list(df['fc'])
amp0    = list(df['amp'])
nsrc    = len(name)

# generate the wavelets, and save to su-format file
for i in range(nsrc):
  print('wavelet: %4d of %4d'%(i+1, nsrc))
  scalel = 1000.0
  scalco = 1000.0
  itr = 1
  wav = wavelet(shape, nt, delta, tshift[i], fc[i], amp0[i])
  fnm = os.path.join('source', name[i]+'.su')
  data = wav.reshape([1, nt])
  Sxz = np.array([X[i], Z[i]])
  Rcoord = np.array([[X[i], Z[i]]])
  write_su(fnm, data, Sxz, Rcoord, delta)

time = np.arange(nt) * delta
plt.plot(time, 0.0*time)
plt.plot(time, wav)

plt.show()

# receiver
os.makedirs(os.path.join(fdio, 'comm'), exist_ok=True)
for evn in name:
  os.makedirs(os.path.join(fdio, evn), exist_ok=True)

Rx = np.arange(801)*20.0 + 500.0
Rz = 20.0 * np.ones_like(Rx) #depth
nrec = len(Rx)
fp = open(frec,'w')
for i in range(nrec):
    tmp=fp.write('%16.2f %10.2f\n'%(Rx[i], Rz[i]))
fp.close()

for i in range(nsrc):
  evn = name[i]
  fou = os.path.join(fdio, evn, 'receiver.dat')
  print("copy file: from %30s to %30s"%(frec, fou))
  copyfile(frec,  fou)

plt.plot(Rx, -Rz,'bo')
plt.plot(Sx, -Sz+1.0,'r*')
plt.show()

################################
#           taper              #
################################
nx, nz = nxzfd
dx, dz = dxzfd
ox, oz = oxzfd
grad_taper_file = 'list/taper.bin'
tdepth = 43 * dz
x = np.arange(nx) * dx + ox
z = np.arange(nz) * dz + oz
taper1d = ttaperL(z, tdepth, 1)

fig1 = plt.figure(dpi=300,figsize=(1,2), constrained_layout=True)
plt.plot(taper1d, -z/1.e3)
plt.title('taper')
plt.xlabel('coeff.')
plt.ylabel('Z (km)')
plt.grid(True)
plt.show()

# order with FD coordinate formats
# taper1d = taper1d[::-1]
taper = np.tile(taper1d,(nx,1)).T 
taper.T.tofile(grad_taper_file) #ximage show
print(taper.shape)

print("Finished...", __file__)