#!/usr/bin/env python3
"""
Usage:
    python *.py
"""
import os
import sys
import numpy as np
import pandas as pd
from util import taper2d

#Input paras
fdir              =sys.argv[1]
foptim            =os.path.join(fdir,'optim.json')
optim             =eval(open(foptim).read())
nx, nz            =optim['nxzfd']
ox, oz            =optim['oxzfd']
dx, dz            =optim['dxzfd']
fdio              =optim['fdio']
fevn_info         =optim['fevn_info']
grad_taper_shot   =optim['grad_taper_shot']
srtradius1        =optim['srtradius1']
srtradius2        =optim['srtradius2']
fgrad             =optim['fgrad']
fout              =fgrad+'_0'

# interpretation
x = np.arange(nx) * dx + ox
z = np.arange(nz) * dz + oz
num = nx * nz

# read source info.
df = pd.read_csv(fevn_info)
events  = list(df['station'])
Sx      = df['x'].to_numpy()
Sz      = df['z'].to_numpy()
nsrc    = len(events)

# sum up
sum1 = np.zeros(num, np.float32)
for isrc in range(nsrc):
  kevnm = events[isrc]
  fnm = os.path.join(fdio, kevnm, 'gradvp_shot1.bin')
  if not os.path.exists(fnm):
    raise ValueError("no file: " + fnm)
  grad = np.fromfile(fnm, np.float32)
  # tapering
  if grad_taper_shot == 1:
    fnm2 = os.path.join(fdio, kevnm, 'receiver.dat')
    Rcoord = np.loadtxt(fnm2)
    ntr, tmp = Rcoord.shape
    posx = np.hstack([Sx[isrc], Rcoord[:, 0]])
    posz = np.hstack([Sz[isrc], Rcoord[:, 1]])
    iscoord = np.zeros([ntr+1, 2], np.int32)
    for i in range(ntr+1):
      iscoord[i, 0] = np.argmin(abs(posx[i] - x))
      iscoord[i, 1] = np.argmin(abs(posz[i] - z))
    iscoord2 = np.unique(iscoord, axis=0)
    taper = taper2d(nx, nz, iscoord2, srtradius1, srtradius2) # nx, nz
    grad *= taper.flatten()
  # sum up
  sum1 += grad

sum1.tofile(fout)
max1=np.abs(sum1).max()/10.

print("\nximage  n1=%d cmap=hsv2 bclip=%f wclip=%f <%s\n"% \
     (nz, max1, -max1, fout))
print("*"*60+"\n*Output file: "+fout)
print("Finished...", __file__)