#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# make directory
os.makedirs('fig', exist_ok=True)

#Input paras
print('Usage:\n \t plot_dispersion.py optim obs/vmap_M0.bin\n')
fdir        =sys.argv[1] # optim
fnm       =sys.argv[2] # obs/vmap_M0.bin
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
fperiods    =optim['fperiods']
fstn        =optim['fstn']
vmin0       =2700
vmax0       =3800

# show figure or not
flag_show = 1
# figure name to save
fignm = 'fig/dispersion.png'
# figure size to save
figsize = [7, 3.5]
# figure resolution to save
figdpi  = 150

#periods
df = pd.read_csv(fperiods)
nper = df.shape[0]
periods = np.array(df['T0'])
# trace indexes
df = pd.read_csv(fstn)
ntraces = df.shape[0]
# data
mod = np.fromfile(fnm, np.float32).reshape([ntraces, nper])
xvel=np.arange(ntraces)+1
yvel=periods

print('\n')
print(' min and max ',mod.min(), mod.max())
print('axis 0:', xvel.min(),xvel.max())
print('axis 1:', yvel.min(),yvel.max())

# media show
fig, ax = plt.subplots(dpi=figdpi,figsize=(figsize[0],figsize[1]))
im = ax.imshow(mod.T, cmap='jet', aspect='auto', vmin=vmin0, vmax=vmax0,extent=[xvel[0],xvel[-1],yvel[-1],yvel[0]])
ax.set_title('Dispersion', fontsize=12)
plt.xlabel('Trace No.', fontsize=10)
plt.ylabel('Period (s)', fontsize=10)

#legend
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
clb = plt.colorbar(im, cax=cax)
clb.ax.set_ylabel('$\mathrm{\mathsf{(m/s)}}$')

plt.savefig(fignm)
if flag_show:
  plt.show()
