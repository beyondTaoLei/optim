#!/usr/bin/env python3
import os
import sys
import h5py
import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

##########################################
print('Usage:\n \t plot_mod.py optim model/inv/inv.vp \n')

fdir        =sys.argv[1] # optim
fnm         =sys.argv[2] # model/true/true.vp
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
nx, ny      =optim['nxyfd']
ox, oy      =optim['oxyfd']
dx, dy      =optim['dxyfd']
vmin0       =1500.0
vmax0       =4500.0

# show figure or not
flag_show    = 1
# figure name to save
fignm = fnm[:-3]+'.png'
# figure size to save
figsize = [6,3]
# figure resolution to save
figdpi  = 150
##########################################
mod    = np.fromfile(fnm, np.float32).reshape([nx, ny])
xvel = np.arange(nx) * dx + ox
yvel = np.arange(ny) * dy + oy

# mod = mod/1e3
print('\n')
print(' min and max ',mod.min(), mod.max())
print('axis 0:', xvel.min(),xvel.max())
print('axis 1:', yvel.min(),yvel.max())

# media show
fig, ax = plt.subplots(dpi=figdpi,figsize=(figsize[0],figsize[1]))
im = ax.imshow(mod.T, cmap='jet', aspect='auto', vmin=vmin0, vmax=vmax0,extent=[xvel[0],xvel[-1],yvel[-1],yvel[0]])
ax.set_title('Vp model', fontsize=17)
plt.xlabel('X', fontsize=15)
plt.ylabel('Y', fontsize=15)

#legend
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
clb = plt.colorbar(im, cax=cax)
clb.ax.set_ylabel('$\mathrm{\mathsf{(m/s)}}$')

plt.savefig(fignm)
if flag_show:
  plt.show()
