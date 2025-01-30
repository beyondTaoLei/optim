#!/usr/bin/env python3
import os,sys
import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# make directory
os.makedirs('fig', exist_ok=True)

#Input paras
print('Usage:\n \t plot_raypath.py optim fdio/sz001/raypath.dat\n')
fdir        =sys.argv[1] # optim
fin         =sys.argv[2] # os.path.join('fdio', evt, fraypath)
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
nx, ny      =optim['nxyfd']
ox, oy      =optim['oxyfd']
dx, dy      =optim['dxyfd']
fmodfd      =optim['fmodfd']
vmin0       =1500
vmax0       =4500

# show figure or not
flag_show    = 1
# figure name to save
fignm = 'fig/raypath.png'
# figure size to save
figsize = [6,3]
# figure resolution to save
figdpi  = 150
##########################################
mod = np.fromfile(fmodfd, np.float32).reshape([nx, ny])
xvel=np.arange(nx)*dx+ox
yvel=np.arange(ny)*dy+oy

#mod = np.fromfile('fdio/sz040/timefield.bin', np.float32).reshape([nx, ny])

raypaths_list = []
npoints_arr = np.loadtxt(fin, max_rows=1)
raypaths = np.loadtxt(fin, skiprows=1)
sum1 = np.cumsum(npoints_arr)
idx = np.hstack([[0],sum1]).astype(np.int32)
for i in range(len(idx)-1):
    raypaths_list.append(raypaths[idx[i]:idx[i+1]])

#plt.show()
fig, ax = plt.subplots(dpi=figdpi,figsize=(figsize[0],figsize[1]))
im = ax.imshow(mod.T, cmap='jet', aspect='auto', vmin=vmin0, vmax=vmax0,extent=[xvel[0],xvel[-1],yvel[-1],yvel[0]])
#im = ax.imshow(mod, interpolation='bilinear', cmap=cmap, aspect='auto', vmax=5, vmin=0,extent=[xvel[0],xvel[-1],yvel[-1],yvel[0]])
ax.set_title('Raypath', fontsize=17)
plt.xlabel('X', fontsize=15)
plt.ylabel('Y', fontsize=15)
for i in range(len(npoints_arr)):
    plt.plot(raypaths_list[i][:,0], raypaths_list[i][:,1],'k-',linewidth=2)

#legend
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
clb = plt.colorbar(im, cax=cax)
clb.ax.set_ylabel('$\mathrm{\mathsf{(m/s)}}$')

plt.savefig(fignm)
if flag_show:
  plt.show()
