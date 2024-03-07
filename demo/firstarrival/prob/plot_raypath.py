#!/usr/bin/env python3
import os,sys
import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm

# make directory
os.makedirs('fig', exist_ok=True)

#Input paras
fdir        =sys.argv[1] # optim
fin         =sys.argv[2] # os.path.join('fdio', evt, fraypath)
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
ox          =optim['o2g']
oy          =optim['o1g']
nx          =optim['n2g']
ny          =optim['n1g']
dx          =optim['d2g']
dy          =optim['d1g']
fmodfd      =optim['fmodfd']

vmin0=1500
vmax0=4500
xvel=np.arange(nx)*dx+ox
yvel=np.arange(ny)*dy+oy

mod = np.fromfile(fmodfd, np.float32).reshape([nx, ny])
#mod = np.fromfile('fdio/sz040/timefield.bin', np.float32).reshape([nx, ny])

raypaths_list = []
npoints_arr = np.loadtxt(fin, max_rows=1)
raypaths = np.loadtxt(fin, skiprows=1)
sum1 = np.cumsum(npoints_arr)
idx = np.hstack([[0],sum1]).astype(np.int32)
for i in range(len(idx)-1):
    raypaths_list.append(raypaths[idx[i]:idx[i+1]])

#plt.show()
cmap=cm.bwr
fig, ax = plt.subplots(figsize=(12,10))
mod=mod.transpose()#for ximage
im = ax.imshow(mod, interpolation='bilinear', cmap=cmap, aspect='auto', vmin=vmin0, vmax=vmax0,extent=[xvel[0],xvel[-1],yvel[-1],yvel[0]])
#im = ax.imshow(mod, interpolation='bilinear', cmap=cmap, aspect='auto', vmax=5, vmin=0,extent=[xvel[0],xvel[-1],yvel[-1],yvel[0]])
ax.set_title('Raytph', fontsize=17)
plt.xlabel('X', fontsize=15)
plt.ylabel('Y', fontsize=15)
for i in range(len(npoints_arr)):
    plt.plot(raypaths_list[i][:,0], raypaths_list[i][:,1],'k-',linewidth=2)
handles, labels = plt.gca().get_legend_handles_labels()
fig.colorbar(im)
fig.savefig(os.path.join('fig','raypath.eps'),format='eps')
plt.show()
print("Finished...", __file__)
